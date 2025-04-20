#include "ns3/ipv4.h"
#include "ns3/packet.h"
#include "ns3/ipv4-header.h"
#include "ns3/pause-header.h"
#include "ns3/flow-id-tag.h"
#include "ns3/boolean.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "switch-node.h"
#include "qbb-net-device.h"
#include "ppp-header.h"
#include "ns3/int-header.h"
#include <cmath>

#include "cmsketch.h"
#include <typeinfo>
#include <cstdint>
#include <sstream>
#include <stdlib.h>
#include <stdexcept> // 包含异常处理

#define DEPTH 7
#define WIDTH 100000
#define CRUX true

namespace ns3 {

TypeId SwitchNode::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::SwitchNode")
    .SetParent<Node> ()
    .AddConstructor<SwitchNode> ()
	.AddAttribute("EcnEnabled",
			"Enable ECN marking.",
			BooleanValue(false),
			MakeBooleanAccessor(&SwitchNode::m_ecnEnabled),
			MakeBooleanChecker())
	.AddAttribute("CcMode",
			"CC mode.",
			UintegerValue(0),
			MakeUintegerAccessor(&SwitchNode::m_ccMode),
			MakeUintegerChecker<uint32_t>())
	.AddAttribute("AckHighPrio",
			"Set high priority for ACK/NACK or not",
			UintegerValue(0),
			MakeUintegerAccessor(&SwitchNode::m_ackHighPrio),
			MakeUintegerChecker<uint32_t>())
	.AddAttribute("MaxRtt",
			"Max Rtt of the network",
			UintegerValue(9000),
			MakeUintegerAccessor(&SwitchNode::m_maxRtt),
			MakeUintegerChecker<uint32_t>())
  ;
  return tid;
}

SwitchNode::SwitchNode(){
	
	/***********************
	 * DNN dcqcn Scheduler
	 ***********************/
	cms_PktSize.cms_init(&cms_PktSize, WIDTH, DEPTH);
	cms_TimeStamp.cms_init(&cms_TimeStamp,WIDTH, DEPTH);
	cms_State.cms_init(&cms_State,WIDTH, DEPTH);
	fflush(stdout);
	fair_factor=7;
	//run_flag = 0;
	/***********************
	 * DNN dcqcn Scheduler
	 ***********************/
	m_ecmpSeed = m_id;
	m_node_type = 1;
	m_mmu = CreateObject<SwitchMmu>();
	for (uint32_t i = 0; i < pCnt; i++)
		for (uint32_t j = 0; j < pCnt; j++)
			for (uint32_t k = 0; k < qCnt; k++)
				m_bytes[i][j][k] = 0;
	for (uint32_t i = 0; i < pCnt; i++)
		m_txBytes[i] = 0;
	for (uint32_t i = 0; i < pCnt; i++)
		m_lastPktSize[i] = m_lastPktTs[i] = 0;
	for (uint32_t i = 0; i < pCnt; i++)
		m_u[i] = 0;
	/*******DNN dcqcn scheduler */
	for (uint32_t i = 0; i < pCnt; i++)
		for (uint32_t k = 0; k < qCnt; k++){
			run_flag[i][k] = 0;
			tag[i][k] = 0;
	}
	/*******DNN dcqcn scheduler */
}

SwitchNode::~SwitchNode(){
}

int SwitchNode::GetOutDev(Ptr<const Packet> p, CustomHeader &ch){
	// look up entries
	auto entry = m_rtTable.find(ch.dip);

	// no matching entry
	if (entry == m_rtTable.end())
		return -1;

	// entry found
	auto &nexthops = entry->second;

	// pick one next hop based on hash
	union {
		uint8_t u8[4+4+2+2];
		uint32_t u32[3];
	} buf;
	buf.u32[0] = ch.sip;
	buf.u32[1] = ch.dip;
	if (ch.l3Prot == 0x6)
		buf.u32[2] = ch.tcp.sport | ((uint32_t)ch.tcp.dport << 16);
	else if (ch.l3Prot == 0x11)
		buf.u32[2] = ch.udp.sport | ((uint32_t)ch.udp.dport << 16);
	else if (ch.l3Prot == 0xFC || ch.l3Prot == 0xFD)
		buf.u32[2] = ch.ack.sport | ((uint32_t)ch.ack.dport << 16);

	uint32_t idx = EcmpHash(buf.u8, 12, m_ecmpSeed) % nexthops.size();
	return nexthops[idx];
}

void SwitchNode::CheckAndSendPfc(uint32_t inDev, uint32_t qIndex){
	Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
	if (m_mmu->CheckShouldPause(inDev, qIndex)){
		device->SendPfc(qIndex, 0);
		m_mmu->SetPause(inDev, qIndex);
	}
}
void SwitchNode::CheckAndSendResume(uint32_t inDev, uint32_t qIndex){
	Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
	if (m_mmu->CheckShouldResume(inDev, qIndex)){
		device->SendPfc(qIndex, 1);
		m_mmu->SetResume(inDev, qIndex);
	}
}

void SwitchNode::SendToDev(Ptr<Packet>p, CustomHeader &ch){
	int idx = GetOutDev(p, ch);
	if (idx >= 0){
		NS_ASSERT_MSG(m_devices[idx]->IsLinkUp(), "The routing table look up should return link that is up");

		// determine the qIndex
		uint32_t qIndex;
		if (ch.l3Prot == 0xFF || ch.l3Prot == 0xFE || (m_ackHighPrio && (ch.l3Prot == 0xFD || ch.l3Prot == 0xFC))){  //QCN or PFC or NACK, go highest priority
			qIndex = 0;
		}else{
			qIndex = (ch.l3Prot == 0x06 ? 1 : ch.udp.pg); // if TCP, put to queue 1
			//std::cout<<"[debug] ch.udp.pg"<<(ch.udp.pg)<<std::endl;
		}

		// admission control
		FlowIdTag t;
		p->PeekPacketTag(t);
		uint32_t inDev = t.GetFlowId();
		if (qIndex != 0){ //not highest priority
			if (m_mmu->CheckIngressAdmission(inDev, qIndex, p->GetSize()) && m_mmu->CheckEgressAdmission(idx, qIndex, p->GetSize())){			// Admission control
				m_mmu->UpdateIngressAdmission(inDev, qIndex, p->GetSize());
				m_mmu->UpdateEgressAdmission(idx, qIndex, p->GetSize());
			}else{
				return; // Drop
			}
			CheckAndSendPfc(inDev, qIndex);
		}
		m_bytes[inDev][idx][qIndex] += p->GetSize();
		m_devices[idx]->SwitchSend(qIndex, p, ch);
	}else
		return; // Drop
}

uint32_t SwitchNode::EcmpHash(const uint8_t* key, size_t len, uint32_t seed) {
  uint32_t h = seed;
  if (len > 3) {
    const uint32_t* key_x4 = (const uint32_t*) key;
    size_t i = len >> 2;
    do {
      uint32_t k = *key_x4++;
      k *= 0xcc9e2d51;
      k = (k << 15) | (k >> 17);
      k *= 0x1b873593;
      h ^= k;
      h = (h << 13) | (h >> 19);
      h += (h << 2) + 0xe6546b64;
    } while (--i);
    key = (const uint8_t*) key_x4;
  }
  if (len & 3) {
    size_t i = len & 3;
    uint32_t k = 0;
    key = &key[i - 1];
    do {
      k <<= 8;
      k |= *key--;
    } while (--i);
    k *= 0xcc9e2d51;
    k = (k << 15) | (k >> 17);
    k *= 0x1b873593;
    h ^= k;
  }
  h ^= len;
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

void SwitchNode::SetEcmpSeed(uint32_t seed){
	m_ecmpSeed = seed;
}

void SwitchNode::AddTableEntry(Ipv4Address &dstAddr, uint32_t intf_idx){
	uint32_t dip = dstAddr.Get();
	m_rtTable[dip].push_back(intf_idx);
}

void SwitchNode::ClearTable(){
	m_rtTable.clear();
}

// This function can only be called in switch mode
bool SwitchNode::SwitchReceiveFromDevice(Ptr<NetDevice> device, Ptr<Packet> packet, CustomHeader &ch){
	SendToDev(packet, ch);
	return true;
}

/***********************
 * DNN dcqcn Scheduler
 ***********************/

union sk_key
{
    uint32_t sk_key;
    char cm_sk_key[8];
};

//将uint32_t和uint16_t转变成string

uint64_t SwitchNode::ld_to_uint64(long double value) {
    // 检查是否为负数
    if (value < 0.0L) {
        throw std::invalid_argument("Negative value cannot be converted to uint64_t");
    }

    // 检查是否超出 uint64_t 范围
    if (value > static_cast<long double>(UINT64_MAX)) {
        throw std::overflow_error("Value exceeds uint64_t maximum");
    }

    // 四舍五入到最近的整数（可选）
    // long double rounded = std::round(value);

    // 直接截断小数部分
    return static_cast<uint64_t>(value);
}

std::string SwitchNode::uint64ToString(uint64_t num){
	return std::to_string(num);
}

std::string SwitchNode::uint32ToString(uint32_t num){
	std::stringstream ss;
	ss<<num;
	return ss.str();
}

// uint32_t SwitchNode::stringToUint64(std::string str){
//     try {
//         // 使用std::stoull将字符串转换为uint64_t
//         return std::stoull(str);
//     } catch (const std::invalid_argument& e) {
//         // 处理无效输入的情况，例如字符串不包含数字
//         throw std::invalid_argument("Invalid argument: string does not represent a valid number");
//     } catch (const std::out_of_range& e) {
//         // 处理超出uint64_t范围的情况
//         throw std::out_of_range("Out of range: number is too large for uint64_t");
//     }
// }
/***********************
 * DNN dcqcn Scheduler
 ***********************/

void SwitchNode::SwitchNotifyDequeue(uint32_t ifIndex, uint32_t qIndex, Ptr<Packet> p){
	FlowIdTag t;
	p->PeekPacketTag(t);
	if (qIndex != 0){
		uint32_t inDev = t.GetFlowId();
		m_mmu->RemoveFromIngressAdmission(inDev, qIndex, p->GetSize());
		m_mmu->RemoveFromEgressAdmission(ifIndex, qIndex, p->GetSize());
		m_bytes[inDev][ifIndex][qIndex] -= p->GetSize();
		if (m_ecnEnabled){

			bool egressCongested = m_mmu->ShouldSendCN(ifIndex, qIndex);
			/***********************
 			 * DNN dcqcn Scheduler
			 * The switch collects job information.
             ***********************/
			
			CustomHeader ch(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
			p->PeekHeader(ch);
			//用四元组来标识一条流
			uint8_t ecn = ch.GetIpv4EcnBits();
			std::string flowId = uint32ToString(ch.sip) + uint32ToString(ch.dip);
			bool ECN_flag = false;
			//std::cout<<"[debug] 13Port= "<<ch.l3Prot<<", FLowID="<<inDev<<std::endl;

			if(ch.l3Prot == 0x6 || ch.l3Prot == 0x11){
				//统计过路的TCP和UDP包的信息
				uint32_t pkt_size = p->GetSize();
				uint32_t payload_size = pkt_size - ch.GetSerializedSize();

				sk_key temp_sk_key;
				uint64_t flow_id = atoll(flowId.c_str());// string to unit64_t
				temp_sk_key.sk_key=flow_id;

				uint64_t pkt_size_counter = cms_PktSize.cms_check_mean(&cms_PktSize,temp_sk_key.cm_sk_key);
				uint64_t time_counter = cms_TimeStamp.cms_check_mean(&cms_TimeStamp,temp_sk_key.cm_sk_key);
				
				uint64_t now = Simulator::Now().GetTimeStep();
				//std::cout<<"[debug] ECN?= "<<(ecn == CustomHeader::EcnType::ECN_CE)<<" and Now="<<time_counter<<std::endl;
				if(time_counter != 0 ){//如果开始统计了
					//判断不同的iteration
					uint64_t timeGap_curr; 
					timeGap_curr = now - time_counter;
					//这里因为没有写覆盖的sketch操作，所以先删除再写入
					cms_TimeStamp.cms_remove_inc(&cms_TimeStamp,temp_sk_key.cm_sk_key,time_counter);////清空时间戳
					cms_TimeStamp.cms_add_inc(&cms_TimeStamp,temp_sk_key.cm_sk_key,now);//写入时间戳
					
					//if(timeGap_curr >= 1e5 && pkt_size_counter >=1e5)//在不同的通信阶段且不是ACK
					uint64_t dnn_totalBytes = cms_State.cms_check_mean(&cms_State,temp_sk_key.cm_sk_key);
					//若在不同的iteration，就更新flow的state
					if(timeGap_curr >= 1e5)
					{
						//更新一次迭代的总字节数
						cms_State.cms_remove_inc(&cms_State,temp_sk_key.cm_sk_key,dnn_totalBytes);////清空上次的计数
						cms_State.cms_add_inc(&cms_State,temp_sk_key.cm_sk_key,pkt_size_counter);//写入这次的计数
						
						//开始新一轮的统计
						cms_PktSize.cms_remove_inc(&cms_PktSize,temp_sk_key.cm_sk_key,pkt_size_counter);
						cms_PktSize.cms_add_inc(&cms_PktSize,temp_sk_key.cm_sk_key,payload_size);
					}else{
						cms_PktSize.cms_add_inc(&cms_PktSize,temp_sk_key.cm_sk_key,payload_size);
					}
					
					uint64_t bytes_send = cms_PktSize.cms_check_mean(&cms_PktSize,temp_sk_key.cm_sk_key);
					
					//进行调度和判断
					if(dnn_totalBytes != 0 && !(ecn == Ipv4Header::CE) && ((payload_size * 8e9 / timeGap_curr) > 1e8)){
						long double run_flag = GetPortRunFlag(ifIndex, qIndex);
						uint64_t remainder_bytes;
						if(dnn_totalBytes > bytes_send)
							remainder_bytes = dnn_totalBytes - bytes_send;
						else remainder_bytes=0;
						/************LAS*************/
						////run_flag is running FlowID
						// if(run_flag == 0)
						// 	run_flag=flow_id;
						// else if(run_flag==flow_id){
						// 	if(remainder_bytes<1000){
						// 		run_flag = 0;
						//		SetPortRunFlag(ifIndex, qIndex,run_flag);
						//}
						// }else{
						// 		ECN_flag=true;
						// }
						/************LAS*************/
						// /************SJF*************/
						// // //run_flag is minremainder_bytes
						// if(remainder_bytes < run_flag && remainder_bytes > 1000){
						// 		run_flag = remainder_bytes;
						// 		SetPortRunFlag(ifIndex, qIndex,run_flag);
						// }else{
						// 	if(remainder_bytes <= 1000 ){
						// 		run_flag = 0;
						// 		SetPortRunFlag(ifIndex, qIndex,run_flag);
						// 	}
						// 	else 
						// 		ECN_flag = true;
						// }
						// /************SJF*************/
						/************SJF*************/
						// //run_flag is minremainder_bytes percent
						// double ratio=static_cast<double>(remainder_bytes)/static_cast<double>(dnn_totalBytes);
						// if(run_flag==0||(ratio< run_flag && remainder_bytes > 1000)){
						// 	run_flag = ratio;
						// 	SetPortRunFlag(ifIndex, qIndex,run_flag);
						// }else{
						// 	if(remainder_bytes <= 1500 ){
						// 		run_flag = 0;
						// 		SetPortRunFlag(ifIndex, qIndex,run_flag);
						// 	}
						// 	//else if((payload_size * 1e9 / timeGap_curr) > 6e9)
						// 	else
						// 		ECN_flag = true;
						// }
						/************SJF*************/
						/************CRUX*************/
						//run_flag is running FlowID							
						// if(run_flag==0||run_flag==flow_id){
						// 	if(remainder_bytes>1000){
						// 		run_flag=flow_id;
						// 		SetPortRunFlag(ifIndex, qIndex,run_flag);
						// 	}
						// 	else{
						// 		run_flag=0;
						// 		AddPortTag(ifIndex, qIndex);
						// 		SetPortRunFlag(ifIndex, qIndex,run_flag);
						// 	}
						// }else {
						// 	sk_key temp_last_key;
						// 	temp_last_key.sk_key=ld_to_uint64(run_flag);
						// 	uint64_t last = cms_State.cms_check_mean(&cms_State,temp_last_key.cm_sk_key);
						// 	if(dnn_totalBytes>last){
						// 		if(remainder_bytes>1000){
						// 			run_flag=flow_id;
						// 			SetPortRunFlag(ifIndex, qIndex,run_flag);
						// 		}
						// 	}
						// 	else{
						// 		ECN_flag=true;
						// 	}
						// }
						/************CRUX*************/
						// if(GetPortTag(ifIndex, qIndex)>fair_factor){
						// 	ECN_flag=(!ECN_flag);
						// }else if(GetPortTag(ifIndex, qIndex)>=10)
						// 	RemovePortTag(ifIndex, qIndex);
						
						Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
						//std::cout<<"[debug] qp->sip="<<ch.sip<<", qp->dip="<<ch.dip<<std::endl;
						if(ECN_flag){
							// std::cout<<"[debug1]: flowID = "<<flowId<<",dnn_totalBytes = "<<dnn_totalBytes<<",bytes_send = "<<bytes_send<<std::endl;
							// std::cout<<"[debug2]: remainder_bytes = "<<remainder_bytes<<",run_flag = "<<run_flag<<std::endl;
							// std::cout<<"[debug3]: ECN_flag = "<<ECN_flag<<",now = "<<now<<std::endl;
							// std::cout<<"[debug4]: timeGap_curr="<<timeGap_curr<<std::endl;
							// std::cout<<"[debug5]: time_counter="<<time_counter<<std::endl;
							// std::cout<<"[debug6]: m_node_GetID= "<<device->GetNode()->GetId()<<std::endl;
							// std::cout<<"[debug7]: inDev="<<inDev<<std::endl;
							// std::cout<<"[debug8]: ifIndex, qIndex="<<ifIndex<<", "<<qIndex<<std::endl;
							// std::cout<<"[debug9]: GetPortRunFlag="<<GetPortRunFlag(ifIndex, qIndex)<<std::endl;
							// std::cout<<std::endl;
						}else{
							// std::cout<<"[debug4]: flowID = "<<flowId<<",dnn_totalBytes = "<<dnn_totalBytes<<",bytes_send = "<<bytes_send<<std::endl;
							// std::cout<<"[debug5]: remainder_bytes = "<<remainder_bytes<<",run_flag = "<<run_flag<<std::endl;
							// std::cout<<"[debug6]: ECN_flag = "<<ECN_flag<<",now = "<<now<<std::endl;
							// std::cout<<"[debug]: pkt_ratio = "<<payload_size * 8e9 / timeGap_curr<<std::endl;
							// std::cout<<std::endl;
						}
					}
					std::cout<<"[debug]: pkt_ratio = "<<payload_size * 8e9 / timeGap_curr<<std::endl;
					std::cout<<"[debug]: inDev="<<inDev<<std::endl;
					std::cout<<"[debug]: ifIndex, qIndex="<<ifIndex<<", "<<qIndex<<std::endl<<std::endl;		

				}else{//第一次统计就直接记录
					cms_PktSize.cms_add_inc(&cms_PktSize,temp_sk_key.cm_sk_key,payload_size);
					cms_TimeStamp.cms_remove_inc(&cms_TimeStamp,temp_sk_key.cm_sk_key,time_counter);////清空时间戳
					cms_TimeStamp.cms_add_inc(&cms_TimeStamp,temp_sk_key.cm_sk_key,now);//写入时间戳
					
				}			
			}
			//if (egressCongested || (ECN_flag && m_mmu->ShouldMarkCN(ifIndex, qIndex))){
			if (egressCongested || (ECN_flag)){
			/***********************
			 * DNN dcqcn Scheduler
			 ***********************/
			//if (egressCongested ){
				PppHeader ppp;
				Ipv4Header h;
				p->RemoveHeader(ppp);
				p->RemoveHeader(h);
				h.SetEcn((Ipv4Header::EcnType)0x03);
				p->AddHeader(h);
				p->AddHeader(ppp);
			}
		}
		//CheckAndSendPfc(inDev, qIndex);
		CheckAndSendResume(inDev, qIndex);
	}
	if (1){
		uint8_t* buf = p->GetBuffer();
		if (buf[PppHeader::GetStaticSize() + 9] == 0x11){ // udp packet
			IntHeader *ih = (IntHeader*)&buf[PppHeader::GetStaticSize() + 20 + 8 + 6]; // ppp, ip, udp, SeqTs, INT
			Ptr<QbbNetDevice> dev = DynamicCast<QbbNetDevice>(m_devices[ifIndex]);
			if (m_ccMode == 3){ // HPCC
				ih->PushHop(Simulator::Now().GetTimeStep(), m_txBytes[ifIndex], dev->GetQueue()->GetNBytesTotal(), dev->GetDataRate().GetBitRate());
			}else if (m_ccMode == 10){ // HPCC-PINT
				uint64_t t = Simulator::Now().GetTimeStep();
				uint64_t dt = t - m_lastPktTs[ifIndex];
				if (dt > m_maxRtt)
					dt = m_maxRtt;
				uint64_t B = dev->GetDataRate().GetBitRate() / 8; //Bps
				uint64_t qlen = dev->GetQueue()->GetNBytesTotal();
				double newU;

				/**************************
				 * approximate calc
				 *************************/
				int b = 20, m = 16, l = 20; // see log2apprx's paremeters
				int sft = logres_shift(b,l);
				double fct = 1<<sft; // (multiplication factor corresponding to sft)
				double log_T = log2(m_maxRtt)*fct; // log2(T)*fct
				double log_B = log2(B)*fct; // log2(B)*fct
				double log_1e9 = log2(1e9)*fct; // log2(1e9)*fct
				double qterm = 0;
				double byteTerm = 0;
				double uTerm = 0;
				if ((qlen >> 8) > 0){
					int log_dt = log2apprx(dt, b, m, l); // ~log2(dt)*fct
					int log_qlen = log2apprx(qlen >> 8, b, m, l); // ~log2(qlen / 256)*fct
					qterm = pow(2, (
								log_dt + log_qlen + log_1e9 - log_B - 2*log_T
								)/fct
							) * 256;
					// 2^((log2(dt)*fct+log2(qlen/256)*fct+log2(1e9)*fct-log2(B)*fct-2*log2(T)*fct)/fct)*256 ~= dt*qlen*1e9/(B*T^2)
				}
				if (m_lastPktSize[ifIndex] > 0){
					int byte = m_lastPktSize[ifIndex];
					int log_byte = log2apprx(byte, b, m, l);
					byteTerm = pow(2, (
								log_byte + log_1e9 - log_B - log_T
								)/fct
							);
					// 2^((log2(byte)*fct+log2(1e9)*fct-log2(B)*fct-log2(T)*fct)/fct) ~= byte*1e9 / (B*T)
				}
				if (m_maxRtt > dt && m_u[ifIndex] > 0){
					int log_T_dt = log2apprx(m_maxRtt - dt, b, m, l); // ~log2(T-dt)*fct
					int log_u = log2apprx(int(round(m_u[ifIndex] * 8192)), b, m, l); // ~log2(u*512)*fct
					uTerm = pow(2, (
								log_T_dt + log_u - log_T
								)/fct
							) / 8192;
					// 2^((log2(T-dt)*fct+log2(u*512)*fct-log2(T)*fct)/fct)/512 = (T-dt)*u/T
				}
				newU = qterm+byteTerm+uTerm;

				#if 0
				/**************************
				 * accurate calc
				 *************************/
				double weight_ewma = double(dt) / m_maxRtt;
				double u;
				if (m_lastPktSize[ifIndex] == 0)
					u = 0;
				else{
					double txRate = m_lastPktSize[ifIndex] / double(dt); // B/ns
					u = (qlen / m_maxRtt + txRate) * 1e9 / B;
				}
				newU = m_u[ifIndex] * (1 - weight_ewma) + u * weight_ewma;
				printf(" %lf\n", newU);
				#endif

				/************************
				 * update PINT header
				 ***********************/
				uint16_t power = Pint::encode_u(newU);
				if (power > ih->GetPower())
					ih->SetPower(power);

				m_u[ifIndex] = newU;
			}
		}
	}
	m_txBytes[ifIndex] += p->GetSize();
	m_lastPktSize[ifIndex] = p->GetSize();
	m_lastPktTs[ifIndex] = Simulator::Now().GetTimeStep();
}

int SwitchNode::logres_shift(int b, int l){
	static int data[] = {0,0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};
	return l - data[b];
}

int SwitchNode::log2apprx(int x, int b, int m, int l){
	int x0 = x;
	int msb = int(log2(x)) + 1;
	if (msb > m){
		x = (x >> (msb - m) << (msb - m));
		#if 0
		x += + (1 << (msb - m - 1));
		#else
		int mask = (1 << (msb-m)) - 1;
		if ((x0 & mask) > (rand() & mask))
			x += 1<<(msb-m);
		#endif
	}
	return int(log2(x) * (1<<logres_shift(b, l)));
}

} /* namespace ns3 */
