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
#include "kv-store.h"
#include "pkt-kv-table.h"
#include <typeinfo>
#include <sstream>
#include <stdexcept> // 包含异常处理

#define DEPTH 7
#define WIDTH 100000

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
	normal_ECN=0;
	my_ECN=0;
	cms_PktSize.cms_init(&cms_PktSize, WIDTH, DEPTH);
	cms_TimeStamp.cms_init(&cms_TimeStamp,WIDTH, DEPTH);
	fflush(stdout);
	run_flag = INT32_MAX;
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
}

SwitchNode::~SwitchNode(){
	std::cout<<"[debug:] normal ECN="<<normal_ECN<<" my_ECN="<<my_ECN<<std::endl;
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

union pktSize_key
{
    uint32_t pktsizekey;
    char cm_PktSize_key[8];
};

//将uint32_t和uint16_t转变成string
std::string SwitchNode::uint32ToString(uint32_t num){
	std::stringstream ss;
	ss<<num;
	return ss.str();
}

std::string SwitchNode::uint16ToString(uint16_t num){
	std::stringstream ss;
	ss<<num;
	return ss.str();
}

uint32_t SwitchNode::stringToUint32(std::string str){
	try{
		return static_cast<uint32_t>(std::stoul(str));
	} catch (const std::out_of_range& e){
		// 如果数值超出了unsigned long的范围，将抛出out_of_range异常
        std::cerr << "Out of range: " << e.what() << '\n';
	}
	return 0; // 如果发生异常，返回0或选择其他错误处理方式
}
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

			std::string flowId = uint32ToString(ch.sip) + uint32ToString(ch.dip);
			bool ECN_flag = false;
			if(ch.sip < ch.dip){
				pktSize_key temp_pktSize_key;
                             
				uint32_t pkt_size = p->GetSize();
				uint32_t payload_size = pkt_size - ch.GetSerializedSize();

				uint32_t flow_id = stringToUint32(flowId);
				temp_pktSize_key.pktsizekey = flow_id;//存储数据包大小的sketch

				uint64_t* pkt_hashes = cms_PktSize.cms_get_hashes_alt(&cms_PktSize,cms_PktSize.depth,temp_pktSize_key.cm_PktSize_key);
				uint32_t pkt_size_counter = cms_PktSize.cms_check_mean(&cms_PktSize,temp_pktSize_key.cm_PktSize_key);
				uint64_t time_size_counter = kvt.get_value(flowId);

				Time now_time = Simulator::Now();
				if(time_size_counter != 0 ){
					uint64_t timeGap_curr;
					timeGap_curr = now_time.GetTimeStep() - time_size_counter;

					kvt.insert_value(flowId,now_time.GetTimeStep());
					if(timeGap_curr >= 1e6)//在不同的通信阶段
					{
						kvs.insert_value(flowId,pkt_size_counter);//更新迭代总字节数
						cms_PktSize.cms_remove_inc_alt(&cms_PktSize,pkt_hashes,cms_PktSize.depth,pkt_size_counter);//清空存储数据包大小的sketch
						cms_PktSize.cms_add_inc_alt(&cms_PktSize,pkt_hashes,cms_PktSize.depth,payload_size);
					}else{
						cms_PktSize.cms_add_inc_alt(&cms_PktSize,pkt_hashes,cms_PktSize.depth,payload_size);
					}

					uint32_t bytes_send = cms_PktSize.cms_check_mean(&cms_PktSize,temp_pktSize_key.cm_PktSize_key);

					uint32_t dnn_totalBytes = kvs.get_value(flowId);
					if(dnn_totalBytes != 0 ){
						uint32_t remainder_bytes;
						if(dnn_totalBytes > bytes_send)
							remainder_bytes = dnn_totalBytes - bytes_send;
						else remainder_bytes=0;
						/************LAS*************/
						////run_flag is running FlowID
						// if(run_flag == INT32_MAX)
						// 	run_flag=flow_id;
						// else if(run_flag==flow_id){
						// 	if(remainder_bytes<1000)
						// 		run_flag = INT32_MAX;
						// }else{
						// 	if((payload_size * 1e9 / timeGap_curr) > 6e9)
						// 		ECN_flag=true;
						// }
						/************LAS*************/
						// /************SJF*************/
						// // //run_flag is minremainder_bytes
						// if(remainder_bytes < run_flag && remainder_bytes > 1000){
						// 	run_flag = remainder_bytes;
						// }else{
						// 	if(remainder_bytes <= 1000 ){
						// 		run_flag = INT32_MAX;
						// 	}
						// 	else if((payload_size * 1e9 / timeGap_curr) > 6e9)
						// 		ECN_flag = true;
						// }
						// /************SJF*************/
						/************SJF*************/
						// //run_flag is minremainder_bytes percent
						double ratio=static_cast<double>(remainder_bytes)/static_cast<double>(dnn_totalBytes);
						if(ratio< run_flag && remainder_bytes > 1000){
							run_flag = ratio;
						}else{
							if(remainder_bytes <= 1000 ){
								run_flag = INT32_MAX;
							}
							//else if((payload_size * 1e9 / timeGap_curr) > 6e9)
							else
								ECN_flag = true;
						}
						/************SJF*************/
						if(ECN_flag){
							// std::cout<<"[debug1]: flowID = "<<flowId<<",dnn_totalBytes = "<<dnn_totalBytes<<",bytes_send = "<<bytes_send<<std::endl;
							// std::cout<<"[debug2]: remainder_bytes = "<<remainder_bytes<<",run_flag = "<<run_flag<<std::endl;
							// std::cout<<"[debug3]: ECN_flag = "<<ECN_flag<<",now_time = "<<now_time.GetTimeStep()<<std::endl;
							// std::cout<<"[debug4]: timeGap_curr="<<timeGap_curr<<std::endl;
							// std::cout<<"[debug5]: time_size_counter="<<time_size_counter<<std::endl;
							// std::cout<<std::endl;
						}else{
							// std::cout<<"[debug4]: flowID = "<<flowId<<",dnn_totalBytes = "<<dnn_totalBytes<<",bytes_send = "<<bytes_send<<std::endl;
							// std::cout<<"[debug5]: remainder_bytes = "<<remainder_bytes<<",run_flag = "<<run_flag<<std::endl;
							// std::cout<<"[debug6]: ECN_flag = "<<ECN_flag<<",now_time = "<<now_time.GetTimeStep()<<std::endl;
							// std::cout<<"[debug]: pkt_ratio = "<<payload_size * 1e9 / timeGap_curr<<std::endl;
							// std::cout<<std::endl;
						}
					}
				}else{
					cms_PktSize.cms_add_inc_alt(&cms_PktSize,pkt_hashes,cms_PktSize.depth,payload_size);
					kvt.insert_value(flowId,now_time.GetTimeStep());
				}
			}
			//if (egressCongested || (ECN_flag && m_mmu->ShouldMarkCN(ifIndex, qIndex))){
			if (egressCongested || (ECN_flag)){
				if(egressCongested) normal_ECN++;
					else my_ECN++;
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
