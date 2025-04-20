#ifndef SWITCH_NODE_H
#define SWITCH_NODE_H

#include <unordered_map>
#include <ns3/node.h>
#include "qbb-net-device.h"
#include "switch-mmu.h"
#include "pint.h"

#include "cmsketch.h"

namespace ns3 {

class Packet;

class SwitchNode : public Node{
	static const uint32_t pCnt = 257;	// Number of ports used
	static const uint32_t qCnt = 8;	// Number of queues/priorities used
	uint32_t m_ecmpSeed;
	std::unordered_map<uint32_t, std::vector<int> > m_rtTable; // map from ip address (u32) to possible ECMP port (index of dev)

	// monitor of PFC
	uint32_t m_bytes[pCnt][pCnt][qCnt]; // m_bytes[inDev][outDev][qidx] is the bytes from inDev enqueued for outDev at qidx
	
	uint64_t m_txBytes[pCnt]; // counter of tx bytes

	uint32_t m_lastPktSize[pCnt];
	uint64_t m_lastPktTs[pCnt]; // ns
	double m_u[pCnt];
	/***********************
	 * DNN dcqcn Scheduler
	 ***********************/
	//设置三个cmsketch 一个是统计数据包的大小 一个是统计时间戳 一个是每次迭代的字节数
	CountMinSketch cms_PktSize;
	CountMinSketch cms_TimeStamp;
	CountMinSketch cms_State;

	//KVSTORE port_flag;
	//uint32_t run_flag;
	//long double run_flag;
	long double run_flag[pCnt][qCnt];
	uint16_t tag[pCnt][qCnt];
	uint16_t fair_factor;

	void SetPortRunFlag(uint32_t port, uint32_t qIndex, long double x){run_flag[port][qIndex] = x;}
	long double GetPortRunFlag(uint32_t port, uint32_t qIndex){return run_flag[port][qIndex];}
	void AddPortTag(uint32_t port, uint32_t qIndex){tag[port][qIndex] += 1;}
	void RemovePortTag(uint32_t port, uint32_t qIndex){tag[port][qIndex] = 0;}
	uint16_t GetPortTag(uint32_t port, uint32_t qIndex){return tag[port][qIndex];}
	/***********************
	 * DNN dcqcn Scheduler
	 ***********************/

protected:
	bool m_ecnEnabled;
	uint32_t m_ccMode;
	uint64_t m_maxRtt;

	uint32_t m_ackHighPrio; // set high priority for ACK/NACK

private:
	int GetOutDev(Ptr<const Packet>, CustomHeader &ch);
	void SendToDev(Ptr<Packet>p, CustomHeader &ch);
	static uint32_t EcmpHash(const uint8_t* key, size_t len, uint32_t seed);
	void CheckAndSendPfc(uint32_t inDev, uint32_t qIndex);
	void CheckAndSendResume(uint32_t inDev, uint32_t qIndex);
public:
	Ptr<SwitchMmu> m_mmu;

	static TypeId GetTypeId (void);
	SwitchNode();
	virtual ~SwitchNode() override;
	void SetEcmpSeed(uint32_t seed);
	void AddTableEntry(Ipv4Address &dstAddr, uint32_t intf_idx);
	void ClearTable();
	bool SwitchReceiveFromDevice(Ptr<NetDevice> device, Ptr<Packet> packet, CustomHeader &ch);
	void SwitchNotifyDequeue(uint32_t ifIndex, uint32_t qIndex, Ptr<Packet> p);

	/***********************
	 * DNN dcqcn Scheduler
	 ***********************/
	std::string uint64ToString(uint64_t num);
	std::string uint32ToString(uint32_t num);
	uint64_t ld_to_uint64(long double value);
	// uint32_t stringToUint64(std::string str);
	/***********************
	 * DNN dcqcn Scheduler
	 ***********************/

	// for approximate calc in PINT
	int logres_shift(int b, int l);
	int log2apprx(int x, int b, int m, int l); // given x of at most b bits, use most significant m bits of x, calc the result in l bits
};

} /* namespace ns3 */

#endif /* SWITCH_NODE_H */
