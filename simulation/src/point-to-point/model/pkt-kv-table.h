/*
 * @Author: Antgo99 616677561@qq.com
 * @Date: 2024-08-19 20:25:21
 * @LastEditors: Antgo99 616677561@qq.com
 * @LastEditTime: 2024-08-21 14:15:30
 * @FilePath: \DNN_Scheduler\simulation\src\point-to-point\model\pkt-kv-table.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef PKT_KV_TABLE_H
#define PKT_KV_TABLE_H

#include <ns3/object.h>
#include <ns3/custom-header.h>
#include <ns3/int-header.h>
#include <ns3/data-rate.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <unordered_map>
#include <utility>
#include <limits>
#include <mutex>

namespace ns3
{
    class TIMEKVTABLE {
        private:
            std::string flow_id;
            uint64_t time_stamp;
            std::unordered_map<std::string,uint64_t> tkv_table;
	    std::mutex write_mutex;
        public:
            //构造函数
            TIMEKVTABLE();

            //插入键值对
            void insert_value(std::string flow_id,uint64_t time_stamp);

            //根据flow id获取pkt_rate
            uint64_t get_value(std::string flow_id);

            //老化机制：按照超过时间戳删除表项
            void remove_key(std::string flow_id);
    };
}

#endif
