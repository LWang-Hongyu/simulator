/*
 * @Author: Antgo99 616677561@qq.com
 * @Date: 2024-08-16 10:56:59
 * @LastEditors: Antgo99 616677561@qq.com
 * @LastEditTime: 2024-08-21 15:23:14
 * @FilePath: \CMdcqcn\simulation\src\point-to-point\model\kv-store.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef KV_STORE_H
#define KV_STORE_H

#include <ns3/object.h>
#include <ns3/custom-header.h>
#include <ns3/int-header.h>
#include <ns3/data-rate.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <unordered_map>

namespace ns3{
    class KVSTORE : public Object{
        private:
            std::string flow_id;
            uint32_t pkt_size;
            std::unordered_map<std::string,uint32_t> kv_table;
        public:
            //构造函数
            KVSTORE();

            //插入键值对
            void insert_value(std::string flow_id,uint32_t pkt_size);

            //根据flow id获取pkt_rate
            uint32_t get_value(std::string flow_id);

            //老化机制：按照超过时间戳删除表项
            void remove_key(std::string flow_id);
    };
}

#endif