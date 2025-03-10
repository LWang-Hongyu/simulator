/*
 * @Author: Antgo99 616677561@qq.com
 * @Date: 2024-08-16 10:57:11
 * @LastEditors: Antgo99 616677561@qq.com
 * @LastEditTime: 2024-08-21 17:43:40
 * @FilePath: \CMdcqcn\simulation\src\point-to-point\model\kv-store.cc
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "kv-store.h"

namespace ns3{
    KVSTORE::KVSTORE(){

    }
    void KVSTORE::insert_value(std::string flow_id,uint32_t pkt_size){
        kv_table[flow_id] = pkt_size;
    }

    uint32_t KVSTORE::get_value(std::string flow_id){
        auto it = kv_table.find(flow_id);
        if(it != kv_table.end()){
            return it->second;
        }
        return 0;
    }

    void KVSTORE::remove_key(std::string flow_id){
        kv_table.erase(flow_id);
    }
}
