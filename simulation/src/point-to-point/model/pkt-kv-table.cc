#include "pkt-kv-table.h"
#include <mutex>

namespace ns3
{
    TIMEKVTABLE::TIMEKVTABLE(){

    }
    void TIMEKVTABLE::insert_value(std::string flow_id,uint64_t time_stamp){
	std::lock_guard<std::mutex> lock(write_mutex);
        tkv_table[flow_id] = time_stamp;
    }

    uint64_t TIMEKVTABLE::get_value(std::string flow_id){
        auto it = tkv_table.find(flow_id);
        if(it != tkv_table.end()){
            return it->second;
        }
        return 0;
    }

    void TIMEKVTABLE::remove_key(std::string flow_id){
        tkv_table.erase(flow_id);
    }
}
