#include <sstream>
#include <memory>

#include "rss.h"
#include "serializer.h"

namespace rss {

class RSS_Serializer {
public:

    static void ToBytes(const RadixStringSpline& rss, std::string *bytes) {
        std::stringstream buffer;

        uint32_t block_num = rss.accumulate_block_kv_num.size();
        buffer.write(reinterpret_cast<const char*>(&block_num),
                    sizeof(uint32_t));
        for(uint32_t i = 0; i < rss.accumulate_block_kv_num.size(); i++) {
            uint32_t kv_num = rss.accumulate_block_kv_num[i];
            buffer.write(reinterpret_cast<const char*>(&kv_num),
                sizeof(uint32_t));
        }

        ToNodeBytes(rss.root_, buffer);

        bytes->append(buffer.str());
    }

    static RadixStringSpline FromBytes(const std::string& bytes) {
        std::istringstream in(bytes);
        RadixStringSpline rss;

        uint32_t block_num;
        in.read(reinterpret_cast<char*>(&block_num), sizeof(uint32_t));
        rss.accumulate_block_kv_num.resize(block_num);

        for(uint32_t i = 0; i < block_num; i++) {
            uint32_t kv_num;
            in.read(reinterpret_cast<char*>(&kv_num), sizeof(uint32_t));
            rss.accumulate_block_kv_num[i] = kv_num;
        }

        rss.root_.reset(FromNodeBytes(in));

        return std::move(rss);
    }

    static void ToNodeBytes(const std::unique_ptr<RadixStringSplineNode> &rsn, std::stringstream &buffer) {
        const size_t rsn_level = rsn->level_;
        buffer.write(reinterpret_cast<const char*>(&rsn_level),
                    sizeof(size_t));

        const size_t error_bound = rsn->error_bound_;
        buffer.write(reinterpret_cast<const char*>(&error_bound),
                    sizeof(uint16_t));

        buffer.write(reinterpret_cast<const char*>(&rsn->begin_ptr_),
                    sizeof(size_t));
        buffer.write(reinterpret_cast<const char*>(&rsn->end_ptr_),
                    sizeof(size_t));
        buffer.write(reinterpret_cast<const char*>(&rsn->start_range_),
                    sizeof(size_t));
        buffer.write(reinterpret_cast<const char*>(&rsn->end_range_),
                    sizeof(size_t));
        
        size_t size = rsn->block_nums_.size();
        buffer.write(reinterpret_cast<const char*>(&size),
                    sizeof(size_t));
        for(uint32_t i = 0; i < rsn->block_nums_.size(); i++) {
            size = rsn->block_nums_[i].size();
            buffer.write(reinterpret_cast<const char*>(&size),
                sizeof(size_t));
            
            for(uint32_t j = 0; j < rsn->block_nums_[i].size(); j++) {
                buffer.write(reinterpret_cast<const char*>(&rsn->block_nums_[i][j].first),
                sizeof(uint64_t));
                buffer.write(reinterpret_cast<const char*>(&rsn->block_nums_[i][j].second),
                sizeof(size_t));
            }

        }
        
        rs::Serializer<uint64_t>::ToBytes(rsn->rs_model_, buffer);  
        
        const size_t rsn_rmap_size = rsn->redirector_map_.size();
        buffer.write(reinterpret_cast<const char*>(&rsn_rmap_size),
                    sizeof(size_t));
        for (auto& item : rsn->redirector_map_) {
            const size_t encoded_prefix = item.first;
            buffer.write(reinterpret_cast<const char*>(&encoded_prefix),
                    sizeof(size_t));
            
            ToNodeBytes(item.second, buffer);
        }

    }

    static RadixStringSplineNode* FromNodeBytes(std::istringstream &in) {

        auto rsn = new RadixStringSplineNode();
        rsn->prefix_size_ = 8;

        in.read(reinterpret_cast<char*>(&rsn->level_), sizeof(size_t));

        in.read(reinterpret_cast<char*>(&rsn->error_bound_), sizeof(uint16_t));

        in.read(reinterpret_cast<char*>(&rsn->begin_ptr_), sizeof(size_t));
        in.read(reinterpret_cast<char*>(&rsn->end_ptr_), sizeof(size_t));
        in.read(reinterpret_cast<char*>(&rsn->start_range_), sizeof(size_t));
        in.read(reinterpret_cast<char*>(&rsn->end_range_), sizeof(size_t));

        size_t block_num_size;
        in.read(reinterpret_cast<char*>(&block_num_size), sizeof(size_t));
        rsn->block_nums_.resize(block_num_size);

        for(size_t i = 0; i < block_num_size; i++) {
            size_t size;
            in.read(reinterpret_cast<char*>(&size), sizeof(size_t));
            rsn->block_nums_[i].resize(size);
            
            for(size_t j = 0; j < size; j++) {
                uint64_t key_prefix;
                size_t block_num;
                in.read(reinterpret_cast<char*>(&key_prefix), sizeof(uint64_t));
                in.read(reinterpret_cast<char*>(&block_num), sizeof(size_t));
                rsn->block_nums_[i][j] = {key_prefix, block_num};
            }
        }

        rsn->rs_model_ = rs::Serializer<uint64_t>::FromBytes(in);


        size_t rmap_size = 0;
        in.read(reinterpret_cast<char*>(&rmap_size), sizeof(size_t));

        for (size_t i = 0; i < rmap_size; i++) {
            size_t encoded_prefix = 0;
            in.read(reinterpret_cast<char*>(&encoded_prefix), sizeof(size_t));
            rsn->redirector_map_[encoded_prefix].reset(FromNodeBytes(in));
        }

        return rsn;
  }

};

}