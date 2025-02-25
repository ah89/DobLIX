#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>

#include "builder.h"
#include "string_serializer.h"

#include "rocksdb/slice.h"

namespace rss {

typedef std::vector<std::tuple<std::string, uint64_t, size_t>> key_vector_t;

template<typename T, typename... Args>
std::unique_ptr<T> make_unique_ptr(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

class RadixStringSplineNode {

public:
    RadixStringSplineNode() = default;

    RadixStringSplineNode(const key_vector_t &keys, const size_t begin_ptr, const size_t end_ptr,
    const size_t start_range, const size_t end_range,
    const uint16_t error_bound, const uint32_t max_block_size, const uint8_t prefix_size, const uint16_t level):
        prefix_size_(prefix_size), error_bound_(error_bound), level_(level), max_block_size_(max_block_size),
        begin_ptr_(begin_ptr), end_ptr_(end_ptr),
        start_range_(start_range), end_range_(end_range) {

	    srand((unsigned) time(NULL));

        buildRadixNode(keys);
    }

    rs::SplineModelResult find(const rocksdb::Slice &key) {
        const auto encoded_prefix = StringSerializer::Touint64_t(key.data(), key.size(), begin_ptr_);
        //std::cout <<"find:begin_ptr_:\t" << begin_ptr_ << "\tencoded_prefix:\t" << encoded_prefix << std::endl;
        const auto redirector_map_node = redirector_map_.find(encoded_prefix);
        if (redirector_map_node == redirector_map_.end()) {
            auto rs_output = rs_model_.GetSpline(encoded_prefix);
            //std::cout <<"find:rs_out_spline: " << rs_output.spline << std::endl;
            rs_output.rss_node = this;
            return std::move(rs_output);
        }
        return std::move(redirector_map_node->second->find(key));
    }

    rs::SearchBound GetBlockSearchBound(const rocksdb::Slice &key, const size_t spline_index, const size_t block_start, const size_t block_end) {
        const auto encoded_prefix = StringSerializer::Touint64_t(key.data(), key.size(), begin_ptr_);
        return rs_model_.GetBlockSearchBound(encoded_prefix, spline_index, block_start, block_end);
    }

    void AddBlockNumToSpline(const size_t spline_num, const rocksdb::Slice& key, const size_t block_num) {
        const auto encoded_prefix_key = StringSerializer::Touint64_t(key.data(), key.size(), begin_ptr_);
        block_nums_[spline_num-1].push_back({encoded_prefix_key, block_num});
    }

    size_t GetBlockNum(const rocksdb::Slice& key, const size_t spline_num) {
        const auto encoded_prefix = StringSerializer::Touint64_t(key.data(), key.size(), begin_ptr_);
        size_t iter_index = std::lower_bound(block_nums_[spline_num-1].begin(), block_nums_[spline_num-1].end(), encoded_prefix, comp_pair_uint64()) - block_nums_[spline_num-1].begin();
        iter_index = std::min(iter_index, block_nums_[spline_num-1].size()-1);
        return block_nums_[spline_num-1][iter_index].second;
    }

    size_t GetPrefixLen() const { return begin_ptr_;}

private:
    uint8_t prefix_size_;
    uint16_t error_bound_;
    uint16_t level_;
    uint32_t max_block_size_;

    size_t begin_ptr_;
    size_t end_ptr_;
    
    size_t start_range_;
    size_t end_range_;

    std::map<uint64_t, std::unique_ptr<RadixStringSplineNode>> redirector_map_;
    rs::RadixSpline<uint64_t> rs_model_;

    std::vector< std::vector< std::pair<uint64_t, size_t>>> block_nums_;

    struct comp_pair_uint64 {
        bool operator()(const std::pair<uint64_t, size_t> &a, const uint64_t & b) {
            return (a.first < b);
        }
        bool operator()(const uint64_t & a,const std::pair<uint64_t, size_t> &b) {
            return (a < b.first);
        }
    };


    size_t getNumRadixBits() const {
        switch(level_) {
            case 0:
                return 18;
            case 1:
                return 12;
            default:
                return 6;
        }
    }

    size_t longestPrefixIndex(const std::string& key1, const std::string& key2, size_t begin_index) {
        const auto min_size = std::min(key1.size(), key2.size());

        size_t index = 0;
        for(index = begin_index; index < min_size; index++) {
            if(key1[index] != key2[index]) {
                break;
            }
        }

        return index;

    }

    void buildRadixNode(const key_vector_t &keys) {
        const std::string min_key_str = std::string(std::get<0>(keys[start_range_]).data(), std::get<0>(keys[start_range_]).size());
        const std::string max_key_str = std::string(std::get<0>(keys[end_range_]).data(), std::get<0>(keys[end_range_]).size());
        //std::cout <<"buildRadixNode,min_key_str: " << StringSerializer::Touint64_t(min_key_str.data(), min_key_str.size(), begin_ptr_) << "\tsize:\t" << min_key_str.size() << "\t" << std::get<0>(keys[start_range_]).size() << std::endl;
        // //std::cout <<"buildRadixNode,max_key_str: " << StringSerializer::Touint64_t(max_key_str.data(), max_key_str.size(), begin_ptr_) << std::endl;

        begin_ptr_ = longestPrefixIndex(min_key_str, max_key_str, begin_ptr_);

        const auto min_key = StringSerializer::Touint64_t(min_key_str.data(), min_key_str.size(), begin_ptr_);
        const auto max_key = StringSerializer::Touint64_t(max_key_str.data(), max_key_str.size(), begin_ptr_);

        const auto num_radix_bits = getNumRadixBits();
        //std::cout << "buildRadixNode: begin_ptr:\t" << begin_ptr_ << "\tmin_key:\t" << min_key << "\tmax_key:\t" << max_key << std::endl;
        rs::Builder<uint64_t> rs_builder(min_key, max_key, end_range_, num_radix_bits, error_bound_, max_block_size_);

        auto cluster_encoded_key = min_key;
        auto cluster_start_idx = start_range_;
        for (auto pos = start_range_; pos <= end_range_; pos++) {
            const auto encoded_key = StringSerializer::Touint64_t(std::get<0>(keys[pos]).data(), std::get<0>(keys[pos]).size(), begin_ptr_);
            //std::cout <<"buildRadixNode,encoded_key: " << encoded_key << std::endl;
            if (encoded_key == cluster_encoded_key) {
                continue;
            }

            const size_t cluster_end_offset = std::get<2>(keys[pos-1]);
            const auto median_pos = cluster_start_idx + (pos - cluster_start_idx) / 2;
            //std::cout << "rss:buildRadixNode:cluster_end_offset:\t" << cluster_end_offset << "\tstart_range:\t" << start_range_ << "\tend_range:\t" << end_range_ << "\tpos:\t" << pos << std::endl;
            rs_builder.AddKey(cluster_encoded_key, median_pos, cluster_end_offset);
            cluster_start_idx = pos;
            cluster_encoded_key = encoded_key;
        }
        const size_t cluster_end_offset = std::get<2>(keys[end_range_]);
        //std::cout << "rss:buildRadixNode:cluster_end_offset:\t" << cluster_end_offset << std::endl;
        rs_builder.AddKey(cluster_encoded_key, cluster_start_idx + (end_range_ - cluster_start_idx) / 2, cluster_end_offset);

        rs_model_ = rs_builder.Finalize();
        block_nums_.resize(rs_model_.spline_points_.size());

        create_child_model_if_needed(keys);
    }

    void create_child_model_if_needed(const key_vector_t &keys) {
        auto prev_encoded_prefix = StringSerializer::Touint64_t(std::get<0>(keys[start_range_]).data(), std::get<0>(keys[start_range_]).size(), begin_ptr_);
        size_t cluster_start_idx = start_range_;
        size_t current_idx = start_range_ + 1;

        while (current_idx <= end_range_) {
            const auto current_encoded_prefix = StringSerializer::Touint64_t(std::get<0>(keys[current_idx]).data(), std::get<0>(keys[current_idx]).size(), begin_ptr_);
            
            //skip to reach the end of prefix
            if(prev_encoded_prefix == current_encoded_prefix && current_idx < end_range_) {
                ++current_idx;
                continue;
            }

            const auto cluster_end_idx = (prev_encoded_prefix == current_encoded_prefix && current_idx == end_range_) ? end_range_ : current_idx - 1;

            //unique prefix
            if (cluster_start_idx == current_idx - 1) {
                prev_encoded_prefix = current_encoded_prefix;
                cluster_start_idx = current_idx;
            }

            //prefix with more than 2E dup
            else if (current_idx - cluster_start_idx > 2 * error_bound_) {
                redirector_map_[prev_encoded_prefix] = make_unique_ptr<RadixStringSplineNode>(keys, begin_ptr_ + prefix_size_, end_ptr_,
                                                        cluster_start_idx, cluster_end_idx, error_bound_, max_block_size_, prefix_size_, level_ + 1);

                prev_encoded_prefix = current_encoded_prefix;
                cluster_start_idx = current_idx;
            }

            // verify model output on the first and last key in the prefix
            else {
                bool model_correct = check_model_correctness(prev_encoded_prefix, cluster_start_idx, current_idx - 1);

                if (!model_correct) {
                    check_model_correctness(prev_encoded_prefix, cluster_start_idx, cluster_end_idx);
                    redirector_map_[prev_encoded_prefix] = make_unique_ptr<RadixStringSplineNode>(keys, begin_ptr_ + prefix_size_,
                                                                end_ptr_, cluster_start_idx, cluster_end_idx,
                                                                error_bound_, max_block_size_, prefix_size_, level_ + 1);
                }

                prev_encoded_prefix = current_encoded_prefix;
                cluster_start_idx = current_idx;
            }

            ++current_idx;
            continue;

        }
    }

    inline bool check_model_correctness(const size_t key, const size_t cluster_start, const size_t cluster_end) {
        const auto model_output = rs_model_.GetSearchBound(key);
        if (model_output.begin <= cluster_start && cluster_end <= model_output.end)
            return true;
        return false;
    }

    void showModelTree() {
        //std::cout << "Bound:\t[" << start_range_ << "," << end_range_ << "]" << std::endl;
        //std::cout << "Height:\t" << size_t(level_) << std::endl;
        //std::cout << "Spline Num:\t" << rs_model_.spline_points_.size() << std::endl;

        for (auto& item : redirector_map_) {
            item.second->showModelTree();
        }
    }

    friend class RadixStringSpline;
    friend class RSS_Serializer;

};

class RadixStringSpline {

public:

    RadixStringSpline() = default;

    RadixStringSpline(const key_vector_t &keys, 
    const uint16_t error_bound, const size_t max_key_length, const uint32_t max_block_size) {
        root_ = make_unique_ptr<RadixStringSplineNode>(keys, 0, max_key_length, 0, keys.size() - 1, error_bound, max_block_size, prefix_size_, 0);
        showModelTree();
    }

    inline rs::SplineModelResult find(const rocksdb::Slice& key) {
        return root_->find(key);
    }

    void showModelTree() {
        root_->showModelTree();
    }

    void add_to_accumulate_block_kv_num(uint32_t kv_num) { accumulate_block_kv_num.push_back(kv_num);}

    uint32_t get_accumulate_block_kv_num(uint32_t block_num) { return accumulate_block_kv_num[block_num];}

private:
    uint8_t prefix_size_ = 8;
    std::vector<uint32_t> accumulate_block_kv_num; 
    std::unique_ptr<RadixStringSplineNode> root_;

    friend class RSS_Serializer;
};

}