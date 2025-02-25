#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>

#include "common.h"

namespace rs {

template <class KeyType>
class Serializer;

// Approximates a cumulative distribution function (CDF) using spline
// interpolation.
    template<class KeyType>
    class RadixSpline {
    public:
        RadixSpline() = default;

        RadixSpline(KeyType min_key, KeyType max_key, size_t num_keys,
                    size_t num_radix_bits, size_t num_shift_bits, size_t max_error, size_t max_key_place,
                    std::vector<uint32_t> radix_table,
                    std::vector<rs::Coord<KeyType>> spline_points, std::vector<size_t> spline_max_errors)
                : min_key_(min_key),
                  max_key_(max_key),
                  num_keys_(num_keys),
                  max_error_(max_error),
                  min_key_place_(spline_points[0].y),
                  max_key_place_(max_key_place),
                  num_radix_bits_(num_radix_bits),
                  num_shift_bits_(num_shift_bits),
                  radix_table_(std::move(radix_table)),
                  spline_points_(std::move(spline_points)),
                  spline_max_errors_(std::move(spline_max_errors)) {}

        // Returns the estimated position of `key`.
        double GetEstimatedPosition(const KeyType key) const {
            // Truncate to entity boundaries.
            if (key <= min_key_) return min_key_place_;
            if (key >= max_key_) return max_key_place_ + 1;

            // Find spline segment with `key` ∈ (spline[index - 1], spline[index]].
            const size_t index = GetSplineSegment(key);
            Coord<KeyType> down = spline_points_[index - 1];
            Coord<KeyType> up = spline_points_[index];

            // Compute slope.
            const double x_diff = up.x - down.x;
            const double y_diff = up.y - down.y;
            const double slope = y_diff / x_diff;

            // Interpolate.
            const double key_diff = key - down.x;
            return std::fma(key_diff, slope, down.y);
        }

        // Returns a search bound [begin, end) around the estimated position.
        SearchBound GetSearchBound(const KeyType key) const {
            const size_t estimate = GetEstimatedPosition(key);
            const size_t begin = (estimate < max_error_) ? 0 : (estimate - max_error_);
            // `end` is exclusive.
            const size_t end = (estimate + max_error_ + 2 > max_key_place_)
                               ? max_key_place_
                               : (estimate + max_error_ + 2);
            return SearchBound{begin, end};
        }

        // Returns the size in bytes.
        size_t GetSize() const {
            return sizeof(*this) + radix_table_.size() * sizeof(uint32_t) +
                   spline_points_.size() * sizeof(Coord<KeyType>);
        }

    // Returns the estimated position of `key`.
    SplineModelResult GetSpline(const KeyType key) {
        size_t spline_index;
        // Truncate to entity boundaries.
        if (key <= min_key_) {
            //std::cout  << "GetSpline:min_key_:\t" << min_key_ << "\tkey:\t" << key << std::endl;
            spline_index = 0;
        } else if (key >= max_key_) {
            //std::cout  << "GetSpline:max_key_" << std::endl;
            spline_index = spline_points_.size() - 1;
        } else {
            // Find spline segment with `key` ∈ (spline[index - 1], spline[index]].
            //std::cout  << "GetSpline:GetSplineSegment" << std::endl;
            spline_index = GetSplineSegment(key);
        }

        spline_index = spline_index == 0 ? 1 : spline_index;
        return std::move(SplineModelResult{nullptr, uint32_t(spline_index)});
    }

    rs::SearchBound GetBlockSearchBound(const KeyType key, const size_t spline_index, const size_t block_start, const size_t block_end) {
        double model_output;
        // Truncate to entity boundaries.
        if (key <= min_key_) {
            model_output = block_start;
        } else if (key >= max_key_) {
            model_output = double(block_end);
        } else {
            Coord<KeyType> down = spline_points_[spline_index - 1];
            Coord<KeyType> up = spline_points_[spline_index];

            // Compute slope.
            const double x_diff = up.x - down.x;
            const double y_diff = up.y - down.y;
            const double slope = y_diff / x_diff;

            // Interpolate.
            const double key_diff = key - down.x;
            model_output = std::fma(key_diff, slope, down.y);
        }

        const size_t estimate = model_output > block_start ? model_output - block_start : 0;
        size_t end_of_indexing = block_end - block_start;

        size_t spline_max_error = max_error_; //[spline_index - 1];

        const size_t begin = (estimate < spline_max_error) ? 0 : (estimate - spline_max_error);
        // `end` is exclusive.
        const size_t end = (estimate + spline_max_error + 2 > end_of_indexing)
                               ? end_of_indexing
                               : (estimate + spline_max_error + 2);

        return rs::SearchBound{begin, end};
    }


    KeyType min_key_;
    KeyType max_key_;
    size_t num_keys_;
    size_t max_error_;
    size_t min_key_place_;
    size_t max_key_place_;
    size_t current_accum_block_size = 0;

    size_t GetBlockSegment(const KeyType key) const {
        // Narrow search range using radix table.
        size_t spline_segment = GetSplineSegment(key);
        
        return 0;
        
    }

    private:
        // Returns the index of the spline point that marks the end of the spline
        // segment that contains the `key`: `key` ∈ (spline[index - 1], spline[index]]
        size_t GetSplineSegment(const KeyType key) const {
            // Narrow search range using radix table.
            const KeyType prefix = (key - min_key_) >> num_shift_bits_;
            assert(prefix + 1 < radix_table_.size());
            const uint32_t begin = radix_table_[prefix];
            const uint32_t end = radix_table_[prefix + 1];

            if (end - begin < 32) {
                // Do linear search over narrowed range.
                uint32_t current = begin;
                while (spline_points_[current].x < key) ++current;
                return current;
            }

            // Do binary search over narrowed range.
            const auto lb = std::lower_bound(
                    spline_points_.begin() + begin, spline_points_.begin() + end, key,
                    [](const Coord<KeyType> &coord, const KeyType new_key) {
                        return coord.x < new_key;
                    });
            return std::distance(spline_points_.begin(), lb);
        }

        size_t num_radix_bits_;
        size_t num_shift_bits_;
        std::vector<uint32_t> radix_table_;
        std::vector<rs::Coord<KeyType>> spline_points_;
        std::vector<size_t> spline_max_errors_;

        template<typename>
        friend
        class Serializer;

        friend class RSS_Serializer;
        friend class rss::RadixStringSplineNode;
    };

}  // namespace rs