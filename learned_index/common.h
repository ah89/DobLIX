#pragma once

#include <cstddef>
#include <cstdint>

//Forward decleration
namespace rss {
  class RadixStringSplineNode;
}

namespace rs {

// A CDF coordinate.
template <class KeyType>
struct Coord {
  KeyType x;
  double y;
  size_t offset=0;
};

struct SearchBound {
  size_t begin;
  size_t end;  // Exclusive.
};

struct SplineModelResult {
    rss::RadixStringSplineNode* rss_node;
    uint32_t spline;
};


}  // namespace rs
