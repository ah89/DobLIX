#pragma once

#include <stdint.h>
#include <string>
#include <cstdint>
#include <cstring>

class StringSerializer {
public:
    static uint32_t reversebytes_uint32t(uint32_t value) {
      return (value & 0x000000FFU) << 24 | (value & 0x0000FF00U) << 8 | 
          (value & 0x00FF0000U) >> 8 | (value & 0xFF000000U) >> 24; 
  }

  static uint64_t reversebytes_uint64t(uint64_t value) {
      uint64_t high_uint64 = uint64_t(reversebytes_uint32t(uint32_t(value)));         
      uint64_t low_uint64 = (uint64_t)reversebytes_uint32t(uint32_t(value >> 32));    
      return (high_uint64 << 32) + low_uint64;
  }

  static uint64_t Touint64_t(const char* key_str, const uint32_t key_size, const uint32_t begin_ptr) {
    if (key_size <= begin_ptr) {
      return 0;
    }
    uint64_t lekey = 0;
    size_t lekey_copy_size = sizeof(lekey);
    size_t remaining_size = key_size - begin_ptr;
    if (remaining_size < lekey_copy_size) {
      lekey_copy_size = remaining_size;
    }
    std::memcpy(&lekey, key_str + begin_ptr, lekey_copy_size);
    lekey =  reversebytes_uint64t(lekey);
    return lekey;
  }

};