//  Copyright (c) 2011-present, Facebook, Inc.  All rights reserved.
//  This source code is licensed under the BSD-style license found in the
//  LICENSE file in the root directory of this source tree. An additional grant
//  of patent rights can be found in the PATENTS file in the same directory.
//
// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "table/block_based_table_builder.h"

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "db/dbformat.h"

#include "rocksdb/cache.h"
#include "rocksdb/comparator.h"
#include "rocksdb/env.h"
#include "rocksdb/filter_policy.h"
#include "rocksdb/flush_block_policy.h"
#include "rocksdb/merge_operator.h"
#include "rocksdb/table.h"

#include "table/block.h"
#include "table/block_based_table_reader.h"
#include "table/block_builder.h"
#include "table/filter_block.h"
#include "table/block_based_filter_block.h"
#include "table/block_based_table_factory.h"
#include "table/full_filter_block.h"
#include "table/format.h"
#include "table/meta_blocks.h"
#include "table/table_builder.h"

#include "learned_index/rss.h"
#include "learned_index/rss_serializer.h"
#include "learned_index/string_serializer.h"

#include "util/string_util.h"
#include "util/coding.h"
#include "util/compression.h"
#include "util/crc32c.h"
#include "util/stop_watch.h"
#include "util/xxhash.h"

namespace rocksdb {

extern const std::string kHashIndexPrefixesBlock;
extern const std::string kHashIndexPrefixesMetadataBlock;

typedef BlockBasedTableOptions::IndexType IndexType;

// The interface for building index.
// Instruction for adding a new concrete IndexBuilder:
//  1. Create a subclass instantiated from IndexBuilder.
//  2. Add a new entry associated with that subclass in TableOptions::IndexType.
//  3. Add a create function for the new subclass in CreateIndexBuilder.
// Note: we can devise more advanced design to simplify the process for adding
// new subclass, which will, on the other hand, increase the code complexity and
// catch unwanted attention from readers. Given that we won't add/change
// indexes frequently, it makes sense to just embrace a more straightforward
// design that just works.
class IndexBuilder {
 public:
  // Index builder will construct a set of blocks which contain:
  //  1. One primary index block.
  //  2. (Optional) a set of metablocks that contains the metadata of the
  //     primary index.
  struct IndexBlocks {
    Slice index_block_contents;

    std::unordered_map<std::string, Slice> meta_blocks;
  };
  explicit IndexBuilder(const InternalKeyComparator* comparator)
      : comparator_(comparator) {}

  virtual ~IndexBuilder() {}

  // Add a new index entry to index block.
  // To allow further optimization, we provide `last_key_in_current_block` and
  // `first_key_in_next_block`, based on which the specific implementation can
  // determine the best index key to be used for the index block.
  // @last_key_in_current_block: this parameter maybe overridden with the value
  //                             "substitute key".
  // @first_key_in_next_block: it will be nullptr if the entry being added is
  //                           the last one in the table
  //
  // REQUIRES: Finish() has not yet been called.
  virtual void AddIndexEntry(std::string* last_key_in_current_block,
                             const Slice* first_key_in_next_block,
                             const BlockHandle& block_handle) = 0;

  // This method will be called whenever a key is added. The subclasses may
  // override OnKeyAdded() if they need to collect additional information.
  virtual void OnKeyAdded(const Slice& key) {}

  virtual void OnKeyAdded(const Slice& key, uint64_t offset) {}


  // Inform the index builder that all entries has been written. Block builder
  // may therefore perform any operation required for block finalization.
  //
  // REQUIRES: Finish() has not yet been called.
  virtual Status Finish(IndexBlocks* index_blocks) = 0;

  // Get the estimated size for index block.
  virtual size_t EstimatedSize() const = 0;

 protected:
  const InternalKeyComparator* comparator_;
};

// This index builder builds space-efficient index block.
//
// Optimizations:
//  1. Made block's `block_restart_interval` to be 1, which will avoid linear
//     search when doing index lookup (can be disabled by setting
//     index_block_restart_interval).
//  2. Shorten the key length for index block. Other than honestly using the
//     last key in the data block as the index key, we instead find a shortest
//     substitute key that serves the same function.
class ShortenedIndexBuilder : public IndexBuilder {
 public:
  explicit ShortenedIndexBuilder(const InternalKeyComparator* comparator,
                                 int index_block_restart_interval)
      : IndexBuilder(comparator),
        index_block_builder_(index_block_restart_interval) {}

  virtual void AddIndexEntry(std::string* last_key_in_current_block,
                             const Slice* first_key_in_next_block,
                             const BlockHandle& block_handle) override {
    if (first_key_in_next_block != nullptr) {
      comparator_->FindShortestSeparator(last_key_in_current_block,
                                         *first_key_in_next_block);
    } else {
      comparator_->FindShortSuccessor(last_key_in_current_block);
    }

    std::string handle_encoding;
    block_handle.EncodeTo(&handle_encoding);
    

    auto s = new Slice(handle_encoding);
    BlockHandle handle;
    auto status = handle.DecodeFrom(s);

    index_block_builder_.Add(*last_key_in_current_block, handle_encoding);
  }

  virtual Status Finish(IndexBlocks* index_blocks) override {
    index_blocks->index_block_contents = index_block_builder_.Finish();
    return Status::OK();
  }

  virtual size_t EstimatedSize() const override {
    return index_block_builder_.CurrentSizeEstimate();
  }

 private:
  BlockBuilder index_block_builder_;
};

// HashIndexBuilder contains a binary-searchable primary index and the
// metadata for secondary hash index construction.
// The metadata for hash index consists two parts:
//  - a metablock that compactly contains a sequence of prefixes. All prefixes
//    are stored consectively without any metadata (like, prefix sizes) being
//    stored, which is kept in the other metablock.
//  - a metablock contains the metadata of the prefixes, including prefix size,
//    restart index and number of block it spans. The format looks like:
//
// +-----------------+---------------------------+---------------------+ <=prefix 1
// | length: 4 bytes | restart interval: 4 bytes | num-blocks: 4 bytes |
// +-----------------+---------------------------+---------------------+ <=prefix 2
// | length: 4 bytes | restart interval: 4 bytes | num-blocks: 4 bytes |
// +-----------------+---------------------------+---------------------+
// |                                                                   |
// | ....                                                              |
// |                                                                   |
// +-----------------+---------------------------+---------------------+ <=prefix n
// | length: 4 bytes | restart interval: 4 bytes | num-blocks: 4 bytes |
// +-----------------+---------------------------+---------------------+
//
// The reason of separating these two metablocks is to enable the efficiently
// reuse the first metablock during hash index construction without unnecessary
// data copy or small heap allocations for prefixes.
class HashIndexBuilder : public IndexBuilder {
 public:
  explicit HashIndexBuilder(const InternalKeyComparator* comparator,
                            const SliceTransform* hash_key_extractor,
                            int index_block_restart_interval)
      : IndexBuilder(comparator),
        primary_index_builder_(comparator, index_block_restart_interval),
        hash_key_extractor_(hash_key_extractor) {}

  virtual void AddIndexEntry(std::string* last_key_in_current_block,
                             const Slice* first_key_in_next_block,
                             const BlockHandle& block_handle) override {
    ++current_restart_index_;
    primary_index_builder_.AddIndexEntry(last_key_in_current_block,
                                        first_key_in_next_block, block_handle);
  }

  virtual void OnKeyAdded(const Slice& key, uint64_t offset) override {
      OnKeyAdded(key);
  }

  virtual void OnKeyAdded(const Slice& key) override {
    auto key_prefix = hash_key_extractor_->Transform(key);
    bool is_first_entry = pending_block_num_ == 0;

    // Keys may share the prefix
    if (is_first_entry || pending_entry_prefix_ != key_prefix) {
      if (!is_first_entry) {
        FlushPendingPrefix();
      }

      // need a hard copy otherwise the underlying data changes all the time.
      // TODO(kailiu) ToString() is expensive. We may speed up can avoid data
      // copy.
      pending_entry_prefix_ = key_prefix.ToString();
      pending_block_num_ = 1;
      pending_entry_index_ = static_cast<uint32_t>(current_restart_index_);
    } else {
      // entry number increments when keys share the prefix reside in
      // different data blocks.
      auto last_restart_index = pending_entry_index_ + pending_block_num_ - 1;
      assert(last_restart_index <= current_restart_index_);
      if (last_restart_index != current_restart_index_) {
        ++pending_block_num_;
      }
    }
  }

  virtual Status Finish(IndexBlocks* index_blocks) override {
    FlushPendingPrefix();
    primary_index_builder_.Finish(index_blocks);
    index_blocks->meta_blocks.insert(
        {kHashIndexPrefixesBlock.c_str(), prefix_block_});
    index_blocks->meta_blocks.insert(
        {kHashIndexPrefixesMetadataBlock.c_str(), prefix_meta_block_});
    return Status::OK();
  }

  virtual size_t EstimatedSize() const override {
    return primary_index_builder_.EstimatedSize() + prefix_block_.size() +
           prefix_meta_block_.size();
  }

 private:
  void FlushPendingPrefix() {
    prefix_block_.append(pending_entry_prefix_.data(),
                         pending_entry_prefix_.size());
    PutVarint32Varint32Varint32(
        &prefix_meta_block_,
        static_cast<uint32_t>(pending_entry_prefix_.size()),
        pending_entry_index_, pending_block_num_);
  }

  ShortenedIndexBuilder primary_index_builder_;
  const SliceTransform* hash_key_extractor_;

  // stores a sequence of prefixes
  std::string prefix_block_;
  // stores the metadata of prefixes
  std::string prefix_meta_block_;

  // The following 3 variables keeps unflushed prefix and its metadata.
  // The details of block_num and entry_index can be found in
  // "block_hash_index.{h,cc}"
  uint32_t pending_block_num_ = 0;
  uint32_t pending_entry_index_ = 0;
  std::string pending_entry_prefix_;

  uint64_t current_restart_index_ = 0;
};

class LearnedIndexBuilder : public IndexBuilder {
   public:
  explicit LearnedIndexBuilder(const InternalKeyComparator* comparator = nullptr) : 
  IndexBuilder(comparator), max_key_place(0) {}
 
  virtual void AddIndexEntry(std::string* last_key_in_current_block,
                             const Slice* first_key_in_next_block,
                             const BlockHandle& block_handle) override {
  }
 
  virtual void OnKeyAdded(const Slice& key) override {
    
  }
 
  virtual void OnKeyAdded(const Slice& key, size_t offset) override {
    const auto& user_key = ExtractUserKey(key);
    // std::cout <<"rss add: " << StringSerializer::Touint64_t(user_key.data(), user_key.size(), 0) << "," << StringSerializer::Touint64_t(user_key.data(), user_key.size(), 8) << std::endl;
    // std::cout << "OnKeyAdded:kv_num:\t" << kv_num << "\toffset:\t" << offset << std::endl;
    spline_data.push_back(std::make_tuple(std::move(user_key.ToString()), kv_num, offset));
    max_key_len = std::max(user_key.size(), max_key_len);
    max_key_place = kv_num;
    ++kv_num;
  }
 
 
 
  virtual Status Finish(IndexBlocks* index_blocks) override {
    rss::RSS_Serializer::ToBytes(rss_learned_model, &serialized_model);
 
    index_blocks->index_block_contents = std::move(serialized_model);
    return Status::OK();
  }
 
  virtual size_t EstimatedSize() const override {
    return 0;
  }
 
  rss::RadixStringSpline* FinalizeRSSModel(const int max_error, const int max_block_size) {
    rss_learned_model = rss::RadixStringSpline(spline_data, max_error, max_key_len, max_block_size);

    return &rss_learned_model;
  }
 
private:
  std::vector<std::tuple<std::string, uint64_t, size_t>> spline_data;
  size_t max_key_place;
  size_t max_key_len = 0;
  size_t kv_num = 0;
  std::string serialized_model;
  rss::RadixStringSpline rss_learned_model;
};

// Without anonymous namespace here, we fail the warning -Wmissing-prototypes
namespace {

// Create a index builder based on its type.
IndexBuilder* CreateIndexBuilder(IndexType type,
                                 const InternalKeyComparator* comparator,
                                 const SliceTransform* prefix_extractor,
                                 int index_block_restart_interval) {
  switch (type) {
    case BlockBasedTableOptions::kBinarySearch: {
      return new ShortenedIndexBuilder(comparator,
                                       index_block_restart_interval);
    }
    case BlockBasedTableOptions::kHashSearch: {
      return new HashIndexBuilder(comparator, prefix_extractor,
                                  index_block_restart_interval);
    }
    case BlockBasedTableOptions::kLearnedIndexModel: {
      return new LearnedIndexBuilder(comparator);
    }

    default: {
      assert(!"Do not recognize the index type ");
      return nullptr;
    }
  }
  // impossible.
  assert(false);
  return nullptr;
}

// Create a index builder based on its type.
FilterBlockBuilder* CreateFilterBlockBuilder(const ImmutableCFOptions& opt,
    const BlockBasedTableOptions& table_opt) {
  if (table_opt.filter_policy == nullptr) return nullptr;

  FilterBitsBuilder* filter_bits_builder =
      table_opt.filter_policy->GetFilterBitsBuilder();
  if (filter_bits_builder == nullptr) {
    return new BlockBasedFilterBlockBuilder(opt.prefix_extractor, table_opt);
  } else {
    return new FullFilterBlockBuilder(opt.prefix_extractor,
                                      table_opt.whole_key_filtering,
                                      filter_bits_builder);
  }
}

bool GoodCompressionRatio(size_t compressed_size, size_t raw_size) {
  // Check to see if compressed less than 12.5%
  return compressed_size < raw_size - (raw_size / 8u);
}

}  // namespace

// format_version is the block format as defined in include/rocksdb/table.h
Slice CompressBlock(const Slice& raw,
                    const CompressionOptions& compression_options,
                    CompressionType* type, uint32_t format_version,
                    const Slice& compression_dict,
                    std::string* compressed_output) {
  if (*type == kNoCompression) {
    return raw;
  }

  // Will return compressed block contents if (1) the compression method is
  // supported in this platform and (2) the compression rate is "good enough".
  switch (*type) {
    case kSnappyCompression:
      if (Snappy_Compress(compression_options, raw.data(), raw.size(),
                          compressed_output) &&
          GoodCompressionRatio(compressed_output->size(), raw.size())) {
        return *compressed_output;
      }
      break;  // fall back to no compression.
    case kZlibCompression:
      if (Zlib_Compress(
              compression_options,
              GetCompressFormatForVersion(kZlibCompression, format_version),
              raw.data(), raw.size(), compressed_output, compression_dict) &&
          GoodCompressionRatio(compressed_output->size(), raw.size())) {
        return *compressed_output;
      }
      break;  // fall back to no compression.
    case kBZip2Compression:
      if (BZip2_Compress(
              compression_options,
              GetCompressFormatForVersion(kBZip2Compression, format_version),
              raw.data(), raw.size(), compressed_output) &&
          GoodCompressionRatio(compressed_output->size(), raw.size())) {
        return *compressed_output;
      }
      break;  // fall back to no compression.
    case kLZ4Compression:
      if (LZ4_Compress(
              compression_options,
              GetCompressFormatForVersion(kLZ4Compression, format_version),
              raw.data(), raw.size(), compressed_output, compression_dict) &&
          GoodCompressionRatio(compressed_output->size(), raw.size())) {
        return *compressed_output;
      }
      break;  // fall back to no compression.
    case kLZ4HCCompression:
      if (LZ4HC_Compress(
              compression_options,
              GetCompressFormatForVersion(kLZ4HCCompression, format_version),
              raw.data(), raw.size(), compressed_output, compression_dict) &&
          GoodCompressionRatio(compressed_output->size(), raw.size())) {
        return *compressed_output;
      }
      break;     // fall back to no compression.
    case kXpressCompression:
      if (XPRESS_Compress(raw.data(), raw.size(),
          compressed_output) &&
          GoodCompressionRatio(compressed_output->size(), raw.size())) {
        return *compressed_output;
      }
      break;
    case kZSTD:
    case kZSTDNotFinalCompression:
      if (ZSTD_Compress(compression_options, raw.data(), raw.size(),
                        compressed_output, compression_dict) &&
          GoodCompressionRatio(compressed_output->size(), raw.size())) {
        return *compressed_output;
      }
      break;     // fall back to no compression.
    default: {}  // Do not recognize this compression type
  }

  // Compression method is not supported, or not good compression ratio, so just
  // fall back to uncompressed form.
  *type = kNoCompression;
  return raw;
}

// kBlockBasedTableMagicNumber was picked by running
//    echo rocksdb.table.block_based | sha1sum
// and taking the leading 64 bits.
// Please note that kBlockBasedTableMagicNumber may also be accessed by other
// .cc files
// for that reason we declare it extern in the header but to get the space
// allocated
// it must be not extern in one place.
const uint64_t kBlockBasedTableMagicNumber = 0x88e241b785f4cff7ull;
// We also support reading and writing legacy block based table format (for
// backwards compatibility)
const uint64_t kLegacyBlockBasedTableMagicNumber = 0xdb4775248b80fb57ull;

// A collector that collects properties of interest to block-based table.
// For now this class looks heavy-weight since we only write one additional
// property.
// But in the foreseeable future, we will add more and more properties that are
// specific to block-based table.
class BlockBasedTableBuilder::BlockBasedTablePropertiesCollector
    : public IntTblPropCollector {
 public:
  explicit BlockBasedTablePropertiesCollector(
      BlockBasedTableOptions::IndexType index_type, bool whole_key_filtering,
      bool prefix_filtering)
      : index_type_(index_type),
        whole_key_filtering_(whole_key_filtering),
        prefix_filtering_(prefix_filtering) {}

  virtual Status InternalAdd(const Slice& key, const Slice& value,
                             uint64_t file_size) override {
    // Intentionally left blank. Have no interest in collecting stats for
    // individual key/value pairs.
    return Status::OK();
  }

  virtual Status Finish(UserCollectedProperties* properties) override {
    std::string val;
    PutFixed32(&val, static_cast<uint32_t>(index_type_));
    properties->insert({BlockBasedTablePropertyNames::kIndexType, val});
    properties->insert({BlockBasedTablePropertyNames::kWholeKeyFiltering,
                        whole_key_filtering_ ? kPropTrue : kPropFalse});
    properties->insert({BlockBasedTablePropertyNames::kPrefixFiltering,
                        prefix_filtering_ ? kPropTrue : kPropFalse});
    return Status::OK();
  }

  // The name of the properties collector can be used for debugging purpose.
  virtual const char* Name() const override {
    return "BlockBasedTablePropertiesCollector";
  }

  virtual UserCollectedProperties GetReadableProperties() const override {
    // Intentionally left blank.
    return UserCollectedProperties();
  }

 private:
  BlockBasedTableOptions::IndexType index_type_;
  bool whole_key_filtering_;
  bool prefix_filtering_;
};

struct BlockBasedTableBuilder::Rep {
  const ImmutableCFOptions ioptions;
  const BlockBasedTableOptions table_options;
  const InternalKeyComparator& internal_comparator;
  WritableFileWriter* file;
  uint64_t offset = 0;
  Status status;
  BlockBuilder data_block;
  BlockBuilder range_del_block;

  std::vector<std::pair<std::string, std::string>> kv_data;
  std::unique_ptr<IndexBuilder> learned_index_builder;


  InternalKeySliceTransform internal_prefix_transform;
  std::unique_ptr<IndexBuilder> index_builder;

  std::string last_key;
  const CompressionType compression_type;
  const CompressionOptions compression_opts;
  // Data for presetting the compression library's dictionary, or nullptr.
  const std::string* compression_dict;
  TableProperties props;

  bool closed = false;  // Either Finish() or Abandon() has been called.
  std::unique_ptr<FilterBlockBuilder> filter_block;
  char compressed_cache_key_prefix[BlockBasedTable::kMaxCacheKeyPrefixSize];
  size_t compressed_cache_key_prefix_size;

  BlockHandle pending_handle;  // Handle to add to index block

  std::string compressed_output;
  std::unique_ptr<FlushBlockPolicy> flush_block_policy;
  uint32_t column_family_id;
  const std::string& column_family_name;

  std::vector<std::unique_ptr<IntTblPropCollector>> table_properties_collectors;

  Rep(const ImmutableCFOptions& _ioptions,
      const BlockBasedTableOptions& table_opt,
      const InternalKeyComparator& icomparator,
      const std::vector<std::unique_ptr<IntTblPropCollectorFactory>>*
          int_tbl_prop_collector_factories,
      uint32_t _column_family_id, WritableFileWriter* f,
      const CompressionType _compression_type,
      const CompressionOptions& _compression_opts,
      const std::string* _compression_dict, const bool skip_filters,
      const std::string& _column_family_name)
      : ioptions(_ioptions),
        table_options(table_opt),
        internal_comparator(icomparator),
        file(f),
        data_block(table_options.block_restart_interval,
                   table_options.use_delta_encoding),
        range_del_block(1),  // TODO(andrewkr): restart_interval unnecessary
        internal_prefix_transform(_ioptions.prefix_extractor),
        index_builder(
            CreateIndexBuilder(table_options.index_type, &internal_comparator,
                               &this->internal_prefix_transform,
                               table_options.index_block_restart_interval)),
        learned_index_builder(
            CreateIndexBuilder(BlockBasedTableOptions::kLearnedIndexModel, &internal_comparator,
                      &this->internal_prefix_transform,
                      table_options.index_block_restart_interval)),
        compression_type(_compression_type),
        compression_opts(_compression_opts),
        compression_dict(_compression_dict),
        filter_block(skip_filters ? nullptr : CreateFilterBlockBuilder(
                                                  _ioptions, table_options)),
        flush_block_policy(
            table_options.flush_block_policy_factory->NewFlushBlockPolicy(
                table_options, data_block)),
        column_family_id(_column_family_id),
        column_family_name(_column_family_name) {
    for (auto& collector_factories : *int_tbl_prop_collector_factories) {
      table_properties_collectors.emplace_back(
          collector_factories->CreateIntTblPropCollector(column_family_id));
    }
    table_properties_collectors.emplace_back(
        new BlockBasedTablePropertiesCollector(
            table_options.index_type, table_options.whole_key_filtering,
            _ioptions.prefix_extractor != nullptr));
  }
};

BlockBasedTableBuilder::BlockBasedTableBuilder(
    const ImmutableCFOptions& ioptions,
    const BlockBasedTableOptions& table_options,
    const InternalKeyComparator& internal_comparator,
    const std::vector<std::unique_ptr<IntTblPropCollectorFactory>>*
        int_tbl_prop_collector_factories,
    uint32_t column_family_id, WritableFileWriter* file,
    const CompressionType compression_type,
    const CompressionOptions& compression_opts,
    const std::string* compression_dict, const bool skip_filters,
    const std::string& column_family_name) {
  BlockBasedTableOptions sanitized_table_options(table_options);
  if (sanitized_table_options.format_version == 0 &&
      sanitized_table_options.checksum != kCRC32c) {
    Log(InfoLogLevel::WARN_LEVEL, ioptions.info_log,
        "Silently converting format_version to 1 because checksum is "
        "non-default");
    // silently convert format_version to 1 to keep consistent with current
    // behavior
    sanitized_table_options.format_version = 1;
  }

  rep_ = new Rep(ioptions, sanitized_table_options, internal_comparator,
                 int_tbl_prop_collector_factories, column_family_id, file,
                 compression_type, compression_opts, compression_dict,
                 skip_filters, column_family_name);

  if (rep_->filter_block != nullptr) {
    rep_->filter_block->StartBlock(0);
  }
  if (table_options.block_cache_compressed.get() != nullptr) {
    BlockBasedTable::GenerateCachePrefix(
        table_options.block_cache_compressed.get(), file->writable_file(),
        &rep_->compressed_cache_key_prefix[0],
        &rep_->compressed_cache_key_prefix_size);
  }
}

BlockBasedTableBuilder::~BlockBasedTableBuilder() {
  assert(rep_->closed);  // Catch errors where caller forgot to call Finish()
  delete rep_;
}

void BlockBasedTableBuilder::Add(const Slice& key, const Slice& value) {
  Rep* r = rep_;
  assert(!r->closed);
  if (!ok()) return;

  r->learned_index_builder->OnKeyAdded(key, r->offset);
  
  size_t shared = 0;
  if(r->kv_data.size() > 0) {
    Slice last_key_piece(r->kv_data.back().first);
    shared = key.difference_offset(last_key_piece);
  }
  
  r->kv_data.push_back({key.ToString(), value.ToString()});

  r->props.num_entries++;
  r->props.raw_key_size += (key.size() - shared);
  r->props.raw_value_size += value.size();

  r->offset = r->props.raw_key_size + r->props.raw_value_size;
}

void BlockBasedTableBuilder::Flush() {
  Rep* r = rep_;
  assert(!r->closed);
  if (!ok()) return;
  if (r->data_block.empty()) return;
  WriteBlock(&r->data_block, &r->pending_handle, true /* is_data_block */);
  if (ok() && !r->table_options.skip_table_builder_flush) {
    r->status = r->file->Flush();
  }
  if (r->filter_block != nullptr) {
    r->filter_block->StartBlock(r->offset);
  }
  r->props.data_size = r->offset;
  ++r->props.num_data_blocks;
}

void BlockBasedTableBuilder::WriteBlock(BlockBuilder* block,
                                        BlockHandle* handle,
                                        bool is_data_block) {
  WriteBlock(block->Finish(), handle, is_data_block);
  block->Reset();
}

void BlockBasedTableBuilder::WriteBlock(const Slice& raw_block_contents,
                                        BlockHandle* handle,
                                        bool is_data_block) {
  // File format contains a sequence of blocks where each block has:
  //    block_data: uint8[n]
  //    type: uint8
  //    crc: uint32
  assert(ok());
  Rep* r = rep_;

  auto type = r->compression_type;
  Slice block_contents;
  bool abort_compression = false;

  StopWatchNano timer(r->ioptions.env,
    ShouldReportDetailedTime(r->ioptions.env, r->ioptions.statistics));

  if (raw_block_contents.size() < kCompressionSizeLimit) {
    Slice compression_dict;
    if (is_data_block && r->compression_dict && r->compression_dict->size()) {
      compression_dict = *r->compression_dict;
    }

    block_contents = CompressBlock(raw_block_contents, r->compression_opts,
                                   &type, r->table_options.format_version,
                                   compression_dict, &r->compressed_output);

    // Some of the compression algorithms are known to be unreliable. If
    // the verify_compression flag is set then try to de-compress the
    // compressed data and compare to the input.
    if (type != kNoCompression && r->table_options.verify_compression) {
      // Retrieve the uncompressed contents into a new buffer
      BlockContents contents;
      Status stat = UncompressBlockContentsForCompressionType(
          block_contents.data(), block_contents.size(), &contents,
          r->table_options.format_version, compression_dict, type,
          r->ioptions);

      if (stat.ok()) {
        bool compressed_ok = contents.data.compare(raw_block_contents) == 0;
        if (!compressed_ok) {
          // The result of the compression was invalid. abort.
          abort_compression = true;
          Log(InfoLogLevel::ERROR_LEVEL, r->ioptions.info_log,
              "Decompressed block did not match raw block");
          r->status =
              Status::Corruption("Decompressed block did not match raw block");
        }
      } else {
        // Decompression reported an error. abort.
        r->status = Status::Corruption("Could not decompress");
        abort_compression = true;
      }
    }
  } else {
    // Block is too big to be compressed.
    abort_compression = true;
  }

  // Abort compression if the block is too big, or did not pass
  // verification.
  if (abort_compression) {
    RecordTick(r->ioptions.statistics, NUMBER_BLOCK_NOT_COMPRESSED);
    type = kNoCompression;
    block_contents = raw_block_contents;
  } else if (type != kNoCompression &&
             ShouldReportDetailedTime(r->ioptions.env,
                                      r->ioptions.statistics)) {
    MeasureTime(r->ioptions.statistics, COMPRESSION_TIMES_NANOS,
                timer.ElapsedNanos());
    MeasureTime(r->ioptions.statistics, BYTES_COMPRESSED,
                raw_block_contents.size());
    RecordTick(r->ioptions.statistics, NUMBER_BLOCK_COMPRESSED);
  }

  WriteRawBlock(block_contents, type, handle);
  r->compressed_output.clear();
}

void BlockBasedTableBuilder::WriteRawBlock(const Slice& block_contents,
                                           CompressionType type,
                                           BlockHandle* handle) {
  Rep* r = rep_;
  StopWatch sw(r->ioptions.env, r->ioptions.statistics, WRITE_RAW_BLOCK_MICROS);
  handle->set_offset(r->offset);
  handle->set_size(block_contents.size());
  r->status = r->file->Append(block_contents);
  if (r->status.ok()) {
    char trailer[kBlockTrailerSize];
    trailer[0] = type;
    char* trailer_without_type = trailer + 1;
    switch (r->table_options.checksum) {
      case kNoChecksum:
        // we don't support no checksum yet
        assert(false);
        // intentional fallthrough
      case kCRC32c: {
        auto crc = crc32c::Value(block_contents.data(), block_contents.size());
        crc = crc32c::Extend(crc, trailer, 1);  // Extend to cover block type
        EncodeFixed32(trailer_without_type, crc32c::Mask(crc));
        break;
      }
      case kxxHash: {
        void* xxh = XXH32_init(0);
        XXH32_update(xxh, block_contents.data(),
                     static_cast<uint32_t>(block_contents.size()));
        XXH32_update(xxh, trailer, 1);  // Extend  to cover block type
        EncodeFixed32(trailer_without_type, XXH32_digest(xxh));
        break;
      }
    }

    r->status = r->file->Append(Slice(trailer, kBlockTrailerSize));
    if (r->status.ok()) {
      r->status = InsertBlockInCache(block_contents, type, handle);
    }
    if (r->status.ok()) {
      r->offset += block_contents.size() + kBlockTrailerSize;
    }
  }
}

Status BlockBasedTableBuilder::status() const {
  return rep_->status;
}

static void DeleteCachedBlock(const Slice& key, void* value) {
  Block* block = reinterpret_cast<Block*>(value);
  delete block;
}

//
// Make a copy of the block contents and insert into compressed block cache
//
Status BlockBasedTableBuilder::InsertBlockInCache(const Slice& block_contents,
                                                  const CompressionType type,
                                                  const BlockHandle* handle) {
  Rep* r = rep_;
  Cache* block_cache_compressed = r->table_options.block_cache_compressed.get();

  if (type != kNoCompression && block_cache_compressed != nullptr) {

    size_t size = block_contents.size();

    std::unique_ptr<char[]> ubuf(new char[size + 1]);
    memcpy(ubuf.get(), block_contents.data(), size);
    ubuf[size] = type;

    BlockContents results(std::move(ubuf), size, true, type);

    Block* block = new Block(std::move(results), kDisableGlobalSequenceNumber);

    // make cache key by appending the file offset to the cache prefix id
    char* end = EncodeVarint64(
                  r->compressed_cache_key_prefix +
                  r->compressed_cache_key_prefix_size,
                  handle->offset());
    Slice key(r->compressed_cache_key_prefix, static_cast<size_t>
              (end - r->compressed_cache_key_prefix));

    // Insert into compressed block cache.
    block_cache_compressed->Insert(key, block, block->usable_size(),
                                   &DeleteCachedBlock);

    // Invalidate OS cache.
    r->file->InvalidateCache(static_cast<size_t>(r->offset), size);
  }
  return Status::OK();
}

Status BlockBasedTableBuilder::Finish() {
  Rep* r = rep_;
  bool empty_data_block = r->data_block.empty();

  auto rss = dynamic_cast<LearnedIndexBuilder*>(r->learned_index_builder.get())->
                     FinalizeRSSModel(128, 4 * 1024);
  
  r->offset = 0;
  r->props.num_entries = 0;
  r->props.raw_key_size = 0;
  r->props.raw_value_size = 0;
 
  rss::RadixStringSplineNode* current_rss_node = nullptr;
  size_t current_spline_num = 0;
  size_t current_block_num = 0;
  size_t prev_kv_num = 0;
 
  for(uint32_t kv_num = 0; kv_num < r->kv_data.size(); kv_num++) {

    auto& kv_entry = r->kv_data[kv_num];
    Slice key(kv_entry.first);
    Slice value(kv_entry.second);
 
    ValueType value_type = ExtractValueType(key);
    if (IsValueType(value_type)) {
 
      if (r->props.num_entries > 0) {
          if (r->internal_comparator.Compare(key, Slice(r->last_key)) <= 0) {
          // We were about to insert keys out of order. Abort.
          Log(InfoLogLevel::ERROR_LEVEL, r->ioptions.info_log,
                          "Out-of-order key insertion into block based table");
          r->status =
              Status::Corruption("Out-of-order key insertion into table");
          return r->status;
        }
      }
 
      auto model_output = rss->find(ExtractUserKey(key));
 
      auto should_flush = (current_rss_node != nullptr) && (current_rss_node != model_output.rss_node || current_spline_num != model_output.spline);
      // std::cout << "Finish:current_rss_node:\t" << current_rss_node << "\tmodel_output.rss_node:\t" << model_output.rss_node << \
        "\tcurrent_spline_num:\t" << current_spline_num << "\tmodel_output.spline:\t" << model_output.spline << std::endl;

      if (should_flush) {
        assert(!r->data_block.empty());
        // std::cout << "Finish:should_flush BLOCK:\t" << current_block_num <<"\tsize:\t" << r->data_block.CurrentSizeEstimate() << "\tkv_num:\t" << (kv_num - prev_kv_num) << std::endl;
        prev_kv_num = kv_num;
        current_rss_node->AddBlockNumToSpline(current_spline_num, key, current_block_num);
        rss->add_to_accumulate_block_kv_num(kv_num);
 
        Flush();
        ++current_block_num;
 
        if (ok()) {
          r->index_builder->AddIndexEntry(&r->last_key, &key, r->pending_handle);
        }
      }
 
      current_rss_node = model_output.rss_node;
      current_spline_num = model_output.spline;
 
 
      if (r->filter_block != nullptr) {
        r->filter_block->Add(ExtractUserKey(key));
      }
 
      r->last_key.assign(key.data(), key.size());
      r->data_block.Add(key, value);
      r->props.num_entries++;
      r->props.raw_key_size += key.size();
      r->props.raw_value_size += value.size();
 
      r->index_builder->OnKeyAdded(key, r->offset + r->data_block.CurrentSizeEstimate());
 
      NotifyCollectTableCollectorsOnAdd(key, value, r->offset,
                                        r->table_properties_collectors,
                                        r->ioptions.info_log);
 
    } else if (value_type == kTypeRangeDeletion) {
      r->range_del_block.Add(key, value);
      ++r->props.num_entries;
      r->props.raw_key_size += key.size();
      r->props.raw_value_size += value.size();
      NotifyCollectTableCollectorsOnAdd(key, value, r->offset,
                                        r->table_properties_collectors,
                                        r->ioptions.info_log);
    } else {
      assert(false);
    }
  }
 
  current_rss_node->AddBlockNumToSpline(current_spline_num, r->kv_data.back().first, current_block_num);
  rss->add_to_accumulate_block_kv_num(r->kv_data.size());
  empty_data_block = r->data_block.empty();
  r->kv_data.clear();

  Flush();
  assert(!r->closed);
  r->closed = true;

  BlockHandle filter_block_handle, metaindex_block_handle, index_block_handle, learned_index_block_handle,
      compression_dict_block_handle, range_del_block_handle;
  // Write filter block
  if (ok() && r->filter_block != nullptr) {
    auto filter_contents = r->filter_block->Finish();
    r->props.filter_size = filter_contents.size();
    WriteRawBlock(filter_contents, kNoCompression, &filter_block_handle);
  }

  // To make sure properties block is able to keep the accurate size of index
  // block, we will finish writing all index entries here and flush them
  // to storage after metaindex block is written.
  if (ok() && !empty_data_block) {
    r->index_builder->AddIndexEntry(
        &r->last_key, nullptr /* no next data block */, r->pending_handle);
  }

  IndexBuilder::IndexBlocks learned_index_blocks;
  auto s = r->learned_index_builder->Finish(&learned_index_blocks);
  if (!s.ok()) {
    return s;
  }

  IndexBuilder::IndexBlocks index_blocks;
  s = r->index_builder->Finish(&index_blocks);
  if (!s.ok()) {
    return s;
  }

  // Write meta blocks and metaindex block with the following order.
  //    1. [meta block: filter]
  //    2. [meta block: properties]
  //    3. [meta block: compression dictionary]
  //    4. [meta block: range deletion tombstone]
  //    5. [metaindex block]
  // write meta blocks
  MetaIndexBuilder meta_index_builder;
  for (const auto& item : index_blocks.meta_blocks) {
    BlockHandle block_handle;
    WriteBlock(item.second, &block_handle, false /* is_data_block */);
    meta_index_builder.Add(item.first, block_handle);
  }

  if (ok()) {
    if (r->filter_block != nullptr) {
      // Add mapping from "<filter_block_prefix>.Name" to location
      // of filter data.
      std::string key;
      if (r->filter_block->IsBlockBased()) {
        key = BlockBasedTable::kFilterBlockPrefix;
      } else {
        key = BlockBasedTable::kFullFilterBlockPrefix;
      }
      key.append(r->table_options.filter_policy->Name());
      meta_index_builder.Add(key, filter_block_handle);
    }

    // Write properties and compression dictionary blocks.
    {
      PropertyBlockBuilder property_block_builder;
      r->props.column_family_id = r->column_family_id;
      r->props.column_family_name = r->column_family_name;
      r->props.filter_policy_name = r->table_options.filter_policy != nullptr ?
          r->table_options.filter_policy->Name() : "";
      r->props.index_size =
          r->index_builder->EstimatedSize() + kBlockTrailerSize;
      r->props.comparator_name = r->ioptions.user_comparator != nullptr
                                     ? r->ioptions.user_comparator->Name()
                                     : "nullptr";
      r->props.merge_operator_name = r->ioptions.merge_operator != nullptr
                                         ? r->ioptions.merge_operator->Name()
                                         : "nullptr";
      r->props.compression_name = CompressionTypeToString(r->compression_type);
      r->props.prefix_extractor_name =
          r->ioptions.prefix_extractor != nullptr
              ? r->ioptions.prefix_extractor->Name()
              : "nullptr";

      std::string property_collectors_names = "[";
      property_collectors_names = "[";
      for (size_t i = 0;
           i < r->ioptions.table_properties_collector_factories.size(); ++i) {
        if (i != 0) {
          property_collectors_names += ",";
        }
        property_collectors_names +=
            r->ioptions.table_properties_collector_factories[i]->Name();
      }
      property_collectors_names += "]";
      r->props.property_collectors_names = property_collectors_names;

      // Add basic properties
      property_block_builder.AddTableProperty(r->props);

      // Add use collected properties
      NotifyCollectTableCollectorsOnFinish(r->table_properties_collectors,
                                           r->ioptions.info_log,
                                           &property_block_builder);

      BlockHandle properties_block_handle;
      WriteRawBlock(
          property_block_builder.Finish(),
          kNoCompression,
          &properties_block_handle
      );
      meta_index_builder.Add(kPropertiesBlock, properties_block_handle);

      // Write compression dictionary block
      if (r->compression_dict && r->compression_dict->size()) {
        WriteRawBlock(*r->compression_dict, kNoCompression,
                      &compression_dict_block_handle);
        meta_index_builder.Add(kCompressionDictBlock,
                               compression_dict_block_handle);
      }
    }  // end of properties/compression dictionary block writing

    if (ok() && !r->range_del_block.empty()) {
      WriteRawBlock(r->range_del_block.Finish(), kNoCompression,
                    &range_del_block_handle);
      meta_index_builder.Add(kRangeDelBlock, range_del_block_handle);
    }  // range deletion tombstone meta block
  }    // meta blocks

  // Write index block
  if (ok()) {
    // flush the meta index block
    WriteRawBlock(meta_index_builder.Finish(), kNoCompression,
                  &metaindex_block_handle);
    WriteBlock(index_blocks.index_block_contents, &index_block_handle,
               false /* is_data_block */);
    WriteRawBlock(learned_index_blocks.index_block_contents, kNoCompression,
                  &learned_index_block_handle);
  }

  // Write footer
  if (ok()) {
    // No need to write out new footer if we're using default checksum.
    // We're writing legacy magic number because we want old versions of RocksDB
    // be able to read files generated with new release (just in case if
    // somebody wants to roll back after an upgrade)
    // TODO(icanadi) at some point in the future, when we're absolutely sure
    // nobody will roll back to RocksDB 2.x versions, retire the legacy magic
    // number and always write new table files with new magic number
    bool legacy = (r->table_options.format_version == 0);
    // this is guaranteed by BlockBasedTableBuilder's constructor
    assert(r->table_options.checksum == kCRC32c ||
           r->table_options.format_version != 0);
    Footer footer(legacy ? kLegacyBlockBasedTableMagicNumber
                         : kBlockBasedTableMagicNumber,
                  r->table_options.format_version);
    footer.set_metaindex_handle(metaindex_block_handle);
    footer.set_index_handle(index_block_handle);
    footer.set_learned_index_handle(learned_index_block_handle);
    footer.set_checksum(r->table_options.checksum);
    std::string footer_encoding;
    footer.EncodeTo(&footer_encoding);
    r->status = r->file->Append(footer_encoding);
    if (r->status.ok()) {
      r->offset += footer_encoding.size();
    }
  }

  return r->status;
}

void BlockBasedTableBuilder::Abandon() {
  Rep* r = rep_;
  assert(!r->closed);
  r->closed = true;
}

uint64_t BlockBasedTableBuilder::NumEntries() const {
  return rep_->props.num_entries;
}

uint64_t BlockBasedTableBuilder::FileSize() const {
  return rep_->offset;
}

bool BlockBasedTableBuilder::NeedCompact() const {
  for (const auto& collector : rep_->table_properties_collectors) {
    if (collector->NeedCompact()) {
      return true;
    }
  }
  return false;
}

TableProperties BlockBasedTableBuilder::GetTableProperties() const {
  TableProperties ret = rep_->props;
  for (const auto& collector : rep_->table_properties_collectors) {
    for (const auto& prop : collector->GetReadableProperties()) {
      ret.readable_properties.insert(prop);
    }
    collector->Finish(&ret.user_collected_properties);
  }
  return ret;
}

const std::string BlockBasedTable::kFilterBlockPrefix = "filter.";
const std::string BlockBasedTable::kFullFilterBlockPrefix = "fullfilter.";
}  // namespace rocksdb
