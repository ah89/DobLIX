//  Copyright (c) 2011-present, Facebook, Inc.  All rights reserved.
//  This source code is licensed under the BSD-style license found in the
//  LICENSE file in the root directory of this source tree. An additional grant
//  of patent rights can be found in the PATENTS file in the same directory.
//
// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#pragma once

#include <stdint.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "rocksdb/options.h"
#include "rocksdb/persistent_cache.h"
#include "rocksdb/statistics.h"
#include "rocksdb/status.h"
#include "rocksdb/table.h"
#include "table/table_properties_internal.h"
#include "table/table_reader.h"
#include "learned_index/radix_spline.h"
#include "util/cf_options.h"
#include "util/coding.h"
#include "util/file_reader_writer.h"

namespace rocksdb {

class Block;
class BlockIter;
class BlockHandle;
class Cache;
class FilterBlockReader;
class BlockBasedFilterBlockReader;
class FullFilterBlockReader;
class Footer;
class InternalKeyComparator;
class Iterator;
class RandomAccessFile;
class TableCache;
class TableReader;
class WritableFile;
struct BlockBasedTableOptions;
struct EnvOptions;
struct ReadOptions;
class GetContext;
class InternalIterator;

using std::unique_ptr;

typedef std::vector<std::pair<std::string, std::string>> KVPairBlock;

// A Table is a sorted map from strings to strings.  Tables are
// immutable and persistent.  A Table may be safely accessed from
// multiple threads without external synchronization.
class BlockBasedTable : public TableReader {
 public:
  static const std::string kFilterBlockPrefix;
  static const std::string kFullFilterBlockPrefix;
  // The longest prefix of the cache key used to identify blocks.
  // For Posix files the unique ID is three varints.
  static const size_t kMaxCacheKeyPrefixSize = kMaxVarint64Length * 3 + 1;

  // Attempt to open the table that is stored in bytes [0..file_size)
  // of "file", and read the metadata entries necessary to allow
  // retrieving data from the table.
  //
  // If successful, returns ok and sets "*table_reader" to the newly opened
  // table.  The client should delete "*table_reader" when no longer needed.
  // If there was an error while initializing the table, sets "*table_reader"
  // to nullptr and returns a non-ok status.
  //
  // @param file must remain live while this Table is in use.
  // @param prefetch_index_and_filter_in_cache can be used to disable
  // prefetching of
  //    index and filter blocks into block cache at startup
  // @param skip_filters Disables loading/accessing the filter block. Overrides
  //    prefetch_index_and_filter_in_cache, so filter will be skipped if both
  //    are set.
  static Status Open(const ImmutableCFOptions& ioptions,
                     const EnvOptions& env_options,
                     const BlockBasedTableOptions& table_options,
                     const InternalKeyComparator& internal_key_comparator,
                     unique_ptr<RandomAccessFileReader>&& file,
                     uint64_t file_size, unique_ptr<TableReader>* table_reader,
                     bool prefetch_index_and_filter_in_cache = true,
                     bool skip_filters = false, int level = -1);

  bool PrefixMayMatch(const Slice& internal_key);

  // Returns a new iterator over the table contents.
  // The result of NewIterator() is initially invalid (caller must
  // call one of the Seek methods on the iterator before using it).
  // @param skip_filters Disables loading/accessing the filter block
  InternalIterator* NewIterator(const ReadOptions&, Arena* arena = nullptr,
                                bool skip_filters = false) override;

  InternalIterator* NewRangeTombstoneIterator(
      const ReadOptions& read_options) override;

  // @param skip_filters Disables loading/accessing the filter block
  Status Get(const ReadOptions& readOptions, const Slice& key,
             GetContext* get_context, bool skip_filters = false) override;

  Status GetUsingLearnedIndex(const ReadOptions& readOptions, const Slice& key,
             GetContext* get_context, bool skip_filters = false);
 
  void GetKeyFromBlock(rocksdb::GetContext* get_context,
                        const rs::SplineModelResult& model_output,
                        rocksdb::FilterBlockReader* filter,
                        const rocksdb::Slice& key,
                        const rocksdb::ReadOptions& read_options,
                        PinnedIteratorsManager* pinned_iters_mgr,
                        bool pin_blocks,
                        rocksdb::Status& s);

  // Pre-fetch the disk blocks that correspond to the key range specified by
  // (kbegin, kend). The call will return error status in the event of
  // IO or iteration error.
  Status Prefetch(const Slice* begin, const Slice* end) override;

  // Given a key, return an approximate byte offset in the file where
  // the data for that key begins (or would begin if the key were
  // present in the file).  The returned value is in terms of file
  // bytes, and so includes effects like compression of the underlying data.
  // E.g., the approximate offset of the last key in the table will
  // be close to the file length.
  uint64_t ApproximateOffsetOf(const Slice& key) override;

  // Returns true if the block for the specified key is in cache.
  // REQUIRES: key is in this table && block cache enabled
  bool TEST_KeyInCache(const ReadOptions& options, const Slice& key);

  // Set up the table for Compaction. Might change some parameters with
  // posix_fadvise
  void SetupForCompaction() override;

  std::shared_ptr<const TableProperties> GetTableProperties() const override;

  size_t ApproximateMemoryUsage() const override;

  // convert SST file to a human readable form
  Status DumpTable(WritableFile* out_file) override;

  void Close() override;

  ~BlockBasedTable();

  bool TEST_filter_block_preloaded() const;
  bool TEST_index_reader_preloaded() const;
  // Implementation of IndexReader will be exposed to internal cc file only.
  class IndexReader;

  static Slice GetCacheKey(const char* cache_key_prefix,
                           size_t cache_key_prefix_size,
                           const BlockHandle& handle, char* cache_key);

  // Retrieve all key value pairs from data blocks in the table.
  // The key retrieved are internal keys.
  Status GetKVPairsFromDataBlocks(std::vector<KVPairBlock>* kv_pair_blocks);

 private:
  template <class TValue>
  struct CachableEntry;

  struct Rep;
  Rep* rep_;
  bool compaction_optimized_;

  class BlockEntryIteratorState;
  // input_iter: if it is not null, update this one and return it as Iterator
  static InternalIterator* NewDataBlockIterator(
      Rep* rep, const ReadOptions& ro, const Slice& index_value,
      BlockIter* input_iter = nullptr);

  static InternalIterator* NewDataBlockIterator(
    Rep* rep, const ReadOptions& ro, const BlockHandle& handle,
    BlockIter* input_iter);
  // If block cache enabled (compressed or uncompressed), looks for the block
  // identified by handle in (1) uncompressed cache, (2) compressed cache, and
  // then (3) file. If found, inserts into the cache(s) that were searched
  // unsuccessfully (e.g., if found in file, will add to both uncompressed and
  // compressed caches if they're enabled).
  //
  // @param block_entry value is set to the uncompressed block if found. If
  //    in uncompressed block cache, also sets cache_handle to reference that
  //    block.
  static Status MaybeLoadDataBlockToCache(
      Rep* rep, const ReadOptions& ro, const BlockHandle& handle,
      Slice compression_dict, CachableEntry<Block>* block_entry);

  // For the following two functions:
  // if `no_io == true`, we will not try to read filter/index from sst file
  // were they not present in cache yet.
  CachableEntry<FilterBlockReader> GetFilter(bool no_io = false) const;

  // Get the iterator from the index reader.
  // If input_iter is not set, return new Iterator
  // If input_iter is set, update it and return it as Iterator
  //
  // Note: ErrorIterator with Status::Incomplete shall be returned if all the
  // following conditions are met:
  //  1. We enabled table_options.cache_index_and_filter_blocks.
  //  2. index is not present in block cache.
  //  3. We disallowed any io to be performed, that is, read_options ==
  //     kBlockCacheTier
  InternalIterator* NewIndexIterator(
      const ReadOptions& read_options, BlockIter* input_iter = nullptr,
      CachableEntry<IndexReader>* index_entry = nullptr);

  // Read block cache from block caches (if set): block_cache and
  // block_cache_compressed.
  // On success, Status::OK with be returned and @block will be populated with
  // pointer to the block as well as its block handle.
  // @param compression_dict Data for presetting the compression library's
  //    dictionary.
  static Status GetDataBlockFromCache(
      const Slice& block_cache_key, const Slice& compressed_block_cache_key,
      Cache* block_cache, Cache* block_cache_compressed,
      const ImmutableCFOptions& ioptions, const ReadOptions& read_options,
      BlockBasedTable::CachableEntry<Block>* block, uint32_t format_version,
      const Slice& compression_dict, size_t read_amp_bytes_per_bit);

  // Put a raw block (maybe compressed) to the corresponding block caches.
  // This method will perform decompression against raw_block if needed and then
  // populate the block caches.
  // On success, Status::OK will be returned; also @block will be populated with
  // uncompressed block and its cache handle.
  //
  // REQUIRES: raw_block is heap-allocated. PutDataBlockToCache() will be
  // responsible for releasing its memory if error occurs.
  // @param compression_dict Data for presetting the compression library's
  //    dictionary.
  static Status PutDataBlockToCache(
      const Slice& block_cache_key, const Slice& compressed_block_cache_key,
      Cache* block_cache, Cache* block_cache_compressed,
      const ReadOptions& read_options, const ImmutableCFOptions& ioptions,
      CachableEntry<Block>* block, Block* raw_block, uint32_t format_version,
      const Slice& compression_dict, size_t read_amp_bytes_per_bit);

  // Calls (*handle_result)(arg, ...) repeatedly, starting with the entry found
  // after a call to Seek(key), until handle_result returns false.
  // May not make such a call if filter policy says that key is not present.
  friend class TableCache;
  friend class BlockBasedTableBuilder;

  void ReadMeta(const Footer& footer);

  // Create a index reader based on the index type stored in the table.
  // Optionally, user can pass a preloaded meta_index_iter for the index that
  // need to access extra meta blocks for index construction. This parameter
  // helps avoid re-reading meta index block if caller already created one.
  Status CreateIndexReader(
      IndexReader** index_reader,
      InternalIterator* preloaded_meta_index_iter = nullptr);

  bool FullFilterKeyMayMatch(const ReadOptions& read_options,
                             FilterBlockReader* filter,
                             const Slice& user_key) const;

  // Read the meta block from sst.
  static Status ReadMetaBlock(Rep* rep, std::unique_ptr<Block>* meta_block,
                              std::unique_ptr<InternalIterator>* iter);

  // Create the filter from the filter block.
  static FilterBlockReader* ReadFilter(Rep* rep);

  static void SetupCacheKeyPrefix(Rep* rep, uint64_t file_size);

  explicit BlockBasedTable(Rep* rep)
      : rep_(rep), compaction_optimized_(false) {}

  // Generate a cache key prefix from the file
  static void GenerateCachePrefix(Cache* cc,
    RandomAccessFile* file, char* buffer, size_t* size);
  static void GenerateCachePrefix(Cache* cc,
    WritableFile* file, char* buffer, size_t* size);

  // Helper functions for DumpTable()
  Status DumpIndexBlock(WritableFile* out_file);
  Status DumpDataBlocks(WritableFile* out_file);
  void DumpKeyValue(const Slice& key, const Slice& value,
                    WritableFile* out_file);

  // No copying allowed
  explicit BlockBasedTable(const TableReader&) = delete;
  void operator=(const TableReader&) = delete;
};

}  // namespace rocksdb
