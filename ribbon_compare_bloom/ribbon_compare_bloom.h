//
// Created by Rain Night on 2023/12/28.
//

#ifndef ROCKSDB_RIBBON_COMPILE_BLOOM_H
#define ROCKSDB_RIBBON_COMPILE_BLOOM_H

#include "ribbon_compare_bloom/hbase_bloom_impl.h"
#include "util/ribbon_test.h"

namespace RIBBON_COMPILE_BLOOM {

using namespace RIBBON_TEST;
using namespace HBASE_BLOOM_NAMESPACE;

using namespace ROCKSDB_NAMESPACE;
using namespace ROCKSDB_NAMESPACE::ribbon;

// HBase Bloom
DEFINE_uint32(chunkByteSizeHint, 128 * 1024, "io.storefile.bloom.block.size");
// Maximum number of times a Bloom filter can be "folded" if oversized
DEFINE_uint32(foldFactor, 7, "io.storefile.bloom.max.fold");
// DEFINE_double(errorRate, 0.01, "io.storefile.bloom.error.rate");
DEFINE_double(errorRate, 0.004, "io.storefile.bloom.error.rate");

// param for test
DEFINE_uint32(testRowLength, 32, "hbase row length for test");

// Generate hbase keys
struct HBaseKeyGen {
  HBaseKeyGen(uint64_t id) : id_(id) {}

  // Prefix (only one required)
  HBaseKeyGen& operator++() {
    ++id_;
    return *this;
  }

  HBaseKeyGen& operator+=(uint64_t i) {
    id_ += i;
    return *this;
  }

  const HBASE_BLOOM_NAMESPACE::Cell operator*() {
    HBASE_BLOOM_NAMESPACE::Cell* cell =
        new HBASE_BLOOM_NAMESPACE::Cell(id_, FLAGS_testRowLength);
    return *cell;
  }

  bool operator==(const HBaseKeyGen& other) const {
    // Same prefix is assumed
    return id_ == other.id_;
  }
  bool operator!=(const HBaseKeyGen& other) const {
    // Same prefix is assumed
    return id_ != other.id_;
  }

  uint64_t id_;
};

struct RibbonCompileBloomDefaultSettings : public DefaultTypesAndSettings {
  static constexpr bool kHomogeneous = false;
  static constexpr bool kUseSmash = false;
  using Key = Cell;
  using KeyGen = HBaseKeyGen;
  static Hash HashFn(const HBASE_BLOOM_NAMESPACE::Cell& key, uint64_t raw_seed) {
    // This version 0.7.2 preview of XXH3 (a.k.a. XXPH3) function does
    // not pass SmallKeyGen tests below without some seed premixing from
    // StandardHasher. See https://github.com/Cyan4973/xxHash/issues/469
    return ROCKSDB_NAMESPACE::Hash64(
        reinterpret_cast<const char*>(key.cellArray), key.GetRowLength(), raw_seed);
  }
};

struct Settings_Coeff128_Homog : public DefaultTypesAndSettings {
  // TODO 目前有两个问题需要解决
  //  1、在这个参数下运行cell测试，求解器会崩溃，每次都返回true，待定位
  //  2、cell运行贼慢，不知道为啥，不应该的
  static constexpr bool kHomogeneous = true;
  static constexpr bool kUseSmash = false;
  using Key = Cell;
  using KeyGen = HBaseKeyGen;
  static Hash HashFn(const HBASE_BLOOM_NAMESPACE::Cell& key, uint64_t raw_seed) {
    // This version 0.7.2 preview of XXH3 (a.k.a. XXPH3) function does
    // not pass SmallKeyGen tests below without some seed premixing from
    // StandardHasher. See https://github.com/Cyan4973/xxHash/issues/469
    return ROCKSDB_NAMESPACE::Hash64(
        reinterpret_cast<const char*>(key.cellArray), key.GetRowLength(), raw_seed);
  }
  static const std::vector<ConstructionFailureChance>& FailureChanceToTest() {
    return kFailureOnlyRare;
  }
};


struct Settings_Coeff128Smash_Homog : public Settings_Coeff128_Homog {
  static constexpr bool kUseSmash = true;
};

struct Settings_Coeff128_Homog2 : public DefaultTypesAndSettings {
  static constexpr bool kHomogeneous = true;
  static constexpr bool kUseSmash = false;
  using Seed = uint32_t;
  using Key = uint32_t;
  using KeyGen = Hash32KeyGenWrapper<StandardKeyGen>;
  static Hash HashFn(const uint32_t& key, uint32_t raw_seed) {
    // This version 0.7.2 preview of XXH3 (a.k.a. XXPH3) function does
    // not pass SmallKeyGen tests below without some seed premixing from
    // StandardHasher. See https://github.com/Cyan4973/xxHash/issues/469
    return ROCKSDB_NAMESPACE::Hash(
        reinterpret_cast<const char*>(&key), 1, raw_seed);
  }
  static const std::vector<ConstructionFailureChance>& FailureChanceToTest() {
    return kFailureOnlyRare;
  }
};

} // namespace

#endif  // ROCKSDB_RIBBON_COMPILE_BLOOM_H
