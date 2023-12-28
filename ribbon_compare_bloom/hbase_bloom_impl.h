//  Copyright (c) 2019-present, Facebook, Inc. All rights reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
//
// Implementation details of various Bloom filter implementations used in
// RocksDB. (DynamicBloom is in a separate file for now because it
// supports concurrent write.)

#pragma once
#include <stddef.h>
#include <stdint.h>

#include <cmath>
#include <random>

#include "port/port.h"  // for PREFETCH
#include "rocksdb/slice.h"
#include "util/hash.h"

#include <iostream>

#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace HBASE_BLOOM_NAMESPACE {
// 随机数 包含low和high
static int GetRandomWithBound(int low, int high) {
  // 创建随机数引擎
  std::random_device rd;  // 随机设备
  std::mt19937 mt(rd());  // 梅森旋转算法

  // 定义随机数分布
  std::uniform_int_distribution<int> dist(low,
                                          high);  // 生成1到50之间的随机整数

  // 生成随机数
  return dist(mt);
}

static char GetRandomReadableChar() {
  int i = GetRandomWithBound(0, 26 * 2 + 10);
  if (i < 26) return (char)('A' + i);
  i -= 26;

  if (i < 26) return (char)('a' + i);
  i -= 26;

  if (i < 10) return (char)('0' + i);
  i -= 10;

  assert(i == 0);
  return '_';
}

static uint8_t* GetFixedOrderedKey(int i, uint32_t rowLen) {
  std::string s = "";
  // 固定长度32的开头
  for (int j = rowLen - 1; j >= 0; --j) {
    if ((i & (1 << j)) == 0)
      s.append("a");
    else
      s.append("b");
  }
  // 为了简单这里改一下 随机生成 固定长度为10 的结尾
  for (int j = 0; j < 10; ++j) {
    s += GetRandomReadableChar();
  }

  uint8_t* res = new uint8_t[s.length()];
  for (int j = 0; j < (int)s.length(); ++j) {
    res[j] = s[j];
  }
  return res;
}

static uint8_t* GetRandomValue(int& vLen) {
  std::string s = "";
  for (int i = 0; i < GetRandomWithBound(1, 2000); ++i) {
    s += (char)(32 + GetRandomWithBound(0, 94));
    vLen++;
  }

  uint8_t* res = new uint8_t[s.length()];
  for (int j = 0; j < (int)s.length(); ++j) {
    res[j] = s[j];
  }
  return res;
}

class Cell {
  short rowLen;
  int rowOffset;
  uint8_t familyLen;
  int familyOffset;
  int qualifierLen;
  int qualifierOffset;
  int valueLen;

 public:
  const uint8_t* cellArray;
  uint8_t* value;
  Cell(uint64_t id, uint32_t rl) { // just for row filter test
    cellArray = GetFixedOrderedKey(id, rl);
    rowLen = rl;
    rowOffset = 0;
  }
  Cell(const std::string& arr, int rLen) : rowLen(rLen) { // test by standard key gen
    rowOffset = 0;
    const char* a = arr.c_str();
    uint8_t* tmp = new uint8_t[rLen]; // 为了让析构函数正常执行，并兼容其他的构造函数，这里我只能new一个出来了，虽然这一步实际没有意义
    for (int i = 0; i < rLen; ++i) {
      tmp[i] = a[i];
    }
    cellArray = tmp;
  }
  Cell(const uint8_t* arr, int rLen, int rOffset, int fLen, int fOffset, int qLen,
       int qOffset, int valueLen, uint8_t* value)
      : rowLen(rLen),
        rowOffset(rOffset),
        familyLen(fLen),
        familyOffset(fOffset),
        qualifierLen(qLen),
        qualifierOffset(qOffset),
        valueLen(valueLen),
        value(value) {
    cellArray = arr;
  }

  ~Cell() {
    delete[] cellArray;
  }

  const uint8_t* GetRowArray();
  int GetRowOffset();
  short GetRowLength() const;

  const uint8_t* GetFamilyArray();
  int GetFamilyOffset();
  uint8_t GetFamilyLength();

  const uint8_t* GetQualifierArray();
  int GetQualifierOffset();
  int GetQualifierLength();

  uint8_t* GetValueArray();
  int GetValueLength();
};

// 模版类
template <typename T>
class HashKey {
 public:
  T t;
  HashKey(T& t) : t(t) {}
  // 获取Hash值
  virtual uint8_t Get(int pos) = 0;
  // 获取Key的长度
  virtual int Length() = 0;
  // 虚函数必须要设置虚析构函数，因为要确保子类资源正确释放
  // 但不要设置纯虚析构函数，别给自己找麻烦，求求了，心好累
  virtual ~HashKey() {};
};

class CellHashKey : public HashKey<Cell> {
 public:
  CellHashKey(Cell& cell) : HashKey<Cell>(cell) {}
  ~CellHashKey() override {};
  uint8_t Get(int pos) override { return 0; };
  int Length() override { return 0; };
};

class RowBloomHashKey : public CellHashKey {
 public:
  RowBloomHashKey(Cell& cell) : CellHashKey(cell){};
  uint8_t Get(int pos) override { return t.cellArray[t.GetRowOffset() + pos]; }

  int Length() override { return t.GetRowLength(); }

  ~RowBloomHashKey(){};
};

class RowColBloomHashKey : public CellHashKey {
 public:
  RowColBloomHashKey(Cell& cell) : CellHashKey(cell){};
  uint8_t Get(int pos) override {
    // Todo
    return t.cellArray[t.GetRowOffset() + pos];
  }

  int Length() override {
    // Todo
    return t.GetRowLength();
  }

  ~RowColBloomHashKey(){};
};

// Bloom Type
enum BLOOM_TYPE { ROW, ROW_COL };

// TODO HBase中chunk添加到queue之后会压缩，也是需要考虑测试的
class HBaseBloomChunk {
 public:
  int maxKeys;
  int keyCount;
  long byteSize;
  int hashCount;
  BLOOM_TYPE bloomType;
  uint8_t* bloom;

  // Constructor
  HBaseBloomChunk(int maxKeys, long byteSize, int hashCount,
                  BLOOM_TYPE bloomType)
      : maxKeys(maxKeys),
        byteSize(byteSize),
        hashCount(hashCount),
        bloomType(bloomType) {
    keyCount = 0;
    bloom = new uint8_t[byteSize];
  }

  ~HBaseBloomChunk() { delete[] bloom; }

  void Add(Cell& cell);
  void SetHashLoc(int hash1, int hash2);
  void Set(long pos);

  // Temporary Query Implementation
  // Hash is uint64_t is assumed
  bool tmpQuery(Cell& cell);

  // for UT
  void ViewBloomUsage();

  // Hash Function
  int MurmurHash(CellHashKey& hashKey, int initVal);
  void Murmur3Hash(){/*todo*/};
  void JenkinsHash(){/*todo*/};
};

// static function for create HBaseBloomChunk
static unsigned int ComputeFoldableByteSize(long bitSize, int foldFactor) {
  long byteSizeLong = (bitSize + 7) / 8;
  int mask = (1 << foldFactor) - 1;
  if ((mask & byteSizeLong) != 0) {
    byteSizeLong >>= foldFactor;
    ++byteSizeLong;
    byteSizeLong <<= foldFactor;
  }
  return (int)byteSizeLong;
}

static long IdealMaxKeys(long bitSize, double errorRate) {
  return (long)(bitSize * (log(2) * log(2) / -log(errorRate)));
}

static int OptimalFunctionCount(int maxKeys, long bitSize) {
  long i = bitSize / maxKeys;
  double result = ceil(log(2) * i);
  return (int)result;
}

static long ComputeMaxKeys(long bitSize, double errorRate, int hashCount) {
  return (long)(-bitSize * 1.0 / hashCount *
                log(1 - exp(log(errorRate) / hashCount)));
}

static HBaseBloomChunk CreateBySize(int byteSizeHint, double errorRate,
                                    int foldFactor) {
  long byteSize = ComputeFoldableByteSize(byteSizeHint * 8L, foldFactor);
  long bitSize = byteSize * 8L;
  int maxKeys = IdealMaxKeys(bitSize, errorRate);
  int hashCount = OptimalFunctionCount(maxKeys, bitSize);
  maxKeys = (int)ComputeMaxKeys(bitSize, errorRate, hashCount);
  return HBaseBloomChunk(maxKeys, byteSize, hashCount, ROW);
}



}  // namespace HBASE_BLOOM_NAMESPACE
