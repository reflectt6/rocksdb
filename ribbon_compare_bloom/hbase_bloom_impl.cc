//  Copyright (c) 2019-present, Facebook, Inc. All rights reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
//
// Implementation details of various Bloom filter implementations used in
// RocksDB. (DynamicBloom is in a separate file for now because it
// supports concurrent write.)

#include "hbase_bloom_impl.h"

namespace HBASE_BLOOM_NAMESPACE {

// Cell
const uint8_t* Cell::GetRowArray() {
  return &cellArray[rowOffset];
}

int Cell::GetRowOffset() {
  return rowOffset;
}
short Cell::GetRowLength() const {
  return rowLen;
}

const uint8_t* Cell::GetFamilyArray() {
  return &cellArray[familyOffset];
}
int Cell::GetFamilyOffset() {
  return familyOffset;
}
uint8_t Cell::GetFamilyLength() {
    return familyLen;
}

const uint8_t* Cell::GetQualifierArray() {
    return &cellArray[qualifierOffset];
}
int Cell::GetQualifierOffset() {
    return qualifierOffset;
}
int Cell::GetQualifierLength() {
    return qualifierLen;
}

uint8_t* Cell::GetValueArray() {
    return value;
}
int Cell::GetValueLength() {
    return valueLen;
}

// HBaseBloomChunk
void HBaseBloomChunk::Add(Cell& cell) {
    int hash1;
    int hash2;
    CellHashKey* hashKey;
    if (bloomType == ROW) {
      hashKey = new RowBloomHashKey(cell);
    } else {
      throw std::exception();
    }
    hash1 = this->MurmurHash(*hashKey, 0);
    hash2 = this->MurmurHash(*hashKey, hash1);
    SetHashLoc(hash1, hash2);
    delete hashKey;
}

bool HBaseBloomChunk::tmpQuery(Cell& cell) {
    int hash1;
    int hash2;
    CellHashKey* hashKey;
    if (bloomType == ROW) {
      hashKey = new RowBloomHashKey(cell);
    } else {
      throw std::exception();
    }
    hash1 = this->MurmurHash(*hashKey, 0);
    hash2 = this->MurmurHash(*hashKey, hash1);

    for (int i = 0; i < hashCount; ++i) {
      long hashLoc = abs((hash1 + i * hash2) % (byteSize * 8));
      int bytePos = (int) hashLoc / 8;
      int bitPos = (int) hashLoc % 8;
      if ((bloom[bytePos] >> (bitPos) & 1) != 1) {
        delete hashKey;
        return false;
      }
    }
    delete hashKey;
    return true;
}

void HBaseBloomChunk::SetHashLoc(int hash1, int hash2) {
    for (int i = 0; i < hashCount; ++i) {
      long hashLoc = abs((hash1 + i * hash2) % (byteSize * 8));
      Set(hashLoc);
    }
    keyCount++;
}

void HBaseBloomChunk::Set(long pos) {
    int bytePos = (int) pos / 8;
    int bitPos = (int) pos % 8;
    bloom[bytePos] |= 1 << (bitPos);
}

void HBaseBloomChunk::ViewBloomUsage() {
    int usage = 0;
    for (int i = 0; i < byteSize; ++i) {
      for (int j = 0; j < 8; ++j) {
        uint8_t a = bloom[i];
        uint8_t b = (1 << j);
        uint8_t c = a & b;
        if (c != 0) {
          usage++;
        }
      }
    }
    double utilizationRate = usage * 0.1 / (byteSize * 8);
    std::cout << "utilizationRate = " << utilizationRate << std::endl;
}

// Hash Function
int HBaseBloomChunk::MurmurHash(CellHashKey& hashKey, int seed) {
    u_int32_t m = 0x5bd1e995;
    int r = 24;
    int length = hashKey.Length();
    u_int32_t h = seed ^ length;

    int len_4 = length >> 2;

    for (int i = 0; i < len_4; i++) {
      int i_4 = (i << 2);
      u_int32_t k = hashKey.Get(i_4 + 3);
      k = k << 8;
      k = k | (hashKey.Get(i_4 + 2) & 0xff);
      k = k << 8;
      k = k | (hashKey.Get(i_4 + 1) & 0xff);
      k = k << 8;
      // noinspection PointlessArithmeticExpression
      k = k | (hashKey.Get(i_4 + 0) & 0xff);
      k *= m;
      k ^= k >> r;
      k *= m;
      h *= m;
      h ^= k;
    }

    // avoid calculating modulo
    int len_m = len_4 << 2;
    int left = length - len_m;
    int i_m = len_m;

    if (left != 0) {
      if (left >= 3) {
        h ^= hashKey.Get(i_m + 2) << 16;
      }
      if (left >= 2) {
        h ^= hashKey.Get(i_m + 1) << 8;
      }
      if (left >= 1) {
        h ^= hashKey.Get(i_m);
      }

      h *= m;
    }

    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;

    return h;
}



}  // namespace HBASE_BLOOM_NAMESPACE
