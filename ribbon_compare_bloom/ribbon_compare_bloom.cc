//  Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#include "ribbon_compare_bloom.h"

using namespace RIBBON_COMPILE_BLOOM;

// 测试用例名称
class CompareTest : public ::testing::Test {};

TEST(CompareTest, TestHBaseBloomChunk) {

  HBaseBloomChunk chunk =
      CreateBySize(FLAGS_chunkByteSizeHint, FLAGS_errorRate, FLAGS_foldFactor);
  uint64_t duration = 0;
  ROCKSDB_NAMESPACE::StopWatchNano timer(ROCKSDB_NAMESPACE::SystemClock::Default().get());
  for (int i = 0; i < 109000; ++i) {
    uint8_t* k = GetFixedOrderedKey(i, FLAGS_testRowLength);
    int vLen = 0;
    uint8_t* v = GetRandomValue(vLen);
    int cfLen = 4;  // 简化
    Cell* cell = new Cell(k, 32, 0, cfLen, 32, 6, 32 + cfLen, vLen, v);

    timer.Start();
    chunk.Add(*cell);
    duration += timer.ElapsedNanos(true);
  }
  chunk.ViewBloomUsage();
  std::cout << "bloom chunk size = " << chunk.byteSize / 1000 << "kb" << std::endl;
  std::cout << "save keys = " << chunk.keyCount << std::endl;
  std::cout << "error rate = " << FLAGS_errorRate << std::endl;
  std::cout << "construction time = " << duration << std::endl;
}

TEST(CompareTest, TestRibbonFilter) {
  IMPORT_RIBBON_TYPES_AND_SETTINGS(RibbonCompileBloomDefaultSettings);
  IMPORT_RIBBON_IMPL_TYPES(RibbonCompileBloomDefaultSettings);
  using KeyGen = RibbonCompileBloomDefaultSettings::KeyGen;

  size_t bytes = 128 * 1024;
  std::unique_ptr<char[]> buf(new char[bytes]);
  InterleavedSoln isoln(buf.get(), bytes);
  SimpleSoln soln;
  Hasher hasher;
  Banding banding;

  KeyGen begin(0);
  KeyGen end(108928);
//  KeyGen end(108928);

  ROCKSDB_NAMESPACE::StopWatchNano timer(ROCKSDB_NAMESPACE::SystemClock::Default().get());

  timer.Start();
  ASSERT_TRUE(banding.ResetAndFindSeedToSolve(/*slots*/ Index{0} - 1, begin, end));
//  ASSERT_TRUE(banding.ResetAndFindSeedToSolve(/*slots*/ 108928 + 128 * 10, begin, end));
//  ASSERT_TRUE(banding.ResetAndFindSeedToSolve(/*slots*/ 217856, begin, end));
  auto duration = timer.ElapsedNanos(true);

  // get solution by backSubst
  timer.Start();
  isoln.BackSubstFrom(banding);
  auto duration2 = timer.ElapsedNanos(true);

  timer.Start();
  soln.BackSubstFrom(banding);
  auto duration3 = timer.ElapsedNanos(true);

  // test query
  for (int i = 0; i < 108928; ++i) {
    KeyGen tmp(i);
    ASSERT_TRUE(isoln.FilterQuery(*tmp, hasher));
    ASSERT_TRUE(soln.FilterQuery(*tmp, hasher));
  }
  int cnt1 = 0;
  int cnt2 = 0;
  for (int i = 108928; i < 110000; ++i) {
    KeyGen tmp(i);
    bool a = isoln.FilterQuery(*tmp, hasher);
    bool v = soln.FilterQuery(*tmp, hasher);
    if (a) cnt1++;
    if (v) cnt2++;
//    ASSERT_FALSE(isoln.FilterQuery(*tmp, hasher));
//    ASSERT_FALSE(soln.FilterQuery(*tmp, hasher));
    short c = 1;
  }

  // And report that in FP rate
  std::cout << "isoln ExpectedFpRate:" << isoln.ExpectedFpRate() << std::endl;
  std::cout << "soln ExpectedFpRate:" << soln.ExpectedFpRate() << std::endl;
  std::cout << "isoln错误个数:" << cnt1 << std::endl;
  std::cout << "soln错误个数:" << cnt2 << std::endl;
  std::cout << "错误率1:" << cnt1/ (110000 - 108928 + 0.1)<< std::endl;
  std::cout << "错误率2:" << cnt2/ (110000 - 108928 + 0.1)<< std::endl;
  std::cout << "banding时间:" << duration << std::endl;
  std::cout << "isoln时间:" << duration2<< std::endl;
  std::cout << "soln时间:" << duration3<< std::endl;

//  std::cout << "核心解矩阵大小为:" << 108928 / 1000 << "kb" << std::endl;
}

TEST(CompareTest, TestRibbonFilter2) {
  IMPORT_RIBBON_TYPES_AND_SETTINGS(RibbonCompileBloomDefaultSettings);
  IMPORT_RIBBON_IMPL_TYPES(RibbonCompileBloomDefaultSettings);
  using KeyGen = RibbonCompileBloomDefaultSettings::KeyGen;

  size_t bytes = 128 * 1024;
  std::unique_ptr<char[]> buf(new char[bytes]);
  InterleavedSoln isoln(buf.get(), bytes);
  SimpleSoln soln;
  Hasher hasher;
  Banding banding;

  KeyGen begin(50);
  KeyGen end(55);

  ASSERT_TRUE(banding.ResetAndFindSeedToSolve(/*slots*/2 * kCoeffBits , begin, end));

  // get solution by backSubst
  isoln.BackSubstFrom(banding);
  soln.BackSubstFrom(banding);

  // test query
  for (int i = 48; i < 52; ++i) {
    KeyGen tmp(i);
    bool a = isoln.FilterQuery(*tmp, hasher);
    bool  b = soln.FilterQuery(*tmp, hasher);
    short c = 1;
  }
  for (int i = 54; i < 56; ++i) {
    KeyGen tmp(i);
    bool a = isoln.FilterQuery(*tmp, hasher);
    bool  b = soln.FilterQuery(*tmp, hasher);
    short c = 1;
  }

  // And report that in FP rate
  std::cout << "isoln ExpectedFpRate:" << isoln.ExpectedFpRate() << std::endl;
  std::cout << "soln ExpectedFpRate:" << soln.ExpectedFpRate() << std::endl;
}

// for Repeat the test results in ribbon_test.cc 使用Cell 复刻原版UT
TEST(CompareTest, RepeatCompactnessAndBacktrackAndFpRate) {
  using TypeParam = Settings_Coeff128_Homog;

  ROCKSDB_NAMESPACE::StopWatchNano timeStub(ROCKSDB_NAMESPACE::SystemClock::Default().get());

  IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
  IMPORT_RIBBON_IMPL_TYPES(TypeParam);
  using KeyGen = TypeParam::KeyGen;
  using ConfigHelper = BandingConfigHelper<TypeParam>;

  if (sizeof(CoeffRow) < 8) {
    ROCKSDB_GTEST_BYPASS("Not fully supported");
    return;
  }

  const auto log2_thoroughness =
      static_cast<uint32_t>(ROCKSDB_NAMESPACE::FloorLog2(FLAGS_thoroughness));

  // 下面给出了中小规模filter的添加key数量范围，原理待研究
  // We are going to choose num_to_add using an exponential distribution,
  // so that we have good representation of small-to-medium filters.
  // Here we just pick some reasonable, practical upper bound based on
  // kCoeffBits or option. // log() 计算以e为底的对数
  const double log_max_add = std::log(
      FLAGS_max_add > 0 ? FLAGS_max_add
                        : static_cast<uint32_t>(kCoeffBits * kCoeffBits) *
                              std::max(FLAGS_thoroughness, uint32_t{32}));

  // This needs to be enough below the minimum number of slots to get a
  // reasonable number of samples with the minimum number of slots.
  const double log_min_add = std::log(0.66 * SimpleSoln::RoundUpNumSlots(1));

  ASSERT_GT(log_max_add, log_min_add);

  const double diff_log_add = log_max_add - log_min_add;

  for (ConstructionFailureChance cs : TypeParam::FailureChanceToTest()) {
    double expected_reseeds;
    switch (cs) {
      default:
        assert(false);
        FALLTHROUGH_INTENDED;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn2:
        fprintf(stderr, "== Failure: 50 percent\n");
        expected_reseeds = 1.0;
        break;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn20:
        fprintf(stderr, "== Failure: 5 percent\n"); //这里应该是写错了，应该是5%，写的是95%。。
        expected_reseeds = 0.053;
        break;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn1000:
        fprintf(stderr, "== Failure: 1/1000\n");
        expected_reseeds = 0.001;
        break;
    }

    uint64_t total_reseeds = 0;
    uint64_t total_singles = 0;
    uint64_t total_single_failures = 0;
    uint64_t total_batch = 0;
    uint64_t total_batch_successes = 0;
    uint64_t total_fp_count = 0;
    uint64_t total_added = 0;
    uint64_t total_expand_trials = 0;
    uint64_t total_expand_failures = 0;
    double total_expand_overhead = 0.0;

    uint64_t soln_query_nanos = 0; // soln查询时间
    uint64_t soln_query_count = 0; // soln查询数量
    uint64_t bloom_query_nanos = 0; // bloom查询的时间
    uint64_t isoln_query_nanos = 0; // isoln查询时间
    uint64_t isoln_query_count = 0; // isoln查询数量

    // Take different samples if you change thoroughness
    ROCKSDB_NAMESPACE::Random32 rnd(FLAGS_thoroughness);

    for (uint32_t i = 0; i < FLAGS_thoroughness; ++i) { // 不太懂这个迭代次数有啥用

      std::cout << std::endl;
      std::cout << "第" << i << "次迭代：" << std::endl;
      timeStub.Start();

      // We are going to choose num_to_add using an exponential distribution
      // as noted above, but instead of randomly choosing them, we generate
      // samples linearly using the golden ratio, which ensures a nice spread
      // even for a small number of samples, and starting with the minimum
      // number of slots to ensure it is tested. // fmod是浮点数取余
      double log_add =
          std::fmod(0.6180339887498948482 * diff_log_add * i, diff_log_add) +
          log_min_add; // 但是可以看到要添加的key的数量是收到迭代次数影响的
      uint32_t num_to_add = static_cast<uint32_t>(std::exp(log_add)); // e的x方

      // Most of the time, test the Interleaved solution storage, but when
      // we do we have to make num_slots a multiple of kCoeffBits. So
      // sometimes we want to test without that limitation.
      bool test_interleaved = (i % 7) != 6; // 迭代次数的作用1：
                                             // 正常来讲num_slots是kCoeffBits的倍数，但有时候我们希望没有这个限制，
                                             // 在迭代次数除7余6的时候，我们取消这个限制

      // Compute num_slots, and re-adjust num_to_add to get as close as possible
      // to next num_slots, to stress that num_slots in terms of construction
      // success. Ensure at least one iteration:
      Index num_slots = Index{0} - 1; // 得到一个Index类型所能表示的最大值
      --num_to_add;
      for (;;) { // 很迷惑，反正通过一些操作，把上面的num_slots重新给了一个较小的值 具体实现在ribbon_config.cc
        Index next_num_slots = SimpleSoln::RoundUpNumSlots(
            ConfigHelper::GetNumSlots(num_to_add + 1, cs));
        if (test_interleaved) {
          next_num_slots = InterleavedSoln::RoundUpNumSlots(next_num_slots);
          // assert idempotent
          EXPECT_EQ(next_num_slots,
                    InterleavedSoln::RoundUpNumSlots(next_num_slots));
        }
        // assert idempotent with InterleavedSoln::RoundUpNumSlots
        EXPECT_EQ(next_num_slots, SimpleSoln::RoundUpNumSlots(next_num_slots));

        if (next_num_slots > num_slots) {
          break;
        }
        num_slots = next_num_slots;
        ++num_to_add;
      }
      assert(num_slots < Index{0} - 1);

      total_added += num_to_add;

      std::string prefix;
      ROCKSDB_NAMESPACE::PutFixed32(&prefix, rnd.Next());

      std::cout << "计算num时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms， 其中 num_to_add = " << num_to_add << ", num_slots = " << num_slots << std::endl;

      // Batch that must be added
      uint32_t addedKeyPrefix = 100000000;
      KeyGen keys_begin(addedKeyPrefix + 0);
      KeyGen keys_end(addedKeyPrefix + num_to_add);
//      std::string added_str = prefix + "added";
//      KeyGen keys_begin(added_str, 0);
//      KeyGen keys_end(added_str, num_to_add);

      // A couple more that will probably be added
      uint32_t moreKeyPrefix = 200000000;
      KeyGen one_more(moreKeyPrefix + 1);
      KeyGen two_more(moreKeyPrefix + 2);
//      KeyGen one_more(prefix + "more", 1);
//      KeyGen two_more(prefix + "more", 2);

      // Batch that may or may not be added
      uint32_t batch_size =
          static_cast<uint32_t>(2.0 * std::sqrt(num_slots - num_to_add));
      if (batch_size < 10U) {
        batch_size = 0;
      }
      uint32_t batchKeyPrefix = 300000000;
      KeyGen batch_begin(batchKeyPrefix + 0);
      KeyGen batch_end(batchKeyPrefix + batch_size);
//      std::string batch_str = prefix + "batch";
//      KeyGen batch_begin(batch_str, 0);
//      KeyGen batch_end(batch_str, batch_size);

      // Batch never (successfully) added, but used for querying FP rate
      uint32_t notKeyPrefix = 400000000;
      KeyGen other_keys_begin(notKeyPrefix + 0);
      KeyGen other_keys_end(notKeyPrefix + FLAGS_max_check);

//      std::string not_str = prefix + "not";
//      KeyGen other_keys_begin(not_str, 0);
//      KeyGen other_keys_end(not_str, FLAGS_max_check);

      double overhead_ratio = 1.0 * num_slots / num_to_add;
      if (FLAGS_verbose) {
        fprintf(stderr, "Adding(%s) %u / %u   Overhead: %g   Batch size: %u\n",
                test_interleaved ? "i" : "s", (unsigned)num_to_add,
                (unsigned)num_slots, overhead_ratio, (unsigned)batch_size);
      }

      // Vary bytes for InterleavedSoln to use number of solution columns
      // from 0 to max allowed by ResultRow type (and used by SimpleSoln).
      // Specifically include 0 and max, and otherwise skew toward max.
      uint32_t max_ibytes =
          static_cast<uint32_t>(sizeof(ResultRow) * num_slots);
      size_t ibytes;
      // 迭代轮数的作用2：设置不同大小的InterleavedSoln
      if (i == 0) {
        ibytes = 0;
      } else if (i == 1) {
        ibytes = max_ibytes;
      } else {
        // Skewed 偏向于更大的随机数
        ibytes =
            std::max(rnd.Uniformish(max_ibytes), rnd.Uniformish(max_ibytes));
      }
      std::unique_ptr<char[]> idata(new char[ibytes]);
      InterleavedSoln isoln(idata.get(), ibytes);

      std::cout << "生成HBaseKey 定义soln时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms， 其中 overhead_ratio = " << overhead_ratio << std::endl;
      SimpleSoln soln;
      Hasher hasher;
      bool first_single;
      bool second_single;
      bool batch_success;
      {
        Banding banding;
        // Traditional solve for a fixed set.
        ASSERT_TRUE(
            banding.ResetAndFindSeedToSolve(num_slots, keys_begin, keys_end));

        //        banding.printCoeRow((int )num_slots);
        //        banding.printResultRow((int )num_slots);

        Index occupied_count = banding.GetOccupiedCount();
        Index more_added = 0;

        std::cout << "banding的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms， 此时添加了num_to_add个key：" << num_to_add << "banding占据数量" << occupied_count<< std::endl;
        if (TypeParam::kHomogeneous || overhead_ratio < 1.01 ||
            batch_size == 0) {
          // 可以看出来backtracking是用来记录banding失败的
          // Homogeneous not compatible with backtracking because add
          // 在开销小于1.01时，我们不应该要求更多
          // doesn't fail. Small overhead ratio too packed to expect more
          first_single = false;
          second_single = false;
          batch_success = false;
        } else {
          // Now to test backtracking, starting with guaranteed fail. By using
          // the keys that will be used to test FP rate, we are then doing an
          // extra check that after backtracking there are no remnants (e.g. in
          // result side of banding) of these entries.
          KeyGen other_keys_too_big_end = other_keys_begin;
          other_keys_too_big_end += num_to_add; // 一般这个num_to_add是大于nums_slots的一半的，所以这里再添加num_to_add一定会超过filter容量
          banding.EnsureBacktrackSize(std::max(num_to_add, batch_size));
          EXPECT_FALSE(banding.AddRangeOrRollBack(other_keys_begin,
                                                  other_keys_too_big_end));// 由于超过容量所以失败
          EXPECT_EQ(occupied_count, banding.GetOccupiedCount());

          // Check that we still have a good chance of adding a couple more
          // individually
          first_single = banding.Add(*one_more);
          second_single = banding.Add(*two_more);
          more_added += (first_single ? 1 : 0) + (second_single ? 1 : 0);
          total_singles += 2U;
          total_single_failures += 2U - more_added;

          // Or as a batch
          batch_success = banding.AddRangeOrRollBack(batch_begin, batch_end);
          ++total_batch;
          if (batch_success) {
            more_added += batch_size;
            ++total_batch_successes;
          }
          EXPECT_LE(banding.GetOccupiedCount(), occupied_count + more_added);
        }

        std::cout << "测试backtracking的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms，是否被跳过："<<(TypeParam::kHomogeneous || overhead_ratio < 1.01 ||
            batch_size == 0 )<< std::endl;

        // Also verify that redundant adds are OK (no effect)
        ASSERT_TRUE(
            banding.AddRange(keys_begin, KeyGen(addedKeyPrefix + num_to_add / 8)));
//        ASSERT_TRUE(
//            banding.AddRange(keys_begin, KeyGen(added_str, num_to_add / 8)));
        EXPECT_LE(banding.GetOccupiedCount(), occupied_count + more_added);

        std::cout << "测试添加冗余key的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms" << std::endl;

        // Now back-substitution
        soln.BackSubstFrom(banding);
        if (test_interleaved) {
          isoln.BackSubstFrom(banding);
        }

        std::cout << "soln求解的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms，是否跳过isoln："<<!test_interleaved<< std::endl;

        Seed reseeds = banding.GetOrdinalSeed();
        total_reseeds += reseeds;

        EXPECT_LE(reseeds, 8 + log2_thoroughness); // 迷惑，其实干嘛要做log运算我是很迷惑的
        if (reseeds > log2_thoroughness + 1) {
          fprintf(
              stderr, "%s high reseeds at %u, %u/%u: %u\n",
              reseeds > log2_thoroughness + 8 ? "ERROR Extremely" : "Somewhat",
              static_cast<unsigned>(i), static_cast<unsigned>(num_to_add),
              static_cast<unsigned>(num_slots), static_cast<unsigned>(reseeds));
        }

        if (reseeds > 0) {
          // "Expand" test: given a failed construction, how likely is it to
          // pass with same seed and more slots. At each step, we increase
          // enough to ensure there is at least one shift within each coeff
          // block.
          ++total_expand_trials;
          Index expand_count = 0;
          Index ex_slots = num_slots;
          banding.SetOrdinalSeed(0);
          for (;; ++expand_count) {
            ASSERT_LE(expand_count, log2_thoroughness);
            ex_slots += ex_slots / kCoeffBits;
            if (test_interleaved) {
              ex_slots = InterleavedSoln::RoundUpNumSlots(ex_slots);
            }
            banding.Reset(ex_slots);
            bool success = banding.AddRange(keys_begin, keys_end);
            if (success) {
              break;
            }
          }
          total_expand_failures += expand_count;
          total_expand_overhead += 1.0 * (ex_slots - num_slots) / num_slots;
        }

        hasher.SetOrdinalSeed(reseeds);
        std::cout << "测试通过小幅度增加空间，是否可以增加banding的成功率的时间（expand）：" << timeStub.ElapsedNanos(true)/1000/1000
                  << "ms，是否被跳过："<<!(reseeds > 0)
                  << ", seed = "<< reseeds<< std::endl;
      }
      // soln and hasher now independent of Banding object

      // Verify keys added
      KeyGen cur = keys_begin;
      while (cur != keys_end) {
        ASSERT_TRUE(soln.FilterQuery(*cur, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*cur, hasher));
        ++cur;
      }
      // We (maybe) snuck these in!
      if (first_single) {
        ASSERT_TRUE(soln.FilterQuery(*one_more, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*one_more, hasher));
      }
      if (second_single) {
        ASSERT_TRUE(soln.FilterQuery(*two_more, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*two_more, hasher));
      }
      if (batch_success) {
        cur = batch_begin;
        while (cur != batch_end) {
          ASSERT_TRUE(soln.FilterQuery(*cur, hasher));
          ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*cur, hasher));
          ++cur;
        }
      }

      std::cout << "查询已添加的key时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms，是否跳过isoln："<<!(test_interleaved) << std::endl;

      // Check FP rate (depends only on number of result bits == solution
      // columns)
      Index fp_count = 0;
      cur = other_keys_begin;
      {
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          bool fp = soln.FilterQuery(*cur, hasher);
          fp_count += fp ? 1 : 0;
          ++cur;
        }
        soln_query_nanos += timer.ElapsedNanos();
        soln_query_count += FLAGS_max_check;
      }
      {
        double expected_fp_count = soln.ExpectedFpRate() * FLAGS_max_check;
        // For expected FP rate, also include false positives due to collisions
        // in Hash value. (Negligible for 64-bit, can matter for 32-bit.)
        double correction =
            FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

        // NOTE: rare violations expected with kHomogeneous
        EXPECT_LE(fp_count,
                  FrequentPoissonUpperBound(expected_fp_count + correction)); // 泊松分布对比误差
        EXPECT_GE(fp_count,
                  FrequentPoissonLowerBound(expected_fp_count + correction));
      }
      total_fp_count += fp_count;

      std::cout << "查询soln 未添加的key，并计算fp 时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms"<< std::endl;

      // And also check FP rate for isoln
      if (test_interleaved) {
        Index ifp_count = 0;
        cur = other_keys_begin;
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          ifp_count += isoln.FilterQuery(*cur, hasher) ? 1 : 0;
          ++cur;
        }
        isoln_query_nanos += timer.ElapsedNanos();
        isoln_query_count += FLAGS_max_check;
        {
          double expected_fp_count = isoln.ExpectedFpRate() * FLAGS_max_check;
          // For expected FP rate, also include false positives due to
          // collisions in Hash value. (Negligible for 64-bit, can matter for
          // 32-bit.)
          double correction =
              FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

          // NOTE: rare violations expected with kHomogeneous
          EXPECT_LE(ifp_count,
                    FrequentPoissonUpperBound(expected_fp_count + correction));

          // FIXME: why sometimes can we slightly "beat the odds"?
          // 这个测试有时会不通过，因此引入了一个因子来降低期望的误差数，但我感觉这是个好事情？
          // (0.95 factor should not be needed)
          EXPECT_GE(ifp_count, FrequentPoissonLowerBound(
                                   0.95 * expected_fp_count + correction));
        }
        // Since the bits used in isoln are a subset of the bits used in soln,
        // it cannot have fewer FPs
        // 从这里看出来 isoln对比soln的优势不在fp，在fp上应该是持平的，所以可能在内存上有优化
        EXPECT_GE(ifp_count, fp_count);
      }

      std::cout << "查询isoln 未添加的key，并计算fp 时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms， 是否跳过"<< !test_interleaved<< std::endl;

      // And compare to Bloom time, for fun
      if (ibytes >= /* minimum Bloom impl bytes*/ 64) {
        Index bfp_count = 0;
        cur = other_keys_begin;
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          uint64_t h = hasher.GetHash(*cur);
          uint32_t h1 = ROCKSDB_NAMESPACE::Lower32of64(h);
          uint32_t h2 = sizeof(Hash) >= 8 ? ROCKSDB_NAMESPACE::Upper32of64(h)
                                          : h1 * 0x9e3779b9;
          bfp_count +=
              ROCKSDB_NAMESPACE::FastLocalBloomImpl::HashMayMatch(
                  h1, h2, static_cast<uint32_t>(ibytes), 6, idata.get())
                  ? 1
                  : 0;
          ++cur;
        }
        bloom_query_nanos += timer.ElapsedNanos();
        // ensure bfp_count is used
        ASSERT_LT(bfp_count, FLAGS_max_check);
      }

      std::cout << "对比bloom的时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms， 此时第" << i <<  "次迭代结束"<< std::endl;
    }

    // "outside" == key not in original set so either negative or false positive
    fprintf(stderr,
            "Simple      outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * soln_query_nanos / soln_query_count);
    fprintf(stderr,
            "Interleaved outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * isoln_query_nanos / isoln_query_count);
    fprintf(stderr,
            "Bloom       outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * bloom_query_nanos / soln_query_count);

    if (TypeParam::kHomogeneous) {
      EXPECT_EQ(total_reseeds, 0U); // 说明kHomogeneous在banding中不会失败，因此不会改变种子的值
    } else {
      double average_reseeds = 1.0 * total_reseeds / FLAGS_thoroughness;
      fprintf(stderr, "Average re-seeds: %g\n", average_reseeds);
      // Values above were chosen to target around 50% chance of encoding
      // success rate (average of 1.0 re-seeds) or slightly better. But 1.15 is
      // also close enough.
      EXPECT_LE(total_reseeds,
                InfrequentPoissonUpperBound(1.15 * expected_reseeds *
                                            FLAGS_thoroughness));
      // Would use 0.85 here instead of 0.75, but
      // TypesAndSettings_Hash32_SmallKeyGen can "beat the odds" because of
      // sequential keys with a small, cheap hash function. We accept that
      // there are surely inputs that are somewhat bad for this setup, but
      // these somewhat good inputs are probably more likely.
      EXPECT_GE(total_reseeds,
                InfrequentPoissonLowerBound(0.75 * expected_reseeds *
                                            FLAGS_thoroughness));
    }

    std::cout << "测试 reseed 是否符合预期的时间（测试逻辑不太懂）：" << timeStub.ElapsedNanos(true)/1000/1000
              << "ms"<< std::endl;

    if (total_expand_trials > 0) {
      double average_expand_failures =
          1.0 * total_expand_failures / total_expand_trials;
      fprintf(stderr, "Average expand failures, and overhead: %g, %g\n",
              average_expand_failures,
              total_expand_overhead / total_expand_trials);
      // Seems to be a generous allowance
      EXPECT_LE(total_expand_failures,
                InfrequentPoissonUpperBound(1.0 * total_expand_trials));
    } else {
      fprintf(stderr, "Average expand failures: N/A\n");
    }

    if (total_singles > 0) {
      double single_failure_rate = 1.0 * total_single_failures / total_singles;
      fprintf(stderr, "Add'l single, failure rate: %g\n", single_failure_rate);
      // A rough bound (one sided) based on nothing in particular
      double expected_single_failures = 1.0 * total_singles /
                                        (sizeof(CoeffRow) == 16 ? 128
                                         : TypeParam::kUseSmash ? 64
                                                                : 32);
      EXPECT_LE(total_single_failures,
                InfrequentPoissonUpperBound(expected_single_failures));
    }

    if (total_batch > 0) {
      // Counting successes here for Poisson to approximate the Binomial
      // distribution.
      // A rough bound (one sided) based on nothing in particular.
      double expected_batch_successes = 1.0 * total_batch / 2;
      uint64_t lower_bound =
          InfrequentPoissonLowerBound(expected_batch_successes);
      fprintf(stderr, "Add'l batch, success rate: %g (>= %g)\n",
              1.0 * total_batch_successes / total_batch,
              1.0 * lower_bound / total_batch);
      EXPECT_GE(total_batch_successes, lower_bound);
    }

    {
      uint64_t total_checked = uint64_t{FLAGS_max_check} * FLAGS_thoroughness;
      double expected_total_fp_count =
          total_checked * std::pow(0.5, 8U * sizeof(ResultRow));
      // For expected FP rate, also include false positives due to collisions
      // in Hash value. (Negligible for 64-bit, can matter for 32-bit.)
      double average_added = 1.0 * total_added / FLAGS_thoroughness;
      expected_total_fp_count +=
          total_checked * ExpectedCollisionFpRate(Hasher(), average_added);

      uint64_t upper_bound =
          InfrequentPoissonUpperBound(expected_total_fp_count);
      uint64_t lower_bound =
          InfrequentPoissonLowerBound(expected_total_fp_count);
      fprintf(stderr, "Average FP rate: %g (~= %g, <= %g, >= %g)\n",
              1.0 * total_fp_count / total_checked,
              expected_total_fp_count / total_checked,
              1.0 * upper_bound / total_checked,
              1.0 * lower_bound / total_checked);
      EXPECT_LE(total_fp_count, upper_bound);
      EXPECT_GE(total_fp_count, lower_bound);
    }

    std::cout << "计算expand、single、batch、avg fp rate是否符合预期的时间（应该没啥计算量） ：" << timeStub.ElapsedNanos(true)/1000/1000
              << "ms， 此时enum["<< cs << "]结束" << std::endl;
  }
}

// 尝试原版测试，但是使用32位的Hash
TEST(CompareTest, RepeatCompactnessAndBacktrackAndFpRate2) {
  using TypeParam = Settings_Coeff128_Homog2;

  ROCKSDB_NAMESPACE::StopWatchNano timeStub(ROCKSDB_NAMESPACE::SystemClock::Default().get());

  IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
  IMPORT_RIBBON_IMPL_TYPES(TypeParam);
  using KeyGen = TypeParam::KeyGen;
  using ConfigHelper = BandingConfigHelper<TypeParam>;

  if (sizeof(CoeffRow) < 8) {
    ROCKSDB_GTEST_BYPASS("Not fully supported");
    return;
  }

  const auto log2_thoroughness =
      static_cast<uint32_t>(ROCKSDB_NAMESPACE::FloorLog2(FLAGS_thoroughness));

  // 下面给出了中小规模filter的添加key数量范围，原理待研究
  // We are going to choose num_to_add using an exponential distribution,
  // so that we have good representation of small-to-medium filters.
  // Here we just pick some reasonable, practical upper bound based on
  // kCoeffBits or option. // log() 计算以e为底的对数
  const double log_max_add = std::log(
      FLAGS_max_add > 0 ? FLAGS_max_add
                        : static_cast<uint32_t>(kCoeffBits * kCoeffBits) *
                              std::max(FLAGS_thoroughness, uint32_t{32}));

  // This needs to be enough below the minimum number of slots to get a
  // reasonable number of samples with the minimum number of slots.
  const double log_min_add = std::log(0.66 * SimpleSoln::RoundUpNumSlots(1));

  ASSERT_GT(log_max_add, log_min_add);

  const double diff_log_add = log_max_add - log_min_add;

  for (ConstructionFailureChance cs : TypeParam::FailureChanceToTest()) {
    double expected_reseeds;
    switch (cs) {
      default:
        assert(false);
        FALLTHROUGH_INTENDED;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn2:
        fprintf(stderr, "== Failure: 50 percent\n");
        expected_reseeds = 1.0;
        break;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn20:
        fprintf(stderr, "== Failure: 5 percent\n"); //这里应该是写错了，应该是5%，写的是95%。。
        expected_reseeds = 0.053;
        break;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn1000:
        fprintf(stderr, "== Failure: 1/1000\n");
        expected_reseeds = 0.001;
        break;
    }

    uint64_t total_reseeds = 0;
    uint64_t total_singles = 0;
    uint64_t total_single_failures = 0;
    uint64_t total_batch = 0;
    uint64_t total_batch_successes = 0;
    uint64_t total_fp_count = 0;
    uint64_t total_added = 0;
    uint64_t total_expand_trials = 0;
    uint64_t total_expand_failures = 0;
    double total_expand_overhead = 0.0;

    uint64_t soln_query_nanos = 0; // soln查询时间
    uint64_t soln_query_count = 0; // soln查询数量
    uint64_t bloom_query_nanos = 0; // bloom查询的时间
    uint64_t isoln_query_nanos = 0; // isoln查询时间
    uint64_t isoln_query_count = 0; // isoln查询数量

    // Take different samples if you change thoroughness
    ROCKSDB_NAMESPACE::Random32 rnd(FLAGS_thoroughness);

    for (uint32_t i = 0; i < FLAGS_thoroughness; ++i) { // 不太懂这个迭代次数有啥用

      std::cout << std::endl;
      std::cout << "第" << i << "次迭代：" << std::endl;
      timeStub.Start();

      // We are going to choose num_to_add using an exponential distribution
      // as noted above, but instead of randomly choosing them, we generate
      // samples linearly using the golden ratio, which ensures a nice spread
      // even for a small number of samples, and starting with the minimum
      // number of slots to ensure it is tested. // fmod是浮点数取余
      double log_add =
          std::fmod(0.6180339887498948482 * diff_log_add * i, diff_log_add) +
          log_min_add; // 但是可以看到要添加的key的数量是收到迭代次数影响的
      uint32_t num_to_add = static_cast<uint32_t>(std::exp(log_add)); // e的x方

      // Most of the time, test the Interleaved solution storage, but when
      // we do we have to make num_slots a multiple of kCoeffBits. So
      // sometimes we want to test without that limitation.
      bool test_interleaved = (i % 7) != 6; // 迭代次数的作用1：
                                             // 正常来讲num_slots是kCoeffBits的倍数，但有时候我们希望没有这个限制，
                                             // 在迭代次数除7余6的时候，我们取消这个限制

      // Compute num_slots, and re-adjust num_to_add to get as close as possible
      // to next num_slots, to stress that num_slots in terms of construction
      // success. Ensure at least one iteration:
      Index num_slots = Index{0} - 1; // 得到一个Index类型所能表示的最大值
      --num_to_add;
      for (;;) { // 很迷惑，反正通过一些操作，把上面的num_slots重新给了一个较小的值 具体实现在ribbon_config.cc
        Index next_num_slots = SimpleSoln::RoundUpNumSlots(
            ConfigHelper::GetNumSlots(num_to_add + 1, cs));
        if (test_interleaved) {
          next_num_slots = InterleavedSoln::RoundUpNumSlots(next_num_slots);
          // assert idempotent
          EXPECT_EQ(next_num_slots,
                    InterleavedSoln::RoundUpNumSlots(next_num_slots));
        }
        // assert idempotent with InterleavedSoln::RoundUpNumSlots
        EXPECT_EQ(next_num_slots, SimpleSoln::RoundUpNumSlots(next_num_slots));

        if (next_num_slots > num_slots) {
          break;
        }
        num_slots = next_num_slots;
        ++num_to_add;
      }
      assert(num_slots < Index{0} - 1);

      total_added += num_to_add;

      std::string prefix;
      ROCKSDB_NAMESPACE::PutFixed32(&prefix, rnd.Next());

      std::cout << "计算num时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms， 其中 num_to_add = " << num_to_add << ", num_slots = " << num_slots << std::endl;

      // Batch that must be added
      //      uint32_t addedKeyPrefix = 100000000;
      //      KeyGen keys_begin(addedKeyPrefix + 0);
      //      KeyGen keys_end(addedKeyPrefix + num_to_add);
      std::string added_str = prefix + "added";
      KeyGen keys_begin(added_str, 0);
      KeyGen keys_end(added_str, num_to_add);

      // A couple more that will probably be added
      //      uint32_t moreKeyPrefix = 200000000;
      //      KeyGen one_more(moreKeyPrefix + 1);
      //      KeyGen two_more(moreKeyPrefix + 2);
      KeyGen one_more(prefix + "more", 1);
      KeyGen two_more(prefix + "more", 2);

      // Batch that may or may not be added
      uint32_t batch_size =
          static_cast<uint32_t>(2.0 * std::sqrt(num_slots - num_to_add));
      if (batch_size < 10U) {
        batch_size = 0;
      }
      //      uint32_t batchKeyPrefix = 300000000;
      //      KeyGen batch_begin(batchKeyPrefix + 0);
      //      KeyGen batch_end(batchKeyPrefix + batch_size);
      std::string batch_str = prefix + "batch";
      KeyGen batch_begin(batch_str, 0);
      KeyGen batch_end(batch_str, batch_size);

      // Batch never (successfully) added, but used for querying FP rate
      //      uint32_t notKeyPrefix = 400000000;
      //      KeyGen other_keys_begin(notKeyPrefix + 0);
      //      KeyGen other_keys_end(notKeyPrefix + FLAGS_max_check);

      std::string not_str = prefix + "not";
      KeyGen other_keys_begin(not_str, 0);
      KeyGen other_keys_end(not_str, FLAGS_max_check);

      double overhead_ratio = 1.0 * num_slots / num_to_add;
      if (FLAGS_verbose) {
        fprintf(stderr, "Adding(%s) %u / %u   Overhead: %g   Batch size: %u\n",
                test_interleaved ? "i" : "s", (unsigned)num_to_add,
                (unsigned)num_slots, overhead_ratio, (unsigned)batch_size);
      }

      // Vary bytes for InterleavedSoln to use number of solution columns
      // from 0 to max allowed by ResultRow type (and used by SimpleSoln).
      // Specifically include 0 and max, and otherwise skew toward max.
      uint32_t max_ibytes =
          static_cast<uint32_t>(sizeof(ResultRow) * num_slots);
      size_t ibytes;
      // 迭代轮数的作用2：设置不同大小的InterleavedSoln
      if (i == 0) {
        ibytes = 0;
      } else if (i == 1) {
        ibytes = max_ibytes;
      } else {
        // Skewed 偏向于更大的随机数
        ibytes =
            std::max(rnd.Uniformish(max_ibytes), rnd.Uniformish(max_ibytes));
      }
      std::unique_ptr<char[]> idata(new char[ibytes]);
      InterleavedSoln isoln(idata.get(), ibytes);

      std::cout << "生成HBaseKey 定义soln时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms， 其中 overhead_ratio = " << overhead_ratio << std::endl;
      SimpleSoln soln;
      Hasher hasher;
      bool first_single;
      bool second_single;
      bool batch_success;
      {
        Banding banding;
        // Traditional solve for a fixed set.
        ASSERT_TRUE(
            banding.ResetAndFindSeedToSolve(num_slots, keys_begin, keys_end));

        //        banding.printCoeRow((int )num_slots);
        //        banding.printResultRow((int )num_slots);

        Index occupied_count = banding.GetOccupiedCount();
        Index more_added = 0;

        std::cout << "banding的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms， 此时添加了num_to_add个key：" << num_to_add << "banding占据数量" << occupied_count<< std::endl;
        if (TypeParam::kHomogeneous || overhead_ratio < 1.01 ||
            batch_size == 0) {
          // 可以看出来backtracking是用来记录banding失败的
          // Homogeneous not compatible with backtracking because add
          // 在开销小于1.01时，我们不应该要求更多
          // doesn't fail. Small overhead ratio too packed to expect more
          first_single = false;
          second_single = false;
          batch_success = false;
        } else {
          // Now to test backtracking, starting with guaranteed fail. By using
          // the keys that will be used to test FP rate, we are then doing an
          // extra check that after backtracking there are no remnants (e.g. in
          // result side of banding) of these entries.
          KeyGen other_keys_too_big_end = other_keys_begin;
          other_keys_too_big_end += num_to_add; // 一般这个num_to_add是大于nums_slots的一半的，所以这里再添加num_to_add一定会超过filter容量
          banding.EnsureBacktrackSize(std::max(num_to_add, batch_size));
          EXPECT_FALSE(banding.AddRangeOrRollBack(other_keys_begin,
                                                  other_keys_too_big_end));// 由于超过容量所以失败
          EXPECT_EQ(occupied_count, banding.GetOccupiedCount());

          // Check that we still have a good chance of adding a couple more
          // individually
          first_single = banding.Add(*one_more);
          second_single = banding.Add(*two_more);
          more_added += (first_single ? 1 : 0) + (second_single ? 1 : 0);
          total_singles += 2U;
          total_single_failures += 2U - more_added;

          // Or as a batch
          batch_success = banding.AddRangeOrRollBack(batch_begin, batch_end);
          ++total_batch;
          if (batch_success) {
            more_added += batch_size;
            ++total_batch_successes;
          }
          EXPECT_LE(banding.GetOccupiedCount(), occupied_count + more_added);
        }

        std::cout << "测试backtracking的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms，是否被跳过："<<(TypeParam::kHomogeneous || overhead_ratio < 1.01 ||
                                                                                                                       batch_size == 0 )<< std::endl;

        // Also verify that redundant adds are OK (no effect)
        //        ASSERT_TRUE(
        //            banding.AddRange(keys_begin, KeyGen(addedKeyPrefix + num_to_add / 8)));
        ASSERT_TRUE(
            banding.AddRange(keys_begin, KeyGen(added_str, num_to_add / 8)));
        EXPECT_LE(banding.GetOccupiedCount(), occupied_count + more_added);

        std::cout << "测试添加冗余key的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms" << std::endl;

        // Now back-substitution
        soln.BackSubstFrom(banding);
        if (test_interleaved) {
          isoln.BackSubstFrom(banding);
        }

        std::cout << "soln求解的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms，是否跳过isoln："<<!test_interleaved<< std::endl;

        Seed reseeds = banding.GetOrdinalSeed();
        total_reseeds += reseeds;

        EXPECT_LE(reseeds, 8 + log2_thoroughness); // 迷惑，其实干嘛要做log运算我是很迷惑的
        if (reseeds > log2_thoroughness + 1) {
          fprintf(
              stderr, "%s high reseeds at %u, %u/%u: %u\n",
              reseeds > log2_thoroughness + 8 ? "ERROR Extremely" : "Somewhat",
              static_cast<unsigned>(i), static_cast<unsigned>(num_to_add),
              static_cast<unsigned>(num_slots), static_cast<unsigned>(reseeds));
        }

        if (reseeds > 0) {
          // "Expand" test: given a failed construction, how likely is it to
          // pass with same seed and more slots. At each step, we increase
          // enough to ensure there is at least one shift within each coeff
          // block.
          ++total_expand_trials;
          Index expand_count = 0;
          Index ex_slots = num_slots;
          banding.SetOrdinalSeed(0);
          for (;; ++expand_count) {
            ASSERT_LE(expand_count, log2_thoroughness);
            ex_slots += ex_slots / kCoeffBits;
            if (test_interleaved) {
              ex_slots = InterleavedSoln::RoundUpNumSlots(ex_slots);
            }
            banding.Reset(ex_slots);
            bool success = banding.AddRange(keys_begin, keys_end);
            if (success) {
              break;
            }
          }
          total_expand_failures += expand_count;
          total_expand_overhead += 1.0 * (ex_slots - num_slots) / num_slots;
        }

        hasher.SetOrdinalSeed(reseeds);
        std::cout << "测试通过小幅度增加空间，是否可以增加banding的成功率的时间（expand）：" << timeStub.ElapsedNanos(true)/1000/1000
                  << "ms，是否被跳过："<<!(reseeds > 0)
                  << ", seed = "<< reseeds<< std::endl;
      }
      // soln and hasher now independent of Banding object

      // Verify keys added
      KeyGen cur = keys_begin;
      while (cur != keys_end) {
        ASSERT_TRUE(soln.FilterQuery(*cur, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*cur, hasher));
        ++cur;
      }
      // We (maybe) snuck these in!
      if (first_single) {
        ASSERT_TRUE(soln.FilterQuery(*one_more, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*one_more, hasher));
      }
      if (second_single) {
        ASSERT_TRUE(soln.FilterQuery(*two_more, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*two_more, hasher));
      }
      if (batch_success) {
        cur = batch_begin;
        while (cur != batch_end) {
          ASSERT_TRUE(soln.FilterQuery(*cur, hasher));
          ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*cur, hasher));
          ++cur;
        }
      }

      std::cout << "查询已添加的key时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms，是否跳过isoln："<<!(test_interleaved) << std::endl;

      // Check FP rate (depends only on number of result bits == solution
      // columns)
      Index fp_count = 0;
      cur = other_keys_begin;
      {
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          bool fp = soln.FilterQuery(*cur, hasher);
          fp_count += fp ? 1 : 0;
          ++cur;
        }
        soln_query_nanos += timer.ElapsedNanos();
        soln_query_count += FLAGS_max_check;
      }
      {
        double expected_fp_count = soln.ExpectedFpRate() * FLAGS_max_check;
        // For expected FP rate, also include false positives due to collisions
        // in Hash value. (Negligible for 64-bit, can matter for 32-bit.)
        double correction =
            FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

        // NOTE: rare violations expected with kHomogeneous
        EXPECT_LE(fp_count,
                  FrequentPoissonUpperBound(expected_fp_count + correction)); // 泊松分布对比误差
        EXPECT_GE(fp_count,
                  FrequentPoissonLowerBound(expected_fp_count + correction));
      }
      total_fp_count += fp_count;

      std::cout << "查询soln 未添加的key，并计算fp 时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms"<< std::endl;

      // And also check FP rate for isoln
      if (test_interleaved) {
        Index ifp_count = 0;
        cur = other_keys_begin;
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          ifp_count += isoln.FilterQuery(*cur, hasher) ? 1 : 0;
          ++cur;
        }
        isoln_query_nanos += timer.ElapsedNanos();
        isoln_query_count += FLAGS_max_check;
        {
          double expected_fp_count = isoln.ExpectedFpRate() * FLAGS_max_check;
          // For expected FP rate, also include false positives due to
          // collisions in Hash value. (Negligible for 64-bit, can matter for
          // 32-bit.)
          double correction =
              FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

          // NOTE: rare violations expected with kHomogeneous
          EXPECT_LE(ifp_count,
                    FrequentPoissonUpperBound(expected_fp_count + correction));

          // FIXME: why sometimes can we slightly "beat the odds"?
          // 这个测试有时会不通过，因此引入了一个因子来降低期望的误差数，但我感觉这是个好事情？
          // (0.95 factor should not be needed)
          EXPECT_GE(ifp_count, FrequentPoissonLowerBound(
                                   0.95 * expected_fp_count + correction));
        }
        // Since the bits used in isoln are a subset of the bits used in soln,
        // it cannot have fewer FPs
        // 从这里看出来 isoln对比soln的优势不在fp，在fp上应该是持平的，所以可能在内存上有优化
        EXPECT_GE(ifp_count, fp_count);
      }

      std::cout << "查询isoln 未添加的key，并计算fp 时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms， 是否跳过"<< !test_interleaved<< std::endl;

      // And compare to Bloom time, for fun
      if (ibytes >= /* minimum Bloom impl bytes*/ 64) {
        Index bfp_count = 0;
        cur = other_keys_begin;
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          uint64_t h = hasher.GetHash(*cur);
          uint32_t h1 = ROCKSDB_NAMESPACE::Lower32of64(h);
          uint32_t h2 = sizeof(Hash) >= 8 ? ROCKSDB_NAMESPACE::Upper32of64(h)
                                          : h1 * 0x9e3779b9;
          bfp_count +=
              ROCKSDB_NAMESPACE::FastLocalBloomImpl::HashMayMatch(
                  h1, h2, static_cast<uint32_t>(ibytes), 6, idata.get())
                  ? 1
                  : 0;
          ++cur;
        }
        bloom_query_nanos += timer.ElapsedNanos();
        // ensure bfp_count is used
        ASSERT_LT(bfp_count, FLAGS_max_check);
      }

      std::cout << "对比bloom的时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms， 此时第" << i <<  "次迭代结束"<< std::endl;
    }

    // "outside" == key not in original set so either negative or false positive
    fprintf(stderr,
            "Simple      outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * soln_query_nanos / soln_query_count);
    fprintf(stderr,
            "Interleaved outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * isoln_query_nanos / isoln_query_count);
    fprintf(stderr,
            "Bloom       outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * bloom_query_nanos / soln_query_count);

    if (TypeParam::kHomogeneous) {
      EXPECT_EQ(total_reseeds, 0U); // 说明kHomogeneous在banding中不会失败，因此不会改变种子的值
    } else {
      double average_reseeds = 1.0 * total_reseeds / FLAGS_thoroughness;
      fprintf(stderr, "Average re-seeds: %g\n", average_reseeds);
      // Values above were chosen to target around 50% chance of encoding
      // success rate (average of 1.0 re-seeds) or slightly better. But 1.15 is
      // also close enough.
      EXPECT_LE(total_reseeds,
                InfrequentPoissonUpperBound(1.15 * expected_reseeds *
                                            FLAGS_thoroughness));
      // Would use 0.85 here instead of 0.75, but
      // TypesAndSettings_Hash32_SmallKeyGen can "beat the odds" because of
      // sequential keys with a small, cheap hash function. We accept that
      // there are surely inputs that are somewhat bad for this setup, but
      // these somewhat good inputs are probably more likely.
      EXPECT_GE(total_reseeds,
                InfrequentPoissonLowerBound(0.75 * expected_reseeds *
                                            FLAGS_thoroughness));
    }

    std::cout << "测试 reseed 是否符合预期的时间（测试逻辑不太懂）：" << timeStub.ElapsedNanos(true)/1000/1000
              << "ms"<< std::endl;

    if (total_expand_trials > 0) {
      double average_expand_failures =
          1.0 * total_expand_failures / total_expand_trials;
      fprintf(stderr, "Average expand failures, and overhead: %g, %g\n",
              average_expand_failures,
              total_expand_overhead / total_expand_trials);
      // Seems to be a generous allowance
      EXPECT_LE(total_expand_failures,
                InfrequentPoissonUpperBound(1.0 * total_expand_trials));
    } else {
      fprintf(stderr, "Average expand failures: N/A\n");
    }

    if (total_singles > 0) {
      double single_failure_rate = 1.0 * total_single_failures / total_singles;
      fprintf(stderr, "Add'l single, failure rate: %g\n", single_failure_rate);
      // A rough bound (one sided) based on nothing in particular
      double expected_single_failures = 1.0 * total_singles /
                                        (sizeof(CoeffRow) == 16 ? 128
                                         : TypeParam::kUseSmash ? 64
                                                                : 32);
      EXPECT_LE(total_single_failures,
                InfrequentPoissonUpperBound(expected_single_failures));
    }

    if (total_batch > 0) {
      // Counting successes here for Poisson to approximate the Binomial
      // distribution.
      // A rough bound (one sided) based on nothing in particular.
      double expected_batch_successes = 1.0 * total_batch / 2;
      uint64_t lower_bound =
          InfrequentPoissonLowerBound(expected_batch_successes);
      fprintf(stderr, "Add'l batch, success rate: %g (>= %g)\n",
              1.0 * total_batch_successes / total_batch,
              1.0 * lower_bound / total_batch);
      EXPECT_GE(total_batch_successes, lower_bound);
    }

    {
      uint64_t total_checked = uint64_t{FLAGS_max_check} * FLAGS_thoroughness;
      double expected_total_fp_count =
          total_checked * std::pow(0.5, 8U * sizeof(ResultRow));
      // For expected FP rate, also include false positives due to collisions
      // in Hash value. (Negligible for 64-bit, can matter for 32-bit.)
      double average_added = 1.0 * total_added / FLAGS_thoroughness;
      expected_total_fp_count +=
          total_checked * ExpectedCollisionFpRate(Hasher(), average_added);

      uint64_t upper_bound =
          InfrequentPoissonUpperBound(expected_total_fp_count);
      uint64_t lower_bound =
          InfrequentPoissonLowerBound(expected_total_fp_count);
      fprintf(stderr, "Average FP rate: %g (~= %g, <= %g, >= %g)\n",
              1.0 * total_fp_count / total_checked,
              expected_total_fp_count / total_checked,
              1.0 * upper_bound / total_checked,
              1.0 * lower_bound / total_checked);
      EXPECT_LE(total_fp_count, upper_bound);
      EXPECT_GE(total_fp_count, lower_bound);
    }

    std::cout << "计算expand、single、batch、avg fp rate是否符合预期的时间（应该没啥计算量） ：" << timeStub.ElapsedNanos(true)/1000/1000
              << "ms， 此时enum["<< cs << "]结束" << std::endl;
  }
}

int main(int argc, char** argv) {
  ROCKSDB_NAMESPACE::port::InstallStackTraceHandler();
  ::testing::InitGoogleTest(&argc, argv);
#ifdef GFLAGS
  ParseCommandLineFlags(&argc, &argv, true);
#endif  // GFLAGS
  return RUN_ALL_TESTS();
}
