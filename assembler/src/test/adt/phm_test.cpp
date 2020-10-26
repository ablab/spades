#include "utils/logger/logger.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils/stl_utils.hpp"
#include "boomphf/BooPHF.h"

#include <iostream>
#include <vector>
#include <algorithm>

#define XXH_INLINE_ALL
#include "xxh/xxhash.h"

#include <gtest/gtest.h>

struct PHMTest : public ::testing::Test {
    class Hasher {
      public:
        std::pair<uint64_t, uint64_t> operator()(uint64_t key) const {
            auto res = XXH3_128bits(&key, sizeof(key));
            return { res.high64, res.low64 };
        }
    };
    
    using phm = boomphf::mphf<Hasher>;
    
    PHMTest() {
        INFO("setup fixture");

        for (size_t i = 0; i < 10005000; ++i)
            vals_.push_back(i);

        phm_.init(vals_.size(),
                  phm::ConflictPolicy::Warning, 1, 0.03, 10);
        phm_.build(boomphf::range(vals_.begin(), vals_.end()));

        INFO("PHM built. Size: " << phm_.size()
             << ", mem size:" << phm_.mem_size()
             << ", collision prob: " << phm_.prob_collision()
             << ", final hash size: " << phm_.last_level_size());
    }
    ~PHMTest() {
        INFO("teardown fixture");
    }

    phm phm_;
    std::vector<uint64_t> vals_;
};

TEST_F(PHMTest, empty_test) {
    ASSERT_GT(phm_.size(), 0);
}

TEST_F(PHMTest, id_test) {
    for (uint64_t entry : vals_) {
        uint64_t idx = phm_.lookup(entry);
        uint64_t invalid = phm::NOT_FOUND;
        EXPECT_NE(idx, invalid);
        EXPECT_LT(idx, phm_.size());
    }
}

TEST_F(PHMTest, collision_test) {
    std::vector<bool> marks(vals_.size());
    for (uint64_t entry : vals_) {
        uint64_t idx = phm_.lookup(entry);
        uint64_t invalid = phm::NOT_FOUND;
        EXPECT_NE(idx, invalid);
        EXPECT_LT(idx, phm_.size());
        EXPECT_NE(marks[idx], true);
        marks[idx] = true;
    }
}

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

GTEST_API_ int main(int argc, char **argv) {
  create_console_logger();

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
