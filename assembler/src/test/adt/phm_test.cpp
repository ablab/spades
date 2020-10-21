#define BOOST_TEST_MODULE phm_test
#include <boost/test/included/unit_test.hpp>

#include "boomphf/BooPHF.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "utils/stl_utils.hpp"

#define XXH_INLINE_ALL
#include "xxh/xxhash.h"

struct F {
    class Hasher {
      public:
        std::pair<uint64_t, uint64_t> operator()(uint64_t key) const {
            auto res = XXH3_128bits(&key, sizeof(key));
            return { res.high64, res.low64 };
        }
    };
    
    using phm = boomphf::mphf<Hasher>;
    
    F() {
        BOOST_TEST_MESSAGE("setup fixture");

        for (size_t i = 0; i < 100500; ++i)
            vals_.push_back(i);

        phm_.init(vals_.size());
        phm_.build(boomphf::range(vals_.begin(), vals_.end()));

        BOOST_TEST_MESSAGE("PHM built. Size: " << phm_.size() << ", mem size:" << phm_.mem_size() << ", collision prob: " << phm_.prob_collision());
    }
    ~F() {
        BOOST_TEST_MESSAGE( "teardown fixture" );
    }

    phm phm_;
    std::vector<uint64_t> vals_;
};

BOOST_FIXTURE_TEST_CASE(empty_test, F) {
  BOOST_CHECK(phm_.size() > 0);
}

BOOST_FIXTURE_TEST_CASE(id_test, F) {
    for (uint64_t entry : vals_) {
        uint64_t idx = phm_.lookup(entry);
        uint64_t invalid = phm::NOT_FOUND;
        BOOST_CHECK_NE(idx, invalid);
        BOOST_CHECK_LT(idx, phm_.size());
    }
}

BOOST_FIXTURE_TEST_CASE(collision_test, F) {
    std::vector<bool> marks(vals_.size());
    for (uint64_t entry : vals_) {
        uint64_t idx = phm_.lookup(entry);
        uint64_t invalid = phm::NOT_FOUND;
        BOOST_CHECK_NE(idx, invalid);
        BOOST_CHECK_LT(idx, phm_.size());
        BOOST_CHECK_NE(marks[idx], true);
        marks[idx] = true;
    }
}
