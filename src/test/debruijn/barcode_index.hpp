#pragma once

#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "common/barcode_index/barcode_index.hpp"

namespace debruijn_graph {
    BOOST_AUTO_TEST_SUITE(barcode_index_tests)

    BOOST_AUTO_TEST_CASE(TestingTest) {
        BOOST_CHECK_EQUAL(1, 1);
    }


    BOOST_AUTO_TEST_SUITE_END()
}