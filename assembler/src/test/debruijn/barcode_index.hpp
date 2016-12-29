#pragma once

#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "projects/tslr_resolver/barcode_mapper.hpp"

namespace debruijn_graph {
    BOOST_AUTO_TEST_SUITE(barcode_index_tests)

    BOOST_AUTO_TEST_CASE(TestingTest) {
        BOOST_CHECK_EQUAL(1, 55);
    }

    BOOST_AUTO_TEST_SUITE_END()
}