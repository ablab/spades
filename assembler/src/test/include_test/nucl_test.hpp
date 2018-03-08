//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <boost/test/unit_test.hpp>
#include "sequence/nucl.hpp"

BOOST_AUTO_TEST_CASE( TestNucl ) {
    BOOST_CHECK_EQUAL('A', nucl(0));
    BOOST_CHECK_EQUAL('C', nucl(1));
    BOOST_CHECK_EQUAL('G', nucl(2));
    BOOST_CHECK_EQUAL('T', nucl(3));
    BOOST_CHECK_EQUAL(0, dignucl('A'));
    BOOST_CHECK_EQUAL(1, dignucl('C'));
    BOOST_CHECK_EQUAL(2, dignucl('G'));
    BOOST_CHECK_EQUAL(3, dignucl('T'));
    BOOST_CHECK_EQUAL(3, complement(0));
    BOOST_CHECK_EQUAL(2, complement(1));
    BOOST_CHECK_EQUAL(1, complement(2));
    BOOST_CHECK_EQUAL(0, complement(3));
    BOOST_CHECK(is_nucl('A'));
    BOOST_CHECK(is_nucl('C'));
    BOOST_CHECK(is_nucl('G'));
    BOOST_CHECK(is_nucl('T'));
    BOOST_CHECK(is_nucl(0));
    BOOST_CHECK(is_nucl(1));
    BOOST_CHECK(is_nucl(2));
    BOOST_CHECK(is_nucl(3));
    BOOST_CHECK(!is_nucl('0'));
    BOOST_CHECK(!is_nucl('1'));
}
