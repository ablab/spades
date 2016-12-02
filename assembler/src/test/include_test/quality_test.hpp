//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <boost/test/unit_test.hpp>
#include "sequence/quality.hpp"

BOOST_AUTO_TEST_CASE ( QualityTest ) {
    Quality q("0123456789");
    BOOST_CHECK_EQUAL('0', q[0]);
    BOOST_CHECK_EQUAL('6', q[6]);
    BOOST_CHECK_EQUAL('9', q[9]);
}
