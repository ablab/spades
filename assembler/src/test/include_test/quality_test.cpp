//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence/quality.hpp"
#include <gtest/gtest.h>

TEST( Quality, Test ) {
    Quality q("0123456789");
    ASSERT_EQ('0', q[0]);
    ASSERT_EQ('6', q[6]);
    ASSERT_EQ('9', q[9]);
}
