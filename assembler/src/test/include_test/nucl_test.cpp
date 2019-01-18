//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence/nucl.hpp"
#include <gtest/gtest.h>

TEST( Nucl, Test ) {
    ASSERT_EQ('A', nucl(0));
    ASSERT_EQ('C', nucl(1));
    ASSERT_EQ('G', nucl(2));
    ASSERT_EQ('T', nucl(3));
    ASSERT_EQ(0, dignucl('A'));
    ASSERT_EQ(1, dignucl('C'));
    ASSERT_EQ(2, dignucl('G'));
    ASSERT_EQ(3, dignucl('T'));
    ASSERT_EQ(3, complement(0));
    ASSERT_EQ(2, complement(1));
    ASSERT_EQ(1, complement(2));
    ASSERT_EQ(0, complement(3));
    ASSERT_TRUE(is_nucl('A'));
    ASSERT_TRUE(is_nucl('C'));
    ASSERT_TRUE(is_nucl('G'));
    ASSERT_TRUE(is_nucl('T'));
    ASSERT_TRUE(is_nucl(0));
    ASSERT_TRUE(is_nucl(1));
    ASSERT_TRUE(is_nucl(2));
    ASSERT_TRUE(is_nucl(3));
    ASSERT_TRUE(!is_nucl('0'));
    ASSERT_TRUE(!is_nucl('1'));
}
