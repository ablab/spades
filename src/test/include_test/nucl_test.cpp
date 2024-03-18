//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence/nucl.hpp"
#include <gtest/gtest.h>

TEST( Nucl, Test ) {
    EXPECT_EQ('A', nucl(0));
    EXPECT_EQ('C', nucl(1));
    EXPECT_EQ('G', nucl(2));
    EXPECT_EQ('T', nucl(3));
    EXPECT_EQ(0, dignucl('A'));
    EXPECT_EQ(1, dignucl('C'));
    EXPECT_EQ(2, dignucl('G'));
    EXPECT_EQ(3, dignucl('T'));
    EXPECT_EQ(3, complement(0));
    EXPECT_EQ(2, complement(1));
    EXPECT_EQ(1, complement(2));
    EXPECT_EQ(0, complement(3));
    EXPECT_TRUE(is_nucl('A'));
    EXPECT_TRUE(is_nucl('C'));
    EXPECT_TRUE(is_nucl('G'));
    EXPECT_TRUE(is_nucl('T'));
    EXPECT_TRUE(is_nucl(0));
    EXPECT_TRUE(is_nucl(1));
    EXPECT_TRUE(is_nucl(2));
    EXPECT_TRUE(is_nucl(3));
    EXPECT_TRUE(!is_nucl('0'));
    EXPECT_TRUE(!is_nucl('1'));
}
