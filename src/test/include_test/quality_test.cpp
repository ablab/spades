//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence/quality.hpp"
#include <gtest/gtest.h>

TEST( Quality, Test ) {
    Quality q("0123456789");
    EXPECT_EQ('0', q[0]);
    EXPECT_EQ('6', q[6]);
    EXPECT_EQ('9', q[9]);
}
