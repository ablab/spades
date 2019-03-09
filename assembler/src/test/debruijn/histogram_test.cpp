//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "paired_info/histogram.hpp"
#include "paired_info/histptr.hpp"

#include <gtest/gtest.h>

using namespace omnigraph::de;

class alignas(2) Counter {
public:
     Counter() { ++count_; }
    ~Counter() { --count_; }

    static int Count() { return count_; }

private:
    static int count_;
};

int Counter::count_ = 0;

TEST(Histogram, StrongWeakPtrBasic) {
    using Ptr = StrongWeakPtr<Counter>;
    //Empty
    Ptr p1;
    EXPECT_FALSE(p1.owning());
    EXPECT_EQ(p1.get(), static_cast<Counter *>(nullptr));
    {
        //Owning
        Ptr p2(new Counter());
        auto ptr = p2.get();
        EXPECT_TRUE(p2.owning());
        EXPECT_EQ(1, Counter::Count());
        //Moving
        Ptr p3 = std::move(p2);
        EXPECT_EQ(p3.get(), ptr);
        EXPECT_EQ(1, Counter::Count());
        //Non-owning
        Ptr p4(p3.get(), false);
        EXPECT_EQ(p4.get(), ptr);
        EXPECT_EQ(1, Counter::Count());
    }
    //Free
    EXPECT_EQ(Counter::Count(), 0);
}
