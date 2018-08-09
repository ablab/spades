//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/test/unit_test.hpp>

#include "paired_info/histogram.hpp"

namespace omnigraph {

namespace de {

class alignas(2) Counter {
public:
     Counter() { ++count_; }
    ~Counter() { --count_; }

    static int Count() { return count_; }

private:
    static int count_;
};

int Counter::count_ = 0;

BOOST_AUTO_TEST_SUITE(histogram_tests)

BOOST_AUTO_TEST_CASE(StrongWeakPtrBasic) {
    using Ptr = StrongWeakPtr<Counter>;
    //Empty
    Ptr p1;
    BOOST_CHECK(!p1.owning());
    BOOST_CHECK_EQUAL(p1.get(), static_cast<Counter *>(nullptr));
    {
        //Owning
        Ptr p2(new Counter());
        auto ptr = p2.get();
        BOOST_CHECK(p2.owning());
        BOOST_CHECK_EQUAL(Counter::Count(), 1);
        //Moving
        Ptr p3 = std::move(p2);
        BOOST_CHECK_EQUAL(p3.get(), ptr);
        BOOST_CHECK_EQUAL(Counter::Count(), 1);
        //Non-owning
        Ptr p4(p3.get(), false);
        BOOST_CHECK_EQUAL(p4.get(), ptr);
        BOOST_CHECK_EQUAL(Counter::Count(), 1);
    }
    //Free
    BOOST_CHECK_EQUAL(Counter::Count(), 0);
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace de

} // namespace omnigraph
