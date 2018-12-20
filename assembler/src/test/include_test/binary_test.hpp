//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/binary.hpp"
#include "utils/stl_utils.hpp"
#include <boost/test/unit_test.hpp>
#include <sstream>

template<typename T>
void SimpleCheck(const T &value = T()) {
    using namespace io::binary;
    std::stringstream str;
    BinWrite(str, value);
    T result;
    BinRead(str, result);
    BOOST_CHECK_EQUAL(value, result);
}

BOOST_AUTO_TEST_CASE(TestBinarySimple) {
    //Integrals
    SimpleCheck((char)rand());
    SimpleCheck((unsigned)rand());
    SimpleCheck((int)rand());
    SimpleCheck((unsigned long long)rand());

    //Floats
    SimpleCheck((float)rand() / 10.0f);
    SimpleCheck((double)rand() / 10.0);

    //Strings
    SimpleCheck<std::string>("A quick brown fox jumps over the lazy dog");
}

struct Boo {int i; float f;};

bool operator==(const Boo &l, const Boo &r) {
    return l.i == r.i && l.f == r.f;
}

std::ostream &operator<<(std::ostream& str, const Boo &b) {
    str << "Boo{" << b.i << "," << b.f << "}";
    return str;
}

BOOST_AUTO_TEST_CASE(TestBinaryContainers) {
    //POD
    SimpleCheck<Boo>({rand(), (float)rand() / 10.0f});

    //Vectors
    std::vector<unsigned> v = {4, 8, 15, 16, 23, 42};
    SimpleCheck(v);

    //Sets
    std::map<std::string, char> s = {{"Finn", 1}, {"Jake", 1000}, {"Bmo", 1000000}};
    SimpleCheck(s);
}

BOOST_AUTO_TEST_CASE(TestBinaryVariadic) {
    using namespace io::binary;

    char c1 = (char)rand(), c2;
    unsigned u1 = (unsigned)rand(), u2;
    float f1 = (float)rand() / 10.0f, f2;
    unsigned long long l1 = (unsigned long long)rand(), l2;

    std::stringstream str;
    BinWrite(str, c1, u1, f1, l1);
    BinRead(str, c2, u2, f2, l2);

    BOOST_CHECK_EQUAL(c1, c2);
    BOOST_CHECK_EQUAL(u1, u2);
    BOOST_CHECK_EQUAL(f1, f2);
    BOOST_CHECK_EQUAL(l1, l2);
}