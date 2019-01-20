//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/binary/binary.hpp"
#include "utils/stl_utils.hpp"
#include <sstream>
#include <gtest/gtest.h>


template<typename T>
void SimpleCheck(const T &value = T()) {
    using namespace io::binary;
    std::stringstream str;
    BinWrite(str, value);
    T result;
    BinRead(str, result);
    ASSERT_EQ(value, result);
}

TEST(Binary, Simple) {
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

struct Boo {int i; float f;
template <typename Archive>
void BinArchive(Archive &ar) {
    ar(i, f);
}
};

bool operator==(const Boo &l, const Boo &r) {
    return l.i == r.i && l.f == r.f;
}

std::ostream &operator<<(std::ostream& str, const Boo &b) {
    str << "Boo{" << b.i << "," << b.f << "}";
    return str;
}

TEST(Binary, Containers) {
    //POD
    SimpleCheck<Boo>({rand(), (float)rand() / 10.0f});

    //Vectors
    std::vector<unsigned> v = {4, 8, 15, 16, 23, 42};
    SimpleCheck(v);

    //Sets
    std::map<std::string, char> s = {{"Finn", 1}, {"Jake", 1000}, {"Bmo", 1000000}};
    SimpleCheck(s);
}

TEST(Binary, Variadic) {
    using namespace io::binary;

    char c1 = (char)rand(), c2;
    unsigned u1 = (unsigned)rand(), u2;
    float f1 = (float)rand() / 10.0f, f2;
    unsigned long long l1 = (unsigned long long)rand(), l2;

    std::stringstream str;
    BinWrite(str, c1, u1, f1, l1);
    BinRead(str, c2, u2, f2, l2);

    ASSERT_EQ(c1, c2);
    ASSERT_EQ(u1, u2);
    ASSERT_EQ(f1, f2);
    ASSERT_EQ(l1, l2);
}
