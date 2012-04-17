//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef TEST_CUCKOOTEST_HPP_
#define TEST_CUCKOOTEST_HPP_

#include "cute/cute.h"
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <vector>
#include <functional>
#include "cuckoo.hpp"

struct Hasher {
  size_t operator()(const int& value, size_t& hash_num, size_t seed) const {
    uint64_t l = 4 * hash_num + 1;
    return (size_t)(l * value % 1000000007);
  }
};

typedef cuckoo<int, int, Hasher, std::equal_to<int> > hm;

// *** Test for: ***
// cuckoo()
// cuckoo(...)
// operator=(...)
// operator[](...)
// find(...)
// size()
// begin()
// end()
// operator==
// operator!=
void TestCuckooCreation() {
  hm map;
  map[1] = 4;
  map[4] = 1;
  map[6];
  ASSERT_EQUAL(3, map.size());
  ASSERT_EQUAL(4, map[1]);
  map[1] = 5;
  ASSERT_EQUAL(5, map[1]);

  hm map2(map);
  ASSERT_EQUAL(3, map2.size());
  ASSERT(map2.find(1) != map2.end());
  ASSERT(map2.find(4) != map2.end());
  ASSERT(map2.find(6) != map2.end());
  ASSERT(map2.find(8) == map2.end());

  hm map3;
  map3 = map;
  ASSERT_EQUAL(3, map3.size());
  ASSERT(map3.find(1) != map3.end());
  ASSERT(map3.find(4) != map3.end());
  ASSERT(map3.find(6) != map3.end());
  ASSERT(map3.find(8) == map3.end());

  ASSERT(map2.begin() != map3.begin());
  ASSERT(map == map);
  ASSERT(map != map2);
}

// *** Test for: ***
// cuckoo(...)
// insert(const key&)
// size()
// find(...)
// erase(const key&)
// erase(iterator)
// begin()
// end()
// iterator()
// iterator::operator++
// iterator::const_iterator(...)
void TestCuckooOperations() {
  hm map(5, 0, 1000, 100, 1.2);
  srand(42);
  for (int i = 0; i < 5000; ++i) {
    int t = rand();
    map.insert(std::make_pair(t, 42));
  }
  size_t size = map.size();
  ASSERT_EQUAL(5000, size);

  size_t n = 0;
  for (hm::iterator it = map.begin(); it != map.end(); ++it) {
    ++n;
  }
  ASSERT_EQUAL(n, size);

  n = 0;
  for (hm::const_iterator it = map.begin(); it != map.end(); ++it) {
    ++n;
  }
  ASSERT_EQUAL(n, size);

  hm::iterator it;
  for (int i = 0; i < 1000; ++i) {
    int t = rand();
    it = map.find(t);
  }
  ASSERT_EQUAL(size, map.size());

  for (int i = 0; i < 1000; ++i) {
    int t = rand();
    n -= map.erase(t);
  }
  ASSERT_EQUAL(n, map.size());

  for (int i = 0; i < 1000; ++i) {
    int t = rand();
    it = map.find(t);
    if (it != map.end()) {
      map.erase(it);
      --n;
    }
  }
  ASSERT_EQUAL(n, map.size());
}

// *** Test for: ***
// cuckoo()
// insert(const value&)
// size()
// empty()
// count()
// equal_range(...)
// begin()
// clear()
// swap(...)
void TestCuckooSize() {
  hm map;
  map.insert(std::make_pair(21, 42));
  ASSERT_EQUAL(1, map.size());
  ASSERT_EQUAL(false, map.empty());
  ASSERT_EQUAL(1, map.count(21));
  ASSERT_EQUAL(0, map.count(22));

  ASSERT_EQUAL(42, map[21]);
  ASSERT(map.equal_range(21).first == map.begin());
  ASSERT(map.equal_range(21).second == map.end());

  hm map_new;
  map.swap(map_new);
  ASSERT_EQUAL(0, map.size());
  ASSERT_EQUAL(1, map_new.size());

  map_new.clear();
  ASSERT_EQUAL(0, map_new.size());
  ASSERT_EQUAL(true, map_new.empty());
}

// *** Test for: ***
// operator[]
// size()
// insert(<range>)
// cuckoo(<range>)
// erase(<range>)
// swap(...)
// size()
// begin()
// end()
void TestCuckooRanges() {
  size_t size = 30;
  hm map1;
  for (size_t i = 1; i <= size; ++i) {
    map1[i * 2 - 1] = i * 2;
  }
  ASSERT_EQUAL(size, map1.size());

  hm map2;
  map2.insert(map1.begin(), map1.end());
  ASSERT_EQUAL(12, map2[11]);
  ASSERT_EQUAL(size, map2.size());

  hm map3(map1.begin(), map1.end());
  ASSERT_EQUAL(12, map3[11]);
  ASSERT_EQUAL(size, map3.size());

  map3.erase(map3.begin(), map3.end());
  ASSERT_EQUAL(0, map3.size());

  map3.swap(map2);
  ASSERT_EQUAL(size, map3.size());
  ASSERT_EQUAL(0, map2.size());
}

// *** Test for: ***
// begin() const
// end() const
// find(...) const
// equal_range(...) const
// count()
// size()
void const_tester(const hm& map) {
  hm::const_iterator it;
  for (it = map.begin(); it != map.end(); ++it) {
    ASSERT_EQUAL((*it).second, (*it).first + 1);
  }
  it = map.find(11);
  ASSERT_EQUAL(12, (*it).second);
  ASSERT(++(map.equal_range(11).first) == map.equal_range(11).second);
  ASSERT_EQUAL(1, map.count(11));
  ASSERT_EQUAL(30, map.size());
}

void TestCuckooConst() {
  hm map;
  for (size_t i = 1; i <= 30; ++i) {
    map[i * 2 - 1] = i * 2;
  }
  const_tester(map);
}

cute::suite CuckooSuite() {
  cute::suite s;
  s.push_back(CUTE(TestCuckooCreation));
  s.push_back(CUTE(TestCuckooOperations));
  s.push_back(CUTE(TestCuckooSize));
  s.push_back(CUTE(TestCuckooRanges));
  s.push_back(CUTE(TestCuckooConst));
  return s;
}

#endif /* TEST_CUCKOOTEST_HPP_ */
