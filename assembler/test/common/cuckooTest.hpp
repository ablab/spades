/* test for cuckoo.hpp */
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "cuckoo.hpp"

struct Hasher {
  size_t operator()(const int& value, size_t& hash_num, size_t seed) const {
    unsigned long l = 4 * hash_num + 1;
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
  ASSERT_EQUAL(map.size(), 3);
  ASSERT_EQUAL(map[1], 4);
  map[1] = 5;
  ASSERT_EQUAL(map[1], 5);

  hm map2(map);
  ASSERT_EQUAL(map2.size(), 3);
  ASSERT(map2.find(1) != map2.end());
  ASSERT(map2.find(4) != map2.end());
  ASSERT(map2.find(6) != map2.end());
  ASSERT(map2.find(8) == map2.end());

  hm map3;
  map3 = map;
  ASSERT_EQUAL(map3.size(), 3); 
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
  hm map(5, 0, 10000, 100, 1.2); 
  srand(42);
  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    map.insert(std::make_pair(t, 42));
  }
  size_t size = map.size();
  ASSERT_EQUAL(size, 99994);
  
  size_t n = 0;
  for (hm::iterator it = map.begin(); it != map.end(); ++it) {
    ++n;
  }
  ASSERT_EQUAL(size, n);

  n = 0;
  for (hm::const_iterator it = map.begin(); it != map.end(); ++it) {
    ++n;
  }
  ASSERT_EQUAL(size, n); 

  hm::iterator it;
  for (int i = 0; i < 10000; ++i) {
    int t = rand();
    it = map.find(t);
  }
  ASSERT_EQUAL(size, map.size());
  
  for (int i = 0; i < 10000; ++i) {
    int t = rand();
    n -= map.erase(t);
  }
  ASSERT_EQUAL(map.size(), n);

  for (int i = 0; i < 10000; ++i) {
    int t = rand();
    it = map.find(t);
    map.erase(it);
    --n;
  }
  ASSERT_EQUAL(map.size(), n);
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
  ASSERT_EQUAL(map.size(), 1);
  ASSERT_EQUAL(map.empty(), false);
  ASSERT_EQUAL(map.count(21), 1);
  ASSERT_EQUAL(map.count(22), 0);

  ASSERT_EQUAL(map[21], 42);
  ASSERT(map.equal_range(21).first == map.begin());
  ASSERT(map.equal_range(21).second == map.end());

  hm map_new;
  map.swap(map_new);
  ASSERT_EQUAL(map.size(), 0);
  ASSERT_EQUAL(map_new.size(), 1);

  map_new.clear();
  ASSERT_EQUAL(map_new.size(), 0);
  ASSERT_EQUAL(map_new.empty(), true);
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
  ASSERT_EQUAL(map1.size(), size);

  hm map2;
  map2.insert(map1.begin(), map1.end());
  ASSERT_EQUAL(map2[11], 12);
  ASSERT_EQUAL(map2.size(), size);

  hm map3(map1.begin(), map1.end());
  ASSERT_EQUAL(map3[11], 12);
  ASSERT_EQUAL(map3.size(), size);

  map3.erase(map3.begin(), map3.end());
  ASSERT_EQUAL(map3.size(), 0);

  map3.swap(map2);
  ASSERT_EQUAL(map3.size(), size);
  ASSERT_EQUAL(map2.size(), 0);
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
  ASSERT_EQUAL((*it).second, 12);
  ASSERT(++(map.equal_range(11).first) == map.equal_range(11).second);
  ASSERT_EQUAL(map.count(11), 1);
  ASSERT_EQUAL(map.size(), 30);
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
