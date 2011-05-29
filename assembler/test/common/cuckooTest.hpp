/* test for cuckoo.hpp */
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "cuckoo.hpp"

struct Hasher {
  size_t operator()(const int& value, size_t& hash_num) const {
    unsigned long l = 4 * hash_num + 1;
    return (size_t)(l * value % 1000000007);
  }
};

typedef cuckoo<int, int, Hasher, std::equal_to<int> > hm; 

void TestCuckoo() {
  srand(42);
  hm map(5, 10000, 100, 1.2); 
  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    map.insert(std::make_pair(t, 42));
  }
  size_t size = map.size();
  ASSERT_EQUAL(size, 99994);
  
  hm::iterator it;
  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    it = map.find(t);
  }
  ASSERT_EQUAL(size, map.size());
  
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

  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    n -= map.erase(t);
  }
  ASSERT_EQUAL(map.size(), n);

  map.clear();
  ASSERT_EQUAL(map.size(), 0);
  ASSERT_EQUAL(map.empty(), true);

  int t = rand();
  map.insert(std::make_pair(t, 42));
  ASSERT_EQUAL(map.size(), 1);
  ASSERT_EQUAL(map.empty(), false);
  ASSERT_EQUAL(map.count(t), 1);
  ASSERT_EQUAL(map.count(t + 1), 0);
}

void TestCuckooCopiing() {
  hm map;
  map[1] = 4;
  map[4] = 1;
  map[6];
  ASSERT_EQUAL(map.size(), 3);

  hm map2;
  map2 = map;
  ASSERT_EQUAL(map2.size(), 3);
  ASSERT(map2.find(1) != map2.end());
  ASSERT(map2.find(4) != map2.end());
  ASSERT(map2.find(6) != map2.end());
  ASSERT(map2.find(8) == map2.end());

  hm map3(map);
  ASSERT_EQUAL(map3.size(), 3); 
  ASSERT(map3.find(1) != map3.end());
  ASSERT(map3.find(4) != map3.end());
  ASSERT(map3.find(6) != map3.end());
  ASSERT(map3.find(8) == map3.end());

  /* this block kills my system!!!
  if (true) {
    hm::iterator it1 = map2.end();
    hm::iterator it2 = map3.end();
    ASSERT(it1 != it2); 
    } */
}

cute::suite CuckooSuite() {
  cute::suite s;
  s.push_back(CUTE(TestCuckoo));
  s.push_back(CUTE(TestCuckooCopiing));
  return s;
}
