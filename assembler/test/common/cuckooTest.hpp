/* test for cuckoo.hpp */
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "cuckoo.hpp"

struct Hasher {
  size_t operator()(int value, int hash_num) {
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
  hm::iterator hmi;
  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    hmi = map.find(t);
  }
  ASSERT_EQUAL(size, map.size());
  size_t n = 0;
  for (hm::iterator it = map.begin(); it != map.end(); ++it) {
    ++n;
  }
  ASSERT_EQUAL(size, n);

  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    map.erase(t);
  }
  map.clear();
  ASSERT_EQUAL(map.size(), 0);
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
  hm map3(map);
  ASSERT_EQUAL(map3.size(), 3); 
}

cute::suite CuckooSuite() {
  cute::suite s;
  s.push_back(CUTE(TestCuckoo));
  s.push_back(CUTE(TestCuckooCopiing));
  return s;
}
