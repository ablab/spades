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
//, 5, 10000, 100, 8, 5> hm;

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
  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    map.erase(t);
  }

  map.clear();
  map[1] = 4;
  map[4] = 1;
  map[6];
  hm map2;
  map2 = map;
  ASSERT_EQUAL(map2.size(), 3);
  hm map3(map);
  ASSERT_EQUAL(map3.size(), 3);
}

cute::suite CuckooSuite() {
  cute::suite s;
  s.push_back(CUTE(TestCuckoo));
  return s;
}
