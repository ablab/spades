#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "trie.hpp"

void TestTrie() {
  srand(42);
  trie<int, int> map; 
  int a[20] = {10, 15, 20, 25, 30, 35, 40, 45, 10, 15, 
               100, 1, 2, 3, 4, 5, 255, 10, 255, 0};
  for (int i = 0; i < 20; ++i) {
    map.insert(std::make_pair(a[i], 42));
  }
  ASSERT_EQUAL(map.size(), 16);
  
  map.clear();
  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    map.insert(std::make_pair(t, 42));
  }
  size_t size = map.size();
  ASSERT_EQUAL(map.size(), 99994);
 
  trie<int, int>::iterator it;
  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    it = map.find(t);
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
  trie<int, int> map2;
  map2 = map;
  ASSERT_EQUAL(map2.size(), 3);
  trie<int, int> map3(map);
  ASSERT_EQUAL(map3.size(), 3); 
}

cute::suite TrieSuite(){
  cute::suite s;
  s.push_back(CUTE(TestTrie));
  return s;
}
