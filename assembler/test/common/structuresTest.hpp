#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <map>
#include <ext/hash_map>
#include <tr1/unordered_map>
#include <google/sparse_hash_map>
#include "cuckoo.hpp"
#include "trie.hpp"
#include "../memory.hpp"

template <class hm>
void TestStructure() {
  const long Number = 1000000;

  timeval tim;
  gettimeofday(&tim, NULL);
  double t1 = tim.tv_sec + ((float)tim.tv_usec/1000000);
  double vm1, rss1;
  process_mem_usage(vm1, rss1);

  srand(42);
  hm map;

  for (int i = 0; i < Number; ++i) {
    int t = rand();
    map.insert(std::make_pair(t, 42));
  }

  gettimeofday(&tim, NULL);
  double t2 = tim.tv_sec + ((float)tim.tv_usec/1000000);
  double vm2, rss2;
  process_mem_usage(vm2, rss2);

  size_t size = map.size();

  typename hm::iterator hmi;
  for (int i = 0; i < Number; ++i) {
    int t = rand();
    hmi = map.find(t);
  }

  gettimeofday(&tim, NULL);
  double t3 = tim.tv_sec + ((float)tim.tv_usec/1000000);

  std::cout << "Memory after inserting " << size << " elements: VM = " << (vm2 - vm1) << " KB; RSS = " << (rss2 - rss1) << " KB.\n";
  std::cout << "Time used for " << Number << " inserts = " << (t2 - t1) << "sec; for " << Number << " finds: " << (t3 - t2) << "sec.\n";
}

void TestStandartMap() {
  TestStructure<std::map<int, int> >();
}

void TestExtHashMap() {
  TestStructure<__gnu_cxx::hash_map<int, int> >();
}

void TestUnorderedMap() {
  TestStructure<std::tr1::unordered_map<int, int> >();
}

void TestSparseHashMap() {
  TestStructure<google::sparse_hash_map<int, int> >();
}

void TestCuckooStructure() {
  TestStructure<cuckoo<int, int, Hasher, std::equal_to<int>, 
    5, 10000, 100, 8, 5> >();
}

//void TestTrieStructure() {
//  TestStructure<trie<int, int> >();
//}

cute::suite StructuresSuite(){
  cute::suite s;
  s.push_back(CUTE(TestStandartMap));
  s.push_back(CUTE(TestUnorderedMap));
  s.push_back(CUTE(TestExtHashMap));
  s.push_back(CUTE(TestSparseHashMap));
  s.push_back(CUTE(TestCuckooStructure));
  //s.push_back(CUTE(TestTrieStructure));
  return s;
}
