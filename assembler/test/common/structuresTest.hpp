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

template <class hm>
void TestStructure() {
  struct rusage ru1;
  getrusage(RUSAGE_SELF, &ru1);
  double t1 = ru1.ru_utime.tv_sec + ((float)ru1.ru_utime.tv_usec/1000000);

  srand(42);
  hm map;

  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    map.insert(std::make_pair(t, 42));
  }

  struct rusage ru2;
  getrusage(RUSAGE_SELF, &ru2);
  double t2 = ru2.ru_utime.tv_sec + ((float)ru2.ru_utime.tv_usec/1000000);

  size_t size = map.size();

  typename hm::iterator hmi;
  for (int i = 0; i < 100000; ++i) {
    int t = rand();
    hmi = map.find(t);
  }

  struct rusage ru3;
  getrusage(RUSAGE_SELF, &ru3);
  double t3 = ru3.ru_utime.tv_sec + ((float)ru3.ru_utime.tv_usec/1000000);

  std::cout << "Number of elements was " << size << std::endl;
  std::cout << "Time for insert is " << (t2 - t1) << std::endl;
  std::cout << "Time for find is " << (t3 - t2) << std::endl;
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
