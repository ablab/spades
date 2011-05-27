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
#include "../memory.hpp"

template <class T>
void TestStructure(T& map, timeval& tim, double& vm, double& rss) {
  const long Number = 1000000;
  srand(42);

  double t1 = tim.tv_sec + ((float)tim.tv_usec/1000000);

  for (int i = 0; i < Number; ++i) {
    int t = rand();
    map.insert(std::make_pair(t, 42));
  }

  gettimeofday(&tim, NULL);
  double t2 = tim.tv_sec + ((float)tim.tv_usec/1000000);

  size_t size = map.size();

  typename T::iterator it;
  for (int i = 0; i < Number; ++i) {
    int t = rand();
    it = map.find(t);
  }

  gettimeofday(&tim, NULL);
  double t3 = tim.tv_sec + ((float)tim.tv_usec/1000000);
  double vm_ = 0;
  double rss_ = 0;
  process_mem_usage(vm_, rss_);
  
  std::cout << "Memory after inserting " << size << 
    " elements: VM = " << (vm_ - vm) << " KB; RSS = " << 
    (rss_ - rss) << " KB.\n";
  std::cout << "Time used for " << Number << " inserts = " << 
    (t2 - t1) << " sec; for " << Number << " finds: " << 
    (t3 - t2) << " sec.\n"; 
}

void TestAllStructures() {
  timeval tim;
  double vm = 0, rss = 0;

  std::cout << "Testing standart map: \n";
  gettimeofday(&tim, NULL);
  process_mem_usage(vm, rss);
  std::map<int, int> s1;
  TestStructure<std::map<int, int> >(s1, tim, vm, rss);

  std::cout << "Testing ext/hash_map: \n";
  gettimeofday(&tim, NULL);
  process_mem_usage(vm, rss);
  __gnu_cxx::hash_map<int, int> s2;
  TestStructure<__gnu_cxx::hash_map<int, int> >(s2, tim, vm, rss);

  std::cout << "Testing tr1/unordered_map: \n";
  gettimeofday(&tim, NULL);
  process_mem_usage(vm, rss);
  std::tr1::unordered_map<int, int> s3;
  TestStructure<std::tr1::unordered_map<int, int> >(s3, tim, vm, rss);

  std::cout << "Testing google/sparse_hash_map: \n";
  gettimeofday(&tim, NULL);
  process_mem_usage(vm, rss);
  google::sparse_hash_map<int, int> s4;
  TestStructure<google::sparse_hash_map<int, int> >(s4, tim, vm, rss);

  std::cout << "Testing cuckoo: \n";
  gettimeofday(&tim, NULL);
  process_mem_usage(vm, rss);
  cuckoo<int, int, Hasher, std::equal_to<int> > s5;
  TestStructure<cuckoo<int, int, Hasher, std::equal_to<int> > >
    (s5, tim, vm, rss);

  /*std::cout << "Testing trie: \n";
  gettimeofday(&tim, NULL);
  process_mem_usage(vm, rss);
  trie<int, int> s6;
  TestStructure<trie<int, int> >(s6, tim, vm, rss); */
}

cute::suite StructuresSuite(){
  cute::suite s;
  s.push_back(CUTE(TestAllStructures));
  return s;
}
