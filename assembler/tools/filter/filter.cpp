#include <iostream>
#include <cstdlib>
#include "seq_filter.hpp"
#include <map>
#include <ext/hash_map>
#include <tr1/unordered_map>
#include <google/sparse_hash_map>
#include "cuckoo.hpp"
#include "../test/memory.hpp"

#define K 31

typedef std::map<Seq<K>, size_t, Seq<K>::less> hm1;
typedef __gnu_cxx::hash_map<Seq<K>, size_t, 
                            Seq<K>::hash, Seq<K>::equal_to> hm2;
typedef std::tr1::unordered_map<Seq<K>, size_t, 
                                Seq<K>::hash, Seq<K>::equal_to> hm3;
typedef google::sparse_hash_map<Seq<K>, size_t, 
                                Seq<K>::hash, Seq<K>::equal_to> hm4;
typedef cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
               Seq<K>::equal_to, 4, 1000, 100, 6, 5> hm5; 

int main(int argc, char** argv) {
  std::string filename = "";
  size_t L = 1;
  size_t m_num = 1;
  bool stat = false;
  double vm1 = 0;
  double rss1 = 0;
  process_mem_usage(vm1, rss1);

  if ((argc < 4) || (argc > 5)) {
		std::cout << "Usage: ./filter <filename> <L> <map number> [--stat]\n"
              << "<map number> is map type:\n"
              << "1 - map, 2 - ext/hash_map, 3 - tr1/unordered_map,\n"
              << "4 - google/sparse_hash_map, 5 - cuckoo\n"
              << "Selects k-mer with amount > L\n" 
              << "--stat output stat info without k-mers" << std::endl;
    return 0;
  }
	if (argc >= 4) {
		filename = argv[1];
    L = atoi(argv[2]);
    m_num = atoi(argv[3]);
  }
  if (argc == 5) {
    if (std::string(argv[4]) == "--stat") {
      stat = true;
    } else {
      std::cout << "Wrong option!\n";
      return 0;
    }
  }
  switch (m_num) {
  case 1: 
    seq_filter<K, hm1>::filter(filename, L, stat);
    break;
  case 2: 
    seq_filter<K, hm2>::filter(filename, L, stat);
    break;
  case 3: 
    seq_filter<K, hm3>::filter(filename, L, stat);
    break;
  case 4: 
    seq_filter<K, hm4>::filter(filename, L, stat);
    break;
  case 5: 
    seq_filter<K, hm5>::filter(filename, L, stat);
    break;
  default:
    std::cout << "Map number is incorrect!\n";
  }

  double vm2 = 0;
  double rss2 = 0;
  process_mem_usage(vm2, rss2);
  if (stat) {
    std::cout << "Memory: " << (vm2 - vm1) << std::endl;
  }
  return 0;
}
