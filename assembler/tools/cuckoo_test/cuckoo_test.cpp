#include <iostream>
#include <cstdlib>
#include "seq_filter.hpp"
#include <map>
#include <ext/hash_map>
#include <tr1/unordered_map>
#include <google/sparse_hash_map>
#include "cuckoo.hpp"

#define K 31

const size_t N = 10;
const size_t ml0 = 100;

//typedef cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
//               Seq<K>::equal_to, 4, 1000, 100, 6, 5> hm5; 

int main(int argc, char** argv) {
  std::string filename = "";
  size_t L = 1;
  bool stat_on = true;
  bool console_on = true;
  bool find_on = false;
  size_t max_loop = 10;
  size_t max_loop_limit = 1000;
  size_t loop_step = 10;

	if (argc == 2) {
		filename = argv[1];
  } else {
		std::cout << "Usage: ./cuckoo_test <filename>" << std::endl;
    return 0;
  }

  for (; max_loop <= max_loop_limit; max_loop += loop_step) {
    seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
      Seq<K>::equal_to> >::filter(filename, L, stat_on, console_on, 
                                  find_on, max_loop);
  }

  return 0;
}
