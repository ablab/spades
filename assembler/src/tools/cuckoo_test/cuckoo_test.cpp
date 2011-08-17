#include <iostream>
#include <cstdlib>
#include "seq_filter_stat.hpp"
#include <map>
#include <ext/hash_map>
#include <tr1/unordered_map>
#include <google/sparse_hash_map>
#include "cuckoo.hpp"

#define K 31

int main(int argc, char** argv) {
  std::string filename = "";
  size_t d = 0;
  size_t lbs = 0;
  double step = 0;
  size_t mld = 1;
  std::string label = "";

	if (argc == 7) {
		filename = argv[1];
    d = atoi(argv[2]);
    lbs = atoi(argv[3]);
    step = atoi(argv[4]) / 10.;
    mld = atoi(argv[5]);
    label = argv[6];
  } else {
		std::cout << "Usage: ./cuckoo_test <filename> <d> " <<
      "<log_2(bucket_size)> <step*10> <max_loop_denom> <label>" << std::endl;
    return 0;
  }

  seq_filter_stat<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to> >::filter(filename, d, lbs, step, mld, label);

  return 0;
}
