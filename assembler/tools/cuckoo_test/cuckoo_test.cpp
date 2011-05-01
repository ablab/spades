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
const size_t ml1 = 200;
const size_t ml2 = 300;
const size_t ml3 = 400;
const size_t ml4 = 500;
const size_t ml5 = 600;
const size_t ml6 = 700;
const size_t ml7 = 800;
const size_t ml8 = 900;
const size_t ml9 = 1000;

typedef cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
               Seq<K>::equal_to, 4, 1000, 100, 6, 5> hm5; 

int main(int argc, char** argv) {
  std::string filename = "";
  size_t L = 1;
  bool stat = true;

	if (argc == 2) {
		filename = argv[1];
  } else {
		std::cout << "Usage: ./cuckoo_test <filename>" << std::endl;
    return 0;
  }

  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml0, 6, 5> >::filter
    (filename, L, stat);
  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml1, 6, 5> >::filter
    (filename, L, stat);
  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml2, 6, 5> >::filter
    (filename, L, stat);
  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml3, 6, 5> >::filter
    (filename, L, stat);
  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml4, 6, 5> >::filter
    (filename, L, stat);
  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml5, 6, 5> >::filter
    (filename, L, stat);
  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml6, 6, 5> >::filter
    (filename, L, stat);
  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml7, 6, 5> >::filter
    (filename, L, stat);
  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml8, 6, 5> >::filter
    (filename, L, stat);
  seq_filter<K, cuckoo<Seq<K>, size_t, Seq<K>::multiple_hash, 
    Seq<K>::equal_to, 4, 100, ml9, 6, 5> >::filter
    (filename, L, stat);

  return 0;
}
