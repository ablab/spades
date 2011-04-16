#include "ireadstream.hpp"
#include "seq.hpp"
#include <iostream>
#include "cuckoo.hpp"

//const std::string filename = "../../src/tools/s_6.first10000_1.fastq.gz";
const std::string filename = "/dev/stdin";
//const std::string output = "/dev/stdout";
const size_t size_ = 5;

/*struct Hasher {
  size_t operator()(Seq<size_> seq, int hash_num) {
    size_t h = 239;
    for (size_t i = 0; i < seq.data_size_; i++) {
      h = ((h << 5) - h) + seq.data_[i];
    }
    size_t h = Seq<size_>::hash(seq);
    unsigned long l = 4 * hash_num + 1;
    return (size_t)(l * h % 1000000007);
  }
  };*/

typedef cuckoo<Seq<size_>, size_t, Seq<size_>::multiple_hash, Seq<size_>::equal_to, 4, 100, 100, 3, 2> hm; 

int main() {
  vector<Read>* reads = ireadstream::readAll(filename, 10000);
  hm map;
  hm::iterator it; 
  //size_t cnt = 5;
  size_t cnt = reads->size();
  for (size_t i = 0; i < cnt; i++) {
    Sequence s = (*reads)[i].getSequence(); 
    Seq<size_> kmer = s.start<size_>();
    it = map.find(kmer);
    if (it == map.end()) {
      map.insert(std::make_pair(kmer, 1));
    } else {
      (*it).second++;
    }
    //std::cout << kmer.str() << " ";
    for (size_t j = size_; j < s.size(); ++j) {
      Seq<size_> next = kmer << s[j];
      kmer = next; 
      it = map.find(kmer);
      if (it == map.end()) {
	map.insert(std::make_pair(kmer, 1));
      } else {
	(*it).second++;
      }
      //std::cout << kmer.str() << " ";
    }
  }
  std::cout << reads->size() << std::endl;
  for (hm::iterator it = map.begin(); it != map.end(); ++it) {
    std::cout << (*it).first.str() << "-" << (*it).second << " "; 
  }
	 std::cout << map.size() << " " << map.length() << std::endl;
  return 0;
}
