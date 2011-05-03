#include "read.hpp"
#include "seq.hpp"
#include "seq_filter.hpp"
#include "cuckoo.hpp"

typedef cuckoo<Seq<4>, size_t, Seq<4>::multiple_hash, 
               Seq<4>::equal_to> hm4; 
typedef cuckoo<Seq<7>, size_t, Seq<7>::multiple_hash, 
               Seq<7>::equal_to> hm7; 

void TestSeqFilter() {
  std::string in = "../../test/data/s_test.fastq.gz";
  //4-mers filtration
  std::vector<Seq<4> > seqs4;
  seqs4 = seq_filter<4, hm4>::filter(in, 0, false, false);
  ASSERT_EQUAL(seqs4.size(), 8);
  seqs4 = seq_filter<4, hm4>::filter(in, 1, false, false);
  ASSERT_EQUAL(seqs4.size(), 7);
  seqs4 = seq_filter<4, hm4>::filter(in, 2, false, false);
  ASSERT_EQUAL(seqs4.size(), 4);
  seqs4 = seq_filter<4, hm4>::filter(in, 5, false, false);
  ASSERT_EQUAL(seqs4.size(), 1);
  seqs4 = seq_filter<4, hm4>::filter(in, 8, false, false);
  ASSERT_EQUAL(seqs4.size(), 0);
  //7-mers filtration
  std::vector<Seq<7> > seqs7;
  seqs7 = seq_filter<7, hm7>::filter(in, 0, false, false);
  ASSERT_EQUAL(seqs7.size(), 8);
  seqs7 = seq_filter<7, hm7>::filter(in, 1, false, false);
  ASSERT_EQUAL(seqs7.size(), 5);
  seqs7 = seq_filter<7, hm7>::filter(in, 2, false, false);
  ASSERT_EQUAL(seqs7.size(), 3);
}

cute::suite SeqFilterSuite(){
  cute::suite s;
  s.push_back(CUTE(TestSeqFilter));
  return s;
}
