#include "read.hpp"
#include "seq.hpp"
#include "seq_filter.hpp"

void TestSeqFilter() {
	Read r1 = Read("1st", "ATGCATGCATGC", "++++++++++");
	Read r2 = Read("2nd", "AAAAAAAAAAAC", "++++++++++");
	Read r3 = Read("3rd", "TGTGTGTGTGTG", "++++++++++");
	std::vector<Read> reads;
	reads.push_back(r1);
	reads.push_back(r2);
	reads.push_back(r3);
	//4-mers filtration
	std::vector<Seq<4> > seqs4;
	seqs4 = seq_filter<4>::filter(reads, 0);
	ASSERT_EQUAL(seqs4.size(), 8);
	seqs4 = seq_filter<4>::filter(reads, 1);
	ASSERT_EQUAL(seqs4.size(), 7);
	seqs4 = seq_filter<4>::filter(reads, 2);
	ASSERT_EQUAL(seqs4.size(), 4);
	seqs4 = seq_filter<4>::filter(reads, 5);
	ASSERT_EQUAL(seqs4.size(), 1);
	seqs4 = seq_filter<4>::filter(reads, 8);
	ASSERT_EQUAL(seqs4.size(), 0);
	//7-mers filtration
	std::vector<Seq<7> > seqs7;
	seqs7 = seq_filter<7>::filter(reads, 0);
	ASSERT_EQUAL(seqs7.size(), 8);
	seqs7 = seq_filter<7>::filter(reads, 1);
	ASSERT_EQUAL(seqs7.size(), 5);
	seqs7 = seq_filter<7>::filter(reads, 2);
	ASSERT_EQUAL(seqs7.size(), 3);
}

cute::suite SeqFilterSuite(){
	cute::suite s;
	s.push_back(CUTE(TestSeqFilter));
	return s;
}
