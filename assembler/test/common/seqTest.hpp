#include "cute.h"
#include "seq.hpp"
#include "sequence.hpp"
#include "nucl.hpp"
#include <string>

void TestSeqSelector() {
	ASSERT_EQUAL('G', nucl(Seq<10>("ACGTACGTAC")[2]));
}

void TestSeqShifts() {
	Seq<10> s("ACGTACGTAC");
	ASSERT_EQUAL("CGTACGTACA", (s << 'A').str());
	ASSERT_EQUAL("CGTACGTACT", (s << 'T').str());
	ASSERT_EQUAL("AACGTACGTA", (s >> 'A').str());
	ASSERT_EQUAL("TACGTACGTA", (s >> 'T').str());
}

void TestSeqStr() {
	Seq<10> s("ACGTACGTAC");
	ASSERT_EQUAL("ACGTACGTAC", s.str());
}

void TestSeqHeadAndTail() {
	Seq<10> s("ACGTACGTAC");
	ASSERT_EQUAL("CGTACGTAC", s.tail<9>().str());
	ASSERT_EQUAL("ACGTACGTA", s.head<9>().str());
}

void TestSeqFromBiggerSeq() {
	Seq<10> s("ACGTACGTAC");
	ASSERT_EQUAL("ACGTA", Seq<5>(s).str());
}

void TestSeqFromType() {
	Sequence s("ACGTACGTAC");
	ASSERT_EQUAL("ACGTA", Seq<5>(s).str());
	ASSERT_EQUAL("GTACG", Seq<5>(s, 2).str());
}

void TestSeqFromCharArray() {
	std::string s = "ACGTACGTAC";
	ASSERT_EQUAL(Seq<10>(s.c_str()).str(), "ACGTACGTAC");
}

void TestSeqReverseComplement() {
	Seq<10> s("ACGTACGTAC");
	ASSERT((!s).str() == "GTACGTACGT");
}

cute::suite SeqSuite(){
	cute::suite s;
	s.push_back(CUTE(TestSeqSelector));
	s.push_back(CUTE(TestSeqStr));
	s.push_back(CUTE(TestSeqShifts));
	s.push_back(CUTE(TestSeqFromCharArray));
	s.push_back(CUTE(TestSeqFromBiggerSeq));
	s.push_back(CUTE(TestSeqFromType));
	s.push_back(CUTE(TestSeqHeadAndTail));
	s.push_back(CUTE(TestSeqReverseComplement));
	return s;
}

