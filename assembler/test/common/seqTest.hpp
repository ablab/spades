#include "cute.h"
#include "seq.hpp"
#include "nucl.hpp"

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


void TestSeqReverseComplement() {
	Seq<10> s("ACGTACGTAC");
	ASSERT((!s).str() == "GTACGTACGT");
}

cute::suite SeqSuite(){
	cute::suite s;
	s.push_back(CUTE(TestSeqSelector));
	s.push_back(CUTE(TestSeqStr));
	s.push_back(CUTE(TestSeqShifts));
	s.push_back(CUTE(TestSeqHeadAndTail));
	s.push_back(CUTE(TestSeqReverseComplement));
	return s;
}

