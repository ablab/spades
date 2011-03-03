#include "cute.h"
#include "seq.hpp"
#include "nucl.hpp"

void TestSelector() {
	ASSERT_EQUAL('G', nucl(Seq<10>("ACGTACGTAC")[2]));
}

void TestShifts() {
	Seq<10> s("ACGTACGTAC");
	ASSERT_EQUAL("CGTACGTACA", (s << 'A').str());
	ASSERT_EQUAL("CGTACGTACT", (s << 'T').str());
	ASSERT_EQUAL("AACGTACGTA", (s >> 'A').str());
	ASSERT_EQUAL("TACGTACGTA", (s >> 'T').str());
}

void TestStr() {
	Seq<10> s("ACGTACGTAC");
	ASSERT_EQUAL("ACGTACGTAC", s.str());
}

void TestHeadAndTail() {
	Seq<10> s("ACGTACGTAC");
	ASSERT_EQUAL("CGTACGTAC", s.tail<9>().str());
	ASSERT_EQUAL("ACGTACGTA", s.head<9>().str());
}


void TestReverseComplement() {
	Seq<10> s("ACGTACGTAC");
	ASSERT((!s).str() == "GTACGTACGT");
}

cute::suite SeqSuite(){
	cute::suite s;
	s.push_back(CUTE(TestSelector));
	s.push_back(CUTE(TestStr));
	s.push_back(CUTE(TestShifts));
	s.push_back(CUTE(TestHeadAndTail));
	s.push_back(CUTE(TestReverseComplement));
	return s;
}

