#include "cute.h"
#include "sequence.hpp"
#include "nucl.hpp"
#include <string>

void TestSelector() {
	ASSERT_EQUAL('G', nucl(Sequence("ACGTACGTAC")[2]));
}

void TestSum() {
	ASSERT_EQUAL("ACGTTGCA", (Sequence("ACGT") + Sequence("TGCA")).str());
}

void TestStr() {
	ASSERT_EQUAL("ACGTACGTAC", Sequence("ACGTACGTAC").str());
}

void TestRefCount() {
	Sequence s("AAAAAAA");
	Sequence s2(s);
	Sequence s3 = !s;
	Sequence s4 = s;
	Sequence ss = s.Subseq(3);
	ASSERTM(s.str() + s2.str() + s3.str() + s4.str()  + ss.str(), true);
}

void TestRefCount2() {
	Sequence *s = new Sequence("AAAAAAA");
	Sequence *s2 = new Sequence(*s);
	Sequence *s3 = new Sequence(!(*s));
	Sequence *ss = new Sequence(s->Subseq(3));
	ASSERTM(s->str() + s2->str() + s3->str() + ss->str(), true);
	delete s;
	delete s2;
	delete s3;
	delete ss;
}


cute::suite SequenceSuite(){
	cute::suite s;
	s.push_back(CUTE(TestSelector));
	s.push_back(CUTE(TestSum));
	s.push_back(CUTE(TestStr));
	s.push_back(CUTE(TestRefCount));
	s.push_back(CUTE(TestRefCount2));
	return s;
}

