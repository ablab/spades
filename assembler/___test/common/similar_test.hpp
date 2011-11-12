#include "sequence.hpp"
#include "cute.h"
#include "vector"

void TestREqualSequences() {
	Sequence a("A");
	Sequence a1("A");
	Sequence b("C");
	ASSERT_EQUAL(0, a.rightSimilar(b, 1));
	ASSERT_EQUAL(1, a.rightSimilar(a1, 1));
	ASSERT_EQUAL(1, a.rightSimilar(a, 1));
}

void TestRSimilarInclude() {
	vector<Sequence*> s;
	s.push_back(new Sequence("AAAAAA"));
	s.push_back(new Sequence("CCCCCC"));
	s.push_back(new Sequence("GGGGGG"));
	s.push_back(new Sequence("TTTTTT"));
	Sequence sA5C("AAAAAACCCCCC");
	Sequence sG5T("GGGGGGTTTTTT");
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			ASSERT_EQUAL(i == j, s[i]->rightSimilar(*s[j], 1));
	ASSERT_EQUAL(false, sA5C.rightSimilar(*s[0], 1));
	ASSERT_EQUAL(true, s[0]->rightSimilar(sA5C, 1));
	ASSERT_EQUAL(false, sG5T.rightSimilar(*s[2], 1));
	ASSERT_EQUAL(true, s[2]->rightSimilar(sG5T, 1));
	ASSERT_EQUAL(true, sA5C.rightSimilar(*s[1], 1));
	ASSERT_EQUAL(false, s[1]->rightSimilar(sA5C, 1));
	ASSERT_EQUAL(true, sG5T.rightSimilar(*s[3], 1));
	ASSERT_EQUAL(false, s[3]->rightSimilar(sG5T, 1));
	for (size_t i = 0; i < s.size(); i++) {
		delete s[i];
	}
}

void TestRShortIntersection() {
	Sequence sGA3("GGGGGGAAA");
	Sequence sA3C("AAACCCCCC");
	for (int i = 1; i < 7; i++) {
		ASSERT_EQUAL(i <= 3, sGA3.rightSimilar(sA3C, i));
	}
}

void TestLEqualSequences() {
	Sequence a("A");
	Sequence a1("A");
	Sequence b("C");
	ASSERT_EQUAL(0, a.leftSimilar(b, 1));
	ASSERT_EQUAL(1, a.leftSimilar(a1, 1));
	ASSERT_EQUAL(1, a.leftSimilar(a, 1));
}

void TestLSimilarInclude() {
	vector<Sequence*> s;
	s.push_back(new Sequence("AAAAAA"));
	s.push_back(new Sequence("CCCCCC"));
	s.push_back(new Sequence("GGGGGG"));
	s.push_back(new Sequence("TTTTTT"));
	Sequence sA5C("AAAAAACCCCCC");
	Sequence sG5T("GGGGGGTTTTTT");
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			ASSERT_EQUAL(i == j, s[i]->leftSimilar(*s[j], 1));
		}
	}
	ASSERT_EQUAL(true, sA5C.leftSimilar(*s[0], 1));
	ASSERT_EQUAL(false, s[0]->leftSimilar(sA5C, 1));
	ASSERT_EQUAL(true, sG5T.leftSimilar(*s[2], 1));
	ASSERT_EQUAL(false, s[2]->leftSimilar(sG5T, 1));
	ASSERT_EQUAL(false, sA5C.leftSimilar(*s[1], 1));
	ASSERT_EQUAL(true, s[1]->leftSimilar(sA5C, 1));
	ASSERT_EQUAL(false, sG5T.leftSimilar(*s[3], 1));
	ASSERT_EQUAL(true, s[3]->leftSimilar(sG5T, 1));
	for (size_t i = 0; i < s.size(); i++) {
		delete s[i];
	}
}

void TestLShortIntersection() {
	Sequence sGA3("GGGGGGAAA");
	Sequence sA3C("AAACCCCCC");
	for (size_t i = 1; i < 7; i++) {
		ASSERT_EQUAL(i <= 3, sA3C.leftSimilar(sGA3, i));
	}
}

cute::suite similarSuite() {
	cute::suite s;
	s.push_back(CUTE(TestREqualSequences));
	s.push_back(CUTE(TestRSimilarInclude));
	s.push_back(CUTE(TestRShortIntersection));
	s.push_back(CUTE(TestLEqualSequences));
	s.push_back(CUTE(TestLSimilarInclude));
	s.push_back(CUTE(TestLShortIntersection));
	return s;
}
