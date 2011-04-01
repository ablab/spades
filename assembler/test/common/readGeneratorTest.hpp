#ifndef READGENERATORTEST_HPP_
#define READGENERATORTEST_HPP_
#include "read_generator.hpp"
#include "cute.h"

void TestAllA() {
	string s = "AAAAAAAAAAAAAAAAAAAAAA";
	ReadGenerator<5, 1, int> r(s, 1);
	strobe_read<5, 1> read;
	int cnt = 0;
	while(!r.eof()) {
		r >> read;
		ASSERT_EQUAL(Seq<5>("AAAAA"), read[0]);
		cnt++;
	}
}

void TestAllAPaired() {
	string s = "AAAAAAAAAAAAAAAAAAAAAAA";
	ReadGenerator<5, 2, int> r(s, 1);
	strobe_read<5, 2> read;
	int cnt = 0;
	while(!r.eof()) {
		r >> read;
		ASSERT_EQUAL(Seq<5>("AAAAA"), read[0]);
		ASSERT_EQUAL(Seq<5>("AAAAA"), read[1]);
		cnt++;
	}
}

void TestAllAPairedIL4() {
	string s = "AAAAAAAAAAAAAAAAAAAAAAA";
	ReadGenerator<5, 2, int> r(s, 1, 4);
	strobe_read<5, 2> read;
	int cnt = 0;
	while(!r.eof()) {
		r >> read;
		ASSERT_EQUAL(Seq<5>("AAAAA"), read[0]);
		ASSERT_EQUAL(Seq<5>("AAAAA"), read[1]);
		cnt++;
	}
}

cute::suite ReadGeneratorSuite(){
	cute::suite s;
	s.push_back(CUTE(TestAllA));
	s.push_back(CUTE(TestAllAPaired));
	s.push_back(CUTE(TestAllAPairedIL4));
	return s;
}

#endif /* READGENERATORTEST_HPP_ */
