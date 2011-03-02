#include "cute.h"
#include "seq.hpp"
#include "nucl.hpp"

void TestSelector() {
	ASSERT(nucl(Seq<5>("ACGTA")[2]) == 'G');
}

void TestShift() {
	ASSERT(Seq<5>::equal_to()(Seq<5>("ACGTA") << unnucl('T'), Seq<5>("CGTAT")));
}

cute::suite SeqSuite(){
	cute::suite s;
	s.push_back(CUTE(TestSelector));
	s.push_back(CUTE(TestShift));
	return s;
}

