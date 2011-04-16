#include "cute.h"
#include "sequence.hpp"
#include "nucl.hpp"
#include <string>
#include "../memory.hpp"
#include <ctime>

void TestSequenceSelector() {
	Sequence s("TTATTAGGGAT");
	ASSERT_EQUAL('G', nucl(Sequence("ACGTACGTAC")[2]));
	ASSERT_EQUAL('A', nucl(Sequence("A")[0]));
}

void TestZeroLengthSequence() {
	Sequence s("");
	ASSERT_EQUAL(0, s.size());
}
void TestSequenceSum() {
	ASSERT_EQUAL("ACG", (Sequence("A") + Sequence("CG")).str());
	ASSERT_EQUAL("ACGTTGCA", (Sequence("ACGT") + Sequence("TGCA")).str());
	ASSERT_EQUAL("ACGTACGTTGCATGCA", (Sequence("ACGTACGT") + Sequence("TGCATGCA")).str());
}

void TestSequenceStr() {
	ASSERT_EQUAL("ACGTACGTAC", Sequence("ACGTACGTAC").str());
	ASSERT_EQUAL("ACG", Sequence("ACG").str());
}

void TestSequenceReverseComplement() {
	Sequence s = Sequence("AACCGGTTAA");
	ASSERT_EQUAL("TTAACCGGTT", (!s).str());
	Sequence s2 = Sequence("ACG");
	ASSERT_EQUAL("CGT", (!s2).str());
}

void TestSequenceRefCount() {
	Sequence s("AAAAAAA");
	Sequence s2(s);
	Sequence s3 = !s;
	Sequence s4 = s;
	Sequence ss = s.Subseq(3);
	ASSERTM(s.str() + s2.str() + s3.str() + s4.str() + ss.str(), true);
}

void TestSequenceRefCount2() {
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

void TestSequenceMemory() {
	time_t now = time(NULL);
	int N = 100000;
	int SIZE = 300;
	vector<Sequence*> vs(N);
	double vm1, rss1;
	double vm2, rss2;
	double vm3, rss3;
	process_mem_usage(vm1, rss1);
	for (int i = 0; i < N; ++i) {
		string s(SIZE,'-');
		for (int j = 0; j < SIZE; ++j) {
			s[j] = nucl(rand() % 4);
		}
		vs[i] = new Sequence(s);
	}
	process_mem_usage(vm2, rss2);
	cout << "Memory after creation for " <<  N << " Sequences of size " << SIZE << ": VM = "<< (vm2 - vm1) << " KB., RSS = "<< (rss2 - rss1) << " KB." << endl;
	for (int i = 0; i < N; ++i) {
		delete vs[i];
		vs[i] = 0;
	}
	process_mem_usage(vm3, rss3);
	cout << "Memory after deletion for " <<  N << " Sequences of size " << SIZE << ": VM = "<< (vm3 - vm1) << " KB., RSS = "<< (rss3 - rss1) << " KB." << endl;
	cout << "Time: " <<  time(NULL) - now << endl;
}

cute::suite SequenceSuite(){
	cute::suite s;
	s.push_back(CUTE(TestSequenceSelector));
	s.push_back(CUTE(TestZeroLengthSequence));
	s.push_back(CUTE(TestSequenceSum));
	s.push_back(CUTE(TestSequenceStr));
	s.push_back(CUTE(TestSequenceReverseComplement));
	s.push_back(CUTE(TestSequenceRefCount));
	s.push_back(CUTE(TestSequenceRefCount2));
//	s.push_back(CUTE(TestSequenceMemory));
	return s;
}

