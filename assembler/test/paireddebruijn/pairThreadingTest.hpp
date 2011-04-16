#ifndef PAIRTHREADINGTEST_HPP_
#define PAIRTHREADINGTEST_HPP_
#include "graphSimplification.hpp"
#include "cute.h"

using namespace paired_assembler;
string sameLetterString(int number, char a) {
	char *tmp = new char[number + 1];
	for (int i = 0; i < number; i++)
		tmp[i] = a;
	tmp[number] = 0;
	string result(tmp);
	delete[] tmp;
	return result;
}

void TestThreadAllA() {
	k = 6;
	l = 6;
	readLength = 10;
	insertLength = 5;
	PairedGraph *g = new PairedGraph();
//	longEdgesMap edges;
	Sequence *s = new Sequence(sameLetterString(15, 'A'));
	for (int i = 0; i < 10; i++) {
		Edge * e = new Edge(s, s, i, i + 1, 15 - k + 1, i);
		g->longEdges.insert(make_pair(i, e));
	}
	g->recreateVerticesInfo(g->VertexCount, g->longEdges);
	PairThreader pt(*g);
	vector<pair<int, Edge *> > v = pt.threadLower(g->longEdges[0]);
	ASSERT_EQUAL(2u, v.size());
	ASSERT_EQUAL(0, v[0].first);
	ASSERT_EQUAL(1, v[0].second->EdgeId);
	ASSERT_EQUAL(10, v[1].first);
	ASSERT_EQUAL(2, v[1].second->EdgeId);
}

cute::suite PairThreadingSuite() {
	cute::suite s;
	s.push_back(CUTE(TestThreadAllA));
	return s;
}

#endif /* PAIRTHREADINGTEST_HPP_ */
