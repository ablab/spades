#include "cute.h"
#include "graphConstruction.hpp"
#include "iostream"

using namespace gvis;

string repeatLetter(char c, int num) {
	string res;
	for (int i = 0; i < num; i++)
		res += c;
	return res;
}

string streamToString(istream &s) {
	string res;
	res.assign((std::istreambuf_iterator<char>(s)),
			std::istreambuf_iterator<char>());
	return res;
}

void TestStoreVertexEmpty() {
	resetVertexCount();
	string s;
	stringstream ss(s);
	GraphPrinter<int> g("oppa", ss);
	verticesMap map;
	Sequence seq(repeatLetter('A', k));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	g.output();
	ASSERT_EQUAL("digraph oppa {\n0 [label=" + repeatLetter('A', k - 1) + "]\n}\n", streamToString(ss));
}

void TestStoreVertexOneEntry() {
	resetVertexCount();
	string s;
	stringstream ss(s);
	GraphPrinter<int> g("oppa", ss);
	verticesMap map;
	Sequence seq(repeatLetter('A', l));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	g.output();
	ASSERT_EQUAL("digraph oppa {\n0 [label=" + repeatLetter('A', k - 1) + "]\n}\n", streamToString(ss));
}

void TestStoreVertexTwoEntries() {
	resetVertexCount();
	string s;
	stringstream ss(s);
	GraphPrinter<int> g("oppa", ss);
	verticesMap map;
	Sequence seqA(repeatLetter('A', l));
	Sequence seqB(repeatLetter('C', l));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seqA));
	ASSERT_EQUAL(1, storeVertex(g, map, 1, &seqB));
	g.output();
	string result = "digraph oppa {\n0 [label=" + repeatLetter('A', k - 1)
			+ "]\n1 [label=" + repeatLetter('A', k - 2) + "C" + "]\n}\n";
	ASSERT_EQUAL(result, streamToString(ss));
}

void TestStoreVertexTwoEntriesWithSameKmer() {
	resetVertexCount();
	string s;
	stringstream ss(s);
	GraphPrinter<int> g("oppa", ss);
	verticesMap map;
	Sequence seqA(repeatLetter('A', l));
	Sequence seqB(repeatLetter('C', l));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seqA));
	ASSERT_EQUAL(1, storeVertex(g, map, 0, &seqB));
	g.output();
	string result = "digraph oppa {\n0 [label=" + repeatLetter('A', k - 1)
			+ "]\n1 [label=" + repeatLetter('A', k - 1) + "]\n}\n";
	ASSERT_EQUAL(result, streamToString(ss));
}

void TestCheckUniqueWayRightEmpty() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	int res = checkUniqueWayRight(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayRightSingle() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	VertexPrototype vp(&seqA, 0);
	edges[1].push_back(&vp);
	cout << edges[0].size() << " " << edges[1].size() << " " << edges[2].size()
			<< endl;
	int res = checkUniqueWayRight(edges, 0, &seqA);
	ASSERT_EQUAL(1, res);
}

void TestCheckUniqueWayRightDoublePair() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	Sequence seqAG(repeatLetter('A', l) + "G");
	Sequence seqAC(repeatLetter('A', l) + "C");
	VertexPrototype vpAC(&seqAC, 0);
	VertexPrototype vpAG(&seqAG, 0);
	edges[1].push_back(&vpAC);
	edges[1].push_back(&vpAG);
	int res = checkUniqueWayRight(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayRightDoubleKmer() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	VertexPrototype vp1(&seqA, 0);
	VertexPrototype vp2(&seqA, 1);
	edges[0].push_back(&vp1);
	edges[1].push_back(&vp2);
	int res = checkUniqueWayRight(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayLeftEmpty() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	int res = checkUniqueWayLeft(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayLeftSingle() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	VertexPrototype vp(&seqA, 0);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vp);
	int res = checkUniqueWayLeft(edges, 0, &seqA);
	ASSERT_EQUAL(1, res);
}

void TestCheckUniqueWayLeftDoublePair() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	Sequence seqAG("G" + repeatLetter('A', l));
	Sequence seqAC("C" + repeatLetter('A', l));
	VertexPrototype vpAC(&seqAC, 0);
	VertexPrototype vpAG(&seqAG, 0);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vpAC);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vpAG);
	int res = checkUniqueWayLeft(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayLeftDoubleKmer() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	VertexPrototype vp1(&seqA, 0);
	VertexPrototype vp2(&seqA, 1);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vp1);
	edges[(ll) 2 << (2 * (k - 1))].push_back(&vp2);
	int res = checkUniqueWayLeft(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestGoUniqueWayRightEmpty() {
	edgesMap edges;
	Sequence * seqA = new Sequence(repeatLetter('A', l));
	ll kMer = 0;
	int res = goUniqueWayRight(edges, kMer, seqA);
	ASSERT_EQUAL(0, res);
}

void TestGoUniqueWayRightSingle() {
	edgesMap edges;
	Sequence * seqA = new Sequence(repeatLetter('A', l));
	VertexPrototype vp(seqA, 0);
	edges[1].push_back(&vp);
	ll kMer = 0;
	int res = goUniqueWayRight(edges, kMer, seqA);
	ASSERT_EQUAL(1, res);
	ASSERT_EQUAL(1, kMer);
	Sequence *result = new Sequence(repeatLetter('A', l - 1));
	ASSERT_EQUAL(*result, *seqA);
}

void TestGoUniqueWayRightDoublePair() {
	edgesMap edges;
	Sequence * seqA = new Sequence(repeatLetter('A', l));
	Sequence * seqAG = new Sequence(repeatLetter('A', l) + "G");
	Sequence * seqAC = new Sequence(repeatLetter('A', l) + "C");
	VertexPrototype vpAC(seqAC, 0);
	VertexPrototype vpAG(seqAG, 0);
	edges[1].push_back(&vpAC);
	edges[1].push_back(&vpAG);
	ll kMer = 0;
	int res = goUniqueWayRight(edges, kMer, seqA);
	ASSERT_EQUAL(0, res);
}

void TestGoUniqueWayRightDoubleKmer() {
	edgesMap edges;
	Sequence *seqA = new Sequence(repeatLetter('A', l));
	VertexPrototype vp1(seqA, 0);
	VertexPrototype vp2(seqA, 1);
	edges[0].push_back(&vp1);
	edges[1].push_back(&vp2);
	ll kMer = 0;
	int res = goUniqueWayRight(edges, kMer, seqA);
	ASSERT_EQUAL(0, res);
}

void TestGoUniqueWayLeftEmpty() {
	edgesMap edges;
	Sequence *seqA = new Sequence(repeatLetter('A', l));
	ll kMer = 0;
	int res = goUniqueWayLeft(edges, kMer, seqA);
	ASSERT_EQUAL(0, res);
}

void TestGoUniqueWayLeftSingle() {
	edgesMap edges;
	Sequence *seqA = new Sequence(repeatLetter('A', l));
	VertexPrototype vp(seqA, 0);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vp);
	ll kMer = 0;
	int res = goUniqueWayLeft(edges, kMer, seqA);
	ASSERT_EQUAL(1, res);
	ASSERT_EQUAL((ll) 1 << (2 * (k - 2)), kMer);
	Sequence *result = new Sequence(repeatLetter('A', l - 1));
	ASSERT_EQUAL(*result, *seqA);
}

void TestGoUniqueWayLeftDoublePair() {
	edgesMap edges;
	Sequence *seqA = new Sequence(repeatLetter('A', l));
	Sequence *seqAG = new Sequence("G" + repeatLetter('A', l));
	Sequence *seqAC = new Sequence("C" + repeatLetter('A', l));
	VertexPrototype vpAC(seqAC, 0);
	VertexPrototype vpAG(seqAG, 0);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vpAC);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vpAG);
	ll kMer = 0;
	int res = goUniqueWayLeft(edges, kMer, seqA);
	ASSERT_EQUAL(0, res);
}

void TestGoUniqueWayLeftDoubleKmer() {
	edgesMap edges;
	Sequence *seqA = new Sequence(repeatLetter('A', l));
	VertexPrototype vp1(seqA, 0);
	VertexPrototype vp2(seqA, 1);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vp1);
	edges[(ll) 2 << (2 * (k - 1))].push_back(&vp2);
	ll kMer = 0;
	int res = goUniqueWayLeft(edges, kMer, seqA);
	ASSERT_EQUAL(0, res);
}

cute::suite CheckStoreVertexSuite() {
	cute::suite s;
	s.push_back(CUTE(TestStoreVertexEmpty));
	s.push_back(CUTE(TestStoreVertexOneEntry));
	s.push_back(CUTE(TestStoreVertexTwoEntries));
	s.push_back(CUTE(TestStoreVertexTwoEntriesWithSameKmer));
	return s;
}

cute::suite CheckUniqueWaySuite() {
	cute::suite s;
	s.push_back(CUTE(TestCheckUniqueWayRightEmpty));
	s.push_back(CUTE(TestCheckUniqueWayRightSingle));
	s.push_back(CUTE(TestCheckUniqueWayRightDoublePair));
	s.push_back(CUTE(TestCheckUniqueWayRightDoubleKmer));
	s.push_back(CUTE(TestCheckUniqueWayLeftEmpty));
	s.push_back(CUTE(TestCheckUniqueWayLeftSingle));
	s.push_back(CUTE(TestCheckUniqueWayLeftDoublePair));
	s.push_back(CUTE(TestCheckUniqueWayLeftDoubleKmer));
	return s;
}

cute::suite GoUniqueWaySuite() {
	cute::suite s;
	s.push_back(CUTE(TestGoUniqueWayRightEmpty));
	s.push_back(CUTE(TestGoUniqueWayRightSingle));
	s.push_back(CUTE(TestGoUniqueWayRightDoublePair));
	s.push_back(CUTE(TestGoUniqueWayRightDoubleKmer));
	s.push_back(CUTE(TestGoUniqueWayLeftEmpty));
	s.push_back(CUTE(TestGoUniqueWayLeftSingle));
	s.push_back(CUTE(TestGoUniqueWayLeftDoublePair));
	s.push_back(CUTE(TestGoUniqueWayLeftDoubleKmer));
	return s;
}
