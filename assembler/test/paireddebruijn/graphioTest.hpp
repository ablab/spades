#include "graphio.hpp"
#include "common.hpp"
#include "pairedGraph.hpp"
#include "cute.h"

void TestIntOutput() {
	DataPrinter dp("test/data/oppa.txt");
	dp.outputInt(1);
	dp.close();
	DataReader dr("test/data/oppa.txt");
	int a;
	dr.readInt(a);
	dr.close();
	ASSERT_EQUAL(1, a);
}

void TestTwoIntOutput() {
	DataPrinter dp("test/data/oppa.txt");
	dp.outputInt(4);
	dp.outputInt(5);
	dp.close();
	DataReader dr("test/data/oppa.txt");
	int a, b;
	dr.readInt(a);
	dr.readInt(b);
	dr.close();
	ASSERT_EQUAL(4, a);
	ASSERT_EQUAL(5, b);
}

void TestIntArrayOutput() {
	int oppa[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 0 };
	DataPrinter dp("test/data/oppa.txt");
	dp.outputIntArray(oppa, 10);
	dp.close();
	DataReader dr("test/data/oppa.txt");
	int result[10];
	dr.readIntArray(result, 10);
	dr.close();
	for(int i = 0; i < 10; i++)
		ASSERT_EQUAL(oppa[i], result[i]);
}

void TestTwoSameIntArrayOutput() {
	int oppa[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 0 };
	DataPrinter dp("test/data/oppa.txt");
	dp.outputIntArray(oppa, 10);
	dp.outputIntArray(oppa, 10);
	dp.close();
	DataReader dr("test/data/oppa.txt");
	int result[10];
	dr.readIntArray(result, 10);
	for(int i = 0; i < 10; i++)
		ASSERT_EQUAL(oppa[i], result[i]);
	dr.readIntArray(result, 10);
	for(int i = 0; i < 10; i++)
		ASSERT_EQUAL(oppa[i], result[i]);
	dr.close();
}

void TestTwoIntArrayOutput() {
	int oppa1[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 0 };
	int oppa2[10] = { 11, 12, 13, 14, 15, 16, 17, 18, 19, 10 };
	DataPrinter dp("test/data/oppa.txt");
	dp.outputIntArray(oppa1, 10);
	dp.outputIntArray(oppa2, 10);
	dp.close();
	DataReader dr("test/data/oppa.txt");
	int result[10];
	dr.readIntArray(result, 10);
	for(int i = 0; i < 10; i++)
		ASSERT_EQUAL(oppa1[i], result[i]);
	dr.readIntArray(result, 10);
	for(int i = 0; i < 10; i++)
		ASSERT_EQUAL(oppa2[i], result[i]);
	dr.close();
}

void Test2DimentionIntArrayOutput() {
	int oppa[3][4] = { {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12} };
	DataPrinter dp("test/data/oppa.txt");
	dp.outputIntArray((int*)oppa, 3, 4);
	dp.close();
	DataReader dr("test/data/oppa.txt");
	int result[3][4];
	dr.readIntArray((int*)result, 3, 4);
	dr.close();
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 4; j++)
			ASSERT_EQUAL(oppa[i][j], result[i][j]);
}

void TestSequenceOutput() {
	Sequence *sequence = new Sequence("ACGACTG");
	DataPrinter dp("test/data/oppa.txt");
	dp.outputSequence(sequence);
	dp.close();
	DataReader dr("test/data/oppa.txt");
	Sequence *result;
	dr.readSequence(result);
	dr.close();
	ASSERT_EQUAL(*sequence, *result);
}

void TestTwoSequenceOutput() {
	Sequence *sequence1 = new Sequence("ACGACTG");
	Sequence *sequence2 = new Sequence("TATATTATATTA");
	DataPrinter dp("test/data/oppa.txt");
	dp.outputSequence(sequence1);
	dp.outputSequence(sequence2);
	dp.close();
	DataReader dr("test/data/oppa.txt");
	Sequence *result;
	dr.readSequence(result);
	ASSERT_EQUAL(*sequence1, *result);
	dr.readSequence(result);
	ASSERT_EQUAL(*sequence2, *result);
	dr.close();
}

void TestLongEdgesMap() {
	longEdgesMap edges;
	Edge *e = new Edge(new Sequence(""), new Sequence(""), 1, 2, 3, 0);
	edges.insert(make_pair(0, e));
	DataPrinter dp("test/data/oppa.txt");
	dp.outputLongEdgesMap(edges);
	dp.close();
	DataReader dr("test/data/oppa.txt");
	edges.clear();
	dr.readLongEdgesMap(edges);
	ASSERT_EQUAL(1u, edges.size());
	ASSERT_EQUAL(*e->upper, *edges[0]->upper);
	ASSERT_EQUAL(*e->lower, *edges[0]->lower);
	ASSERT_EQUAL(e->EdgeId, edges[0]->EdgeId);
	ASSERT_EQUAL(e->FromVertex, edges[0]->FromVertex);
	ASSERT_EQUAL(e->length, edges[0]->length);
	ASSERT_EQUAL(e->ToVertex, edges[0]->ToVertex);
	dr.close();
}

void TestSave() {
	PairedGraph *graph = new PairedGraph();
	longEdgesMap *edges = new longEdgesMap() ;
	int a;
	int b;
	save("test/data/oppa.txt", *graph, *edges, a, b);
}

cute::suite GraphioSuite() {
	cute::suite s;
	s.push_back(CUTE(TestIntOutput));
	s.push_back(CUTE(TestTwoIntOutput));
	s.push_back(CUTE(TestIntArrayOutput));
	s.push_back(CUTE(TestTwoSameIntArrayOutput));
	s.push_back(CUTE(TestTwoIntArrayOutput));
	s.push_back(CUTE(Test2DimentionIntArrayOutput));
	s.push_back(CUTE(TestSequenceOutput));
	s.push_back(CUTE(TestTwoSequenceOutput));
	s.push_back(CUTE(TestSave));
	s.push_back(CUTE(TestLongEdgesMap));
	return s;
}

