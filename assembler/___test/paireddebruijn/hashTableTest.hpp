#include "cute.h"
#include "seq.hpp"
#include "nucl.hpp"
#include "constructHashTable.hpp"
#include "graphio.hpp"

void TestEmptyDecompress() {
	ASSERT_EQUAL(decompress(0, 0), "");
}

void TestSingleNucleotideDecompress() {
	ASSERT_EQUAL("A", decompress(0, 1));
	ASSERT_EQUAL("C", decompress(1, 1));
	ASSERT_EQUAL("G", decompress(2, 1));
	ASSERT_EQUAL("T", decompress(3, 1));
}

void TestDoubleNucleotideDecompress() {
	ASSERT_EQUAL("AA", decompress(0, 2));
	ASSERT_EQUAL("AC", decompress(1, 2));
	ASSERT_EQUAL("AG", decompress(2, 2));
	ASSERT_EQUAL("AT", decompress(3, 2));
	ASSERT_EQUAL("CA", decompress(4 + 0, 2));
	ASSERT_EQUAL("CC", decompress(4 + 1, 2));
	ASSERT_EQUAL("CG", decompress(4 + 2, 2));
	ASSERT_EQUAL("CT", decompress(4 + 3, 2));
	ASSERT_EQUAL("GA", decompress(8 + 0, 2));
	ASSERT_EQUAL("GC", decompress(8 + 1, 2));
	ASSERT_EQUAL("GG", decompress(8 + 2, 2));
	ASSERT_EQUAL("GT", decompress(8 + 3, 2));
	ASSERT_EQUAL("TA", decompress(12 + 0, 2));
	ASSERT_EQUAL("TC", decompress(12 + 1, 2));
	ASSERT_EQUAL("TG", decompress(12 + 2, 2));
	ASSERT_EQUAL("TT", decompress(12 + 3, 2));
}

void TestMultyNucleotideDecompress() {
	ASSERT_EQUAL("CATGCA", decompress(124132, 6));
}

template<typename tType>
void checkArraysEqual(tType *a, tType *b, int length) {
	for (int i = 0; i < length; i++) {
		ASSERT_EQUAL(a[i], b[i]);
	}
}

void TestCodeReadAllA() {
	char *code = new char[readLength + 1];
	char *input = new char[readLength + 1];
	char *result = new char[readLength + 1];
	fill_n(input, readLength, 'A');
	fill_n(result, readLength, 0);
	codeRead(input, code);
	checkArraysEqual<char> (result, code, readLength);
	delete[] input;
	delete[] code;
}

void TestCodeReadAllACGT() {
	char *code = new char[readLength + 1];
	char *input = new char[readLength + 1];
	char *result = new char[readLength + 1];
	for (int i = 0; i < readLength; i++)
		if (i % 4 == 0) {
			input[i] = 'A';
			result[i] = 0;
		} else if (i % 4 == 1) {
			input[i] = 'C';
			result[i] = 1;
		} else if (i % 4 == 2) {
			input[i] = 'G';
			result[i] = 2;
		} else if (i % 4 == 3) {
			input[i] = 'T';
			result[i] = 3;
		}
	codeRead(input, code);
	checkArraysEqual<char> (result, code, readLength);
	delete[] input;
	delete[] code;
}

void TestExtractMerLen2() {
	char *input = new char[2];
	for (input[0] = 0; input[0] < 4; input[0]++)
		for (input[1] = 0; input[1] < 4; input[1]++)
			ASSERT_EQUAL(input[0] * 4 + input[1], extractMer(input, 0, 2));
	delete[] input;
}

void TestExtractMerLen3FromLen20() {
	int length = 20;
	char *input = new char[length];
	for (int i = 0; i < length; i++)
		input[i] = (100l << i) % 101 % 4;
	for (int i = 0; i + 2 < length; i++)
		ASSERT_EQUAL(input[i] * 16 + input[i + 1] * 4 + input[i + 2], extractMer(input, i, 3));
	delete[] input;
}

void TestExtractMerAllA() {
	char *input = new char[100];
	fill_n(input, readLength, 0);
	for (int merLength = 1; merLength < 100; merLength++) {
		for (int shift = 0; shift + merLength < 100; shift++) {
			ll result = extractMer(input, shift, merLength);
			ASSERT_EQUAL(0l, result);
		}
	}
	delete[] input;
}

void TestAddPairToTableSingleEntry() {
	myMap map;
	addPairToTable(map, 1, 2);
	ASSERT_EQUAL(1u, map.size());
	for (myMap::iterator it = map.begin(); it != map.end(); ++it) {
		ASSERT_EQUAL(1, it->first);
		ASSERT_EQUAL(1u, it->second.size());
		ASSERT_EQUAL(2, it->second[0]);
	}
}

void TestAddPairToTableTwoEntries() {
	myMap map;
	addPairToTable(map, 1, 2);
	addPairToTable(map, 3, 4);
	ASSERT_EQUAL(2u, map.size());
	ASSERT_EQUAL(1u, map[1].size());
	ASSERT_EQUAL(1u, map[3].size());
	ASSERT_EQUAL(2, map[1][0]);
	ASSERT_EQUAL(4, map[3][0]);
}

void TestAddPairToTableCoinsidingEntries() {
	myMap map;
	for(int i = 0; i < 10; i++)
		addPairToTable(map, 1, i);
	ASSERT_EQUAL(1u, map.size());
	for(int i = 0; i < 10; i++)
		ASSERT_EQUAL(i, map[1][i]);
}

cute::suite HashTableSuite() {
	cute::suite s;
	s.push_back(CUTE(TestEmptyDecompress));
	s.push_back(CUTE(TestSingleNucleotideDecompress));
	s.push_back(CUTE(TestDoubleNucleotideDecompress));
	s.push_back(CUTE(TestMultyNucleotideDecompress));
	s.push_back(CUTE(TestCodeReadAllA));
	s.push_back(CUTE(TestCodeReadAllACGT));
	s.push_back(CUTE(TestExtractMerLen2));
	s.push_back(CUTE(TestExtractMerLen3FromLen20));
	s.push_back(CUTE(TestExtractMerAllA));
	s.push_back(CUTE(TestAddPairToTableSingleEntry));
	s.push_back(CUTE(TestAddPairToTableTwoEntries));
	s.push_back(CUTE(TestAddPairToTableCoinsidingEntries));
	return s;
}

