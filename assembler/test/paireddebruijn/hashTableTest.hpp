#include "cute.h"
#include "seq.hpp"
#include "nucl.hpp"
#include "constructHashTable.hpp"

void TestEmptyDecompress() {
	ASSERT_EQUAL(decompress(0, 0), "");
}

void TestSingleNucleotideDecompress() {
	ASSERT_EQUAL(decompress(0, 1), "A");
	ASSERT_EQUAL(decompress(1, 1), "C");
	ASSERT_EQUAL(decompress(2, 1), "G");
	ASSERT_EQUAL(decompress(3, 1), "T");
}

void TestDoubleNucleotideDecompress() {
	ASSERT_EQUAL(decompress(0, 2), "AA");
	ASSERT_EQUAL(decompress(1, 2), "AC");
	ASSERT_EQUAL(decompress(2, 2), "AG");
	ASSERT_EQUAL(decompress(3, 2), "AT");
	ASSERT_EQUAL(decompress(4 + 0, 2), "CA");
	ASSERT_EQUAL(decompress(4 + 1, 2), "CC");
	ASSERT_EQUAL(decompress(4 + 2, 2), "CG");
	ASSERT_EQUAL(decompress(4 + 3, 2), "CT");
	ASSERT_EQUAL(decompress(8 + 0, 2), "GA");
	ASSERT_EQUAL(decompress(8 + 1, 2), "GC");
	ASSERT_EQUAL(decompress(8 + 2, 2), "GG");
	ASSERT_EQUAL(decompress(8 + 3, 2), "GT");
	ASSERT_EQUAL(decompress(12 + 0, 2), "TA");
	ASSERT_EQUAL(decompress(12 + 1, 2), "TC");
	ASSERT_EQUAL(decompress(12 + 2, 2), "TG");
	ASSERT_EQUAL(decompress(12 + 3, 2), "TT");
}

void TestMultyNucleotideDecompress() {
	ASSERT_EQUAL(decompress(124132, 6), "CATGCA");
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
	checkArraysEqual<char>(result, code, readLength);
	delete input;
	delete code;
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
	for (int i = 0; i < readLength; i++)
		ASSERT_EQUAL(result[0], code[0]);
	checkArraysEqual<char>(result, code, readLength);
	delete input;
	delete code;
}



cute::suite HashTableSuite() {
	cute::suite s;
	s.push_back(CUTE(TestEmptyDecompress));
	s.push_back(CUTE(TestSingleNucleotideDecompress));
	s.push_back(CUTE(TestDoubleNucleotideDecompress));
	s.push_back(CUTE(TestMultyNucleotideDecompress));
	s.push_back(CUTE(TestCodeReadAllA));
	s.push_back(CUTE(TestCodeReadAllACGT));
	return s;
}

