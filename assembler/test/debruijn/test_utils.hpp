/*
 * test_utils.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: sergey
 */

#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

#include "read.hpp"
#include "ireadstream.hpp"
#include "read_generator.hpp"

#define SUBSTR_LENGTH 1000
#define COVERAGE 1
#define R 35
//#define K 15
#define filename "./data/input/MG1655-K12.fasta.gz"
#define readsfilename "./data/input/s_6.first100000_1.fastq.gz"

vector<Read> GenerateReadsWithMistakes() {
	LOGGER("d.test_utils");
	INFO("Reading " << filename);

	ireadstream stream(filename);
	Read r;
	stream >> r;

	vector<Read> reads;
	INFO("Closing " << filename);
	stream.close();
	INFO("Generating reads for substring of length " << SUBSTR_LENGTH << " and coverage " << COVERAGE);

	ReadGenerator<R> gen(r.getSequenceString().substr(0, SUBSTR_LENGTH), COVERAGE);
//	gen.setErrorProbability(2);

	while (!gen.eof()) {
		Read read;
		gen >> read;
//		cout << read[0] << endl;
		reads.push_back(read);
	}

	INFO("Reads generated");
	return reads;
}

vector<Read> ReadFromFile() {
	LOGGER("d.test_utils");
	INFO("Reading " << readsfilename);

	vector<Read> reads;
	ireadstream stream(readsfilename);

	while (!stream.eof()) {
		Read r;
		stream >> r;
		reads.push_back(r);
	}

	INFO("Closing " << readsfilename);
	stream.close();

	return reads;

}

template<typename T>
struct PairHash {
	size_t operator()(pair<T, T> p) const {
		return hash<T> ()(p.first) + hash<T> ()(p.second);
	}
};

template<typename T>
struct PairLess {
	bool operator()(pair<T, T> p1, pair<T, T> p2) const {
		return less<T> ()(p1.first, p2.first) ? true : (less<T> ()(p2.first,
				p1.first) ? false : less<T> ()(p1.second, p2.second));
	}
};

std::string complement(const std::string& s) {
	return (!Sequence(s)).str();
}

vector<Read> MakeReads(string *ss, size_t count) {
	vector<Read> ans;
	for (size_t i = 0; i < count; ++i) {
		Read r("", *ss, "");
		ss++;
		ans.push_back(r);
	}
	return ans;
}

#endif /* TEST_UTILS_HPP_ */
