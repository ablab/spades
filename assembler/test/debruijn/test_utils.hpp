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
//#include "read_generator.hpp"

#define ECOLI_FILE "./data/input/MG1655-K12.fasta.gz"
#define QUAKE_CROPPED_10_3 make_pair(make_pair("./data/input/s_6.first1000_1.fastq.gz", "./data/input/s_6.first1000_2.fastq.gz"), 10000)
#define QUAKE_CROPPED_10_4 make_pair(make_pair("./data/input/s_6.first10000_1.fastq.gz","./data/input/s_6.first10000_2.fastq.gz"), 10000)
#define QUAKE_CROPPED_10_5 make_pair(make_pair("./data/input/s_6.first100000_1.fastq.gz","./data/input/s_6.first100000_2.fastq.gz"), 100000)
#define QUAKE_CROPPED_4_10_5 make_pair(make_pair("./data/input/s_6.first400000_1.fastq.gz","./data/input/s_6.first400000_2.fastq.gz"), 400000)

//#define filename "./data/input/MG1655-K12.fasta.gz"
#define readsfilename "./data/input/s_6.first100000_1.fastq.gz"

//vector<Read> GenerateReadsWithMistakes(const string& file_name) {
//	LOGGER("d.test_utils");
//	INFO("Reading " << file_name);
//
//	ireadstream stream(file_name);
//	Read r;
//	stream >> r;
//
//	vector<Read> reads;
//	INFO("Closing " << file_name);
//	stream.close();
//	INFO("Generating reads for substring of length " << SUBSTR_LENGTH << " and coverage " << COVERAGE);
//
//	ReadGenerator<R> gen(r.getSequenceString().substr(0, SUBSTR_LENGTH), COVERAGE);
////	gen.setErrorProbability(2);
//
//	while (!gen.eof()) {
//		Read read;
//		gen >> read;
////		cout << read[0] << endl;
//		reads.push_back(read);
//	}
//
//	INFO("Reads generated");
//	return reads;
//}

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
