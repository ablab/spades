#ifndef READGENERATOR_HPP_
#define READGENERATOR_HPP_
#include "iostream"
#include "fstream"
#include "string"
#include <cstdlib>
#include "strobe_read.hpp"
#include "vector"
#include "nucl.hpp"
#define MAX_PROBABILITY 100
using namespace std;

template<int size, int cnt = 1, typename T = char>
class ReadGenerator {
private:
	//	vector<ifaststream*> ifs_;
	string genome;
	int currentPosition;
	int maxPosition;
	int readNumber;
	int current;
	int insertLength;
	int errorProbability; // Probability in percents
	int errorDistribution[3];
	int insertError_;
	bool readingStarted_;
	int noErrorProbability(int p, int length) {
		return MAX_PROBABILITY - p * length + p * p * length * (length - 1) / 2
				/ MAX_PROBABILITY;
	}

public:
	ReadGenerator(const char *fileName, int coverage, int insert = 0) {
		insertLength = insert;
		ifstream s;
		s.open(fileName);
		s >> genome;
		maxPosition = genome.size() - size * cnt - insertLength * (cnt - 1);
		readNumber = (maxPosition + 1) * coverage;
		current = 0;
		setErrorProbability(0);
		setMaxInsertLengthError(0);
		readingStarted_ = false;
	}

	void close() {
	}

	void setRandSeed(unsigned int randSeed) {
		srand(randSeed);
	}

	void setErrorProbability(int probability) {
		if (readingStarted_) {
			cerr << "can not change generator parameters while reading" << endl;
			assert(1);
		}
		errorProbability = probability;
		errorDistribution[0] = noErrorProbability(probability, size);
		errorDistribution[1] = probability * size * noErrorProbability(
				probability, size) / MAX_PROBABILITY;
		errorDistribution[2] = MAX_PROBABILITY - errorDistribution[0]
				- errorDistribution[1];
	}

	void setMaxInsertLengthError(int insertError) {
		if (readingStarted_) {
			cerr << "can not change generator parameters while reading" << endl;
			assert(1);
		}
		insertError_ = insertError;
		currentPosition = insertError_;
	}

	ReadGenerator& operator>>(strobe_read<size, cnt, T> &sr) {
		if (eof()) {
			return *this;
		}
		if (!readingStarted_) {
			read_ahead();
			readingStarted_ = true;
		}
		sr = next_sr_;
		read_ahead();
		return *this;
	}

	inline bool is_open() const {
		return true;
	}

	inline bool eof() const {
		return current >= readNumber;
	}

	vector<strobe_read<size, cnt, T> >* readAll(int number = -1) {
		vector<strobe_read<size, cnt, T> > *v = new vector<strobe_read<size,
				cnt, T> > ();
		strobe_read<size, cnt, T> sr;
		while (!eof() && number--) {
			this->operator>>(sr);
			v->push_back(sr);
		}
		return v;
	}

private:
	strobe_read<size, cnt, T> next_sr_;

	inline void read_ahead() {
		while (!eof() && !read(next_sr_)) {
			;
		}
	}

	inline int numberOfErrors() {
		int dice = rand() % MAX_PROBABILITY;
		if (dice < errorDistribution[0])
			return 0;
		if (dice < errorDistribution[1] + errorDistribution[0])
			return 1;
		return 2;
	}

	void introduceErrors(string &tmp) {
		int errors = numberOfErrors();
		for (int j = 0; j < errors; j++) {
			int pos = rand() % size;
			tmp[pos] = nucl(dignucl(tmp[pos]) ^ (1 + (rand() % 3)));
		}
	}

	bool read(strobe_read<size, cnt, T> &sr) {
		if (!is_open() || eof()) {
			return false;
		}
		int p = currentPosition;
		for (int i = 0; i < cnt; i++) {
			string readString = genome.substr(
					p + rand() % (2 * insertError_ + 1) - insertError_, size);
			introduceErrors(readString);
			sr.put(i, readString);
			p += size + insertLength;
		}
		current++;
		currentPosition++;
		if (currentPosition > maxPosition - insertError_)
			currentPosition = insertError_;
		return true;
	}
};

#endif /* READGENERATOR_HPP_ */
