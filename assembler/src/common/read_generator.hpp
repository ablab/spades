#ifndef READGENERATOR_HPP_
#define READGENERATOR_HPP_
#include "iostream"
#include "fstream"
#include "string"
#include "strobe_read.hpp"
#include <cstdlib>
#include "vector"
#include "nucl.hpp"

#define MAX_PROBABILITY 10000
using namespace std;

class SmoothPositionChooser {
public:
	static int choosePosition(int currentReadNumber, int readNumber,
			int minPosition, int maxPosition) {
		return minPosition + (maxPosition - minPosition)
				* (long long) currentReadNumber / (readNumber - 1);
	}
};

class RandomPositionChooser {
public:
	static int choosePosition(int currentReadNumber, int readNumber,
			int minPosition, int maxPosition) {
		return minPosition + rand() % (maxPosition - minPosition + 1);
	}
};

template<int size, int cnt = 1, typename T = int,
		typename PositionChooser = SmoothPositionChooser>
class ReadGenerator {
private:
	//	vector<ifaststream*> ifs_;
	string genome_;
	int minPosition_;
	int maxPosition_;
	int coverage_;
	int readNumber_;
	int currentReadNumber_;
	int insertLength_;
	int errorProbability_; // Probability in percents
	int errorDistribution_[3];
	int insertError_;
	bool readingStarted_;
	double noErrorProbability(int p, int length) {
		return MAX_PROBABILITY - p * length + p * p * length * (length - 1) / 2
				/ MAX_PROBABILITY;
	}

public:
	void initParameters(int coverage, int insert) {
		insertLength_ = insert;
		coverage_ = coverage;
		readNumber_ = coverage_ * genome_.size() / (size * cnt);
		currentReadNumber_ = 0;
		readingStarted_ = false;
		setErrorProbability(0);
		setMaxInsertLengthError(0);
	}

	ReadGenerator(istream is, int coverage, int insert = 0) {
		is >> genome_;
		initParameters(coverage, insert);
	}

	ReadGenerator(string genome, int coverage, int insert = 0) {
		genome_ = genome;
		initParameters(coverage, insert);
	}

	void close() {
	}

	void setRandSeed(unsigned int randSeed) {
		srand(randSeed);
	}

	void setErrorProbability(double probability) {
		setErrorProbability((int) (probability * MAX_PROBABILITY));
	}

	void setErrorProbability(int probability) {
		if (readingStarted_) {
			cerr << "can not change generator parameters while reading" << endl;
			assert(1);
		}
		errorProbability_ = probability;
		errorDistribution_[0] = noErrorProbability(probability, size);
		errorDistribution_[1] = probability * size * noErrorProbability(
				probability, size) / MAX_PROBABILITY;
		errorDistribution_[2] = MAX_PROBABILITY - errorDistribution_[0]
				- errorDistribution_[1];
	}

	void setMaxInsertLengthError(int insertError) {
		if (readingStarted_) {
			cerr << "can not change generator parameters while reading" << endl;
			assert(1);
		}
		insertError_ = insertError / 2;
		maxPosition_ = genome_.size() - size * cnt - insertLength_ * (cnt - 1)
				- insertError_;
		minPosition_ = insertError_;
	}

	ReadGenerator& operator>>(strobe_read<size, cnt, T> &sr) {
		if (eof()) {
			return *this;
		}
		if (!readingStarted_) {
			readingStarted_ = true;
			read_ahead();
		}
		sr = next_sr_;
		read_ahead();
		return *this;
	}

	//todo think about interface
	ReadGenerator& operator>>(Read &r) {
		if (eof()) {
			return *this;
		}
		if (!readingStarted_) {
			readingStarted_ = true;
			read_ahead();
		}
		r = Read("", next_sr_[0].str(), "");
		read_ahead();
		return *this;
	}

	inline bool is_open() const {
		return true;
	}

	inline bool eof() const {
		return currentReadNumber_ >= readNumber_;
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
		if (dice < errorDistribution_[0])
			return 0;
		if (dice < errorDistribution_[1] + errorDistribution_[0])
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
		//				cout << currentReadNumber_ << " " << readNumber_ << endl;
		if (!is_open() || eof()) {
			return false;
		}
		int p = PositionChooser::choosePosition(currentReadNumber_,
				readNumber_, minPosition_, maxPosition_);
		cout << p << endl;
		for (int i = 0; i < cnt; i++) {
			int positionError = rand() % (2 * insertError_ + 1) - insertError_;
			string readString = genome_.substr(p + positionError, size);
			introduceErrors(readString);
			sr.put(i, readString);
			p += size + insertLength_;
		}
		currentReadNumber_++;
		return true;
	}
};

//todo reduce copy-paste (graphio.hpp)
Sequence readGenome(istream &is) {
	SequenceBuilder sb;
	string buffer;
	while(!is.eof()){
		is >> buffer;
		sb.append(Sequence(buffer));
	}
	return sb.BuildSequence();
}

Sequence readGenomeFromFile(const string &fileName) {
	ifstream is;
	is.open(fileName.c_str());
	Sequence result(readGenome(is));
	is.close();
	return result;
}
template<typename PositionChooser>
void generateReads(string fileName, string genomeFileName, int insertLength,
		int coverage, double errorProbability, int maxInsertLengthError) {
	ofstream os;
	os.open(fileName.c_str());
	Sequence genome(readGenomeFromFile(genomeFileName));
	stringstream ss;
	ss << genome;
	ReadGenerator<100, 2, int, PositionChooser> gen(ss.str(), coverage,
			insertLength);
	gen.setErrorProbability(errorProbability);
	gen.setMaxInsertLengthError(maxInsertLengthError);
	strobe_read<100, 2> readPair;
	while (!gen.eof()) {
		gen >> readPair;
		os << readPair[0] << " " << readPair[1] << endl;
	}
	os.close();
}

#endif /* READGENERATOR_HPP_ */
