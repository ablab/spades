#ifndef READGENERATOR_HPP_
#define READGENERATOR_HPP_
#include "iostream"
#include "fstream"
#include "string"
#include "strobe_read.hpp"
#include <cstdlib>
#include "vector"
#include "nucl.hpp"
#include "read.hpp"

/**
 * @file read_generator.hpp
 */

#define MAX_PROBABILITY 10000
using namespace std;

class SmoothPositionChooser {
public:
	static int choosePosition(int currentReadNumber, int readNumber,
			int minPosition, int maxPosition) {
		if (readNumber == 1) {
			return (minPosition + maxPosition) / 2;
		} else
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

template<typename PositionChooser = SmoothPositionChooser>
class ReadGenerator {
private:
	//	vector<ifaststream*> ifs_;
	string genome_;
	size_t size_;
	size_t cnt_;
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
	bool closed_;

	double ErrorFreeProbability(int p, int length) {
		return MAX_PROBABILITY - p * length + p * p * length * (length - 1) / 2
				/ MAX_PROBABILITY;
	}

	void inner_set_error_probability(int probability) {
		if (readingStarted_) {
			cerr << "can not change generator parameters while reading" << endl;
			assert(1);
		}
		errorProbability_ = probability;
		errorDistribution_[0] = ErrorFreeProbability(probability, size_);
		errorDistribution_[1] = probability * size_ * ErrorFreeProbability(
				probability, size_) / MAX_PROBABILITY;
		errorDistribution_[2] = MAX_PROBABILITY - errorDistribution_[0]
				- errorDistribution_[1];
	}

	vector<Read> next_sr_;

	void ReadAhead() {
		if (is_open() && !eof()) {
			int p = PositionChooser::choosePosition(currentReadNumber_,
					readNumber_, minPosition_, maxPosition_);
			for (int i = 0; i < cnt_; i++) {
				int positionError = rand() % (2 * insertError_ + 1)
						- insertError_;
				string read = genome_.substr(p + positionError, size_);
				IntroduceErrors(read);
				next_sr_[i] = Read("", read, "");
				p += size_ + insertLength_;
			}
			currentReadNumber_++;
		}
	}

	int NumberOfErrors() {
		int dice = rand() % MAX_PROBABILITY;
		if (dice < errorDistribution_[0])
			return 0;
		if (dice < errorDistribution_[1] + errorDistribution_[0])
			return 1;
		return 2;
	}

	void IntroduceErrors(string &tmp) {
		int errors = NumberOfErrors();
		for (int j = 0; j < errors; j++) {
			int pos = rand() % size_;
			tmp[pos] = nucl(dignucl(tmp[pos]) ^ (1 + (rand() % 3)));
		}
	}

	void InitParameters(int coverage, int insert) {
		insertLength_ = insert;
		coverage_ = coverage;
		readNumber_ = coverage_ * genome_.size() / (size_ * cnt_);
		currentReadNumber_ = 0;
		readingStarted_ = false;
		set_error_probability(0);
		set_max_insert_length_error(0);
	}

public:

	ReadGenerator(size_t count, size_t size, istream is, int coverage,
			int insert = 0) :
		cnt_(count), size_(size), closed_(false) {
		is >> genome_;
		InitParameters(coverage, insert);
	}

	ReadGenerator(size_t count, size_t size, string genome, int coverage,
			int insert = 0) :
		cnt_(count), size_(size), closed_(false) {
		genome_ = genome;
		InitParameters(coverage, insert);
	}

	void close() {
		closed_ = true;
	}

	void set_rand_seed(unsigned int randSeed) {
		srand(randSeed);
	}

	void set_error_probability(double probability) {
		inner_set_error_probability((int) (probability * MAX_PROBABILITY));
	}

	void set_max_insert_length_error(int insertError) {
		if (readingStarted_) {
			cerr << "can not change generator parameters while reading" << endl;
			assert(1);
		}
		insertError_ = insertError / 2;
		maxPosition_ = genome_.size() - size_ * cnt_ - insertLength_ * (cnt_ - 1)
				- insertError_;
		minPosition_ = insertError_;
	}

	ReadGenerator& operator>>(vector<Read>& reads) {
		if (eof()) {
			return *this;
		}
		if (!readingStarted_) {
			readingStarted_ = true;
			ReadAhead();
		}
		reads = next_sr_;
		ReadAhead();
		return *this;
	}

	ReadGenerator& operator>>(PairedRead &p_r) {
		assert(cnt_ ==2);
		if (eof()) {
			return *this;
		}
		if (!readingStarted_) {
			readingStarted_ = true;
			ReadAhead();
		}
		p_r = PairedRead(next_sr_[0], next_sr_[1], insertLength_ + 2 * size_);
		ReadAhead();
		return *this;
	}

	//todo think about interface
	ReadGenerator& operator>>(Read &r) {
		assert(cnt_ == 1);

		if (eof()) {
			return *this;
		}
		if (!readingStarted_) {
			readingStarted_ = true;
			ReadAhead();
		}
		r = next_sr_[0];
		ReadAhead();
		return *this;
	}

	bool is_open() const {
		return !closed_;
	}

	bool eof() const {
		return currentReadNumber_ >= readNumber_;
	}

//	vector<strobe_read<size, cnt, T> >* readAll(int number = -1) {
//		vector<strobe_read<size, cnt, T> > *v = new vector<strobe_read<size,
//				cnt, T> > ();
//		strobe_read<size, cnt, T> sr;
//		while (!eof() && number--) {
//			this->operator>>(sr);
//			v->push_back(sr);
//		}
//		return v;
//	}

};

//todo reduce copy-paste (graphio.hpp)
Sequence ReadGenome(istream &is) {
	SequenceBuilder sb;
	string buffer;
	while (!is.eof()) {
		is >> buffer;
		sb.append(Sequence(buffer));
	}
	return sb.BuildSequence();
}

Sequence ReadGenomeFromFile(const string &fileName) {
	ifstream is;
	is.open(fileName.c_str());
	Sequence result(ReadGenome(is));
	is.close();
	return result;
}

//template<typename PositionChooser>
//void GenerateReads(string fileName, string genomeFileName, int insertLength,
//		int coverage, double errorProbability, int maxInsertLengthError) {
//	ofstream os;
//	os.open(fileName.c_str());
//	Sequence genome(ReadGenomeFromFile(genomeFileName));
//	stringstream ss;
//	ss << genome;
//	ReadGenerator<PositionChooser> gen(ss.str(), coverage,
//			insertLength);
//	gen.setErrorProbability(errorProbability);
//	gen.set_max_insert_length_error(maxInsertLengthError);
//	strobe_read<100, 2> readPair;
//	while (!gen.eof()) {
//		gen >> readPair;
//		os << readPair[0] << " " << readPair[1] << endl;
//	}
//	os.close();
//}

#endif /* READGENERATOR_HPP_ */
