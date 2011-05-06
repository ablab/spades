/**
 * @file read_generator.hpp
 */
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

namespace read_generator {
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
	size_t cnt_;
	size_t size_;
	int min_position_;
	int max_position_;
	int coverage_;
	int read_count_;
	int current_read_number_;
	int insertLength_;
	int error_probability_; // Probability in percents
	int error_distribution_[3];
	int insert_error_;
	bool reading_started_;
	bool closed_;

	double ErrorFreeProbability(int p, int length) {
		return MAX_PROBABILITY - p * length + p * p * length * (length - 1) / 2
				/ MAX_PROBABILITY;
	}

	void inner_set_error_probability(int probability) {
		if (reading_started_) {
			cerr << "can not change generator parameters while reading" << endl;
			assert(1);
		}
		error_probability_ = probability;
		error_distribution_[0] = ErrorFreeProbability(probability, size_);
		error_distribution_[1] = probability * size_ * ErrorFreeProbability(
				probability, size_) / MAX_PROBABILITY;
		error_distribution_[2] = MAX_PROBABILITY - error_distribution_[0]
				- error_distribution_[1];
	}

	vector<Read> next_sr_;

	void ReadAhead() {
		if (is_open() && !eof()) {
			int p = PositionChooser::choosePosition(current_read_number_,
					read_count_, min_position_, max_position_);
			for (size_t i = 0; i < cnt_; i++) {
				int position_error = rand() % (2 * insert_error_ + 1)
						- insert_error_;
				string read = genome_.substr(p + position_error, size_);
				IntroduceErrors(read);
				next_sr_[i] = Read("", read, "");
				p += size_ + insertLength_;
			}
			current_read_number_++;
		}
	}

	int NumberOfErrors() {
		int dice = rand() % MAX_PROBABILITY;
		if (dice < error_distribution_[0])
			return 0;
		if (dice < error_distribution_[1] + error_distribution_[0])
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
		read_count_ = coverage_ * genome_.size() / (size_ * cnt_);
		current_read_number_ = 0;
		reading_started_ = false;
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

	void reset() {
		reading_started_ = false;
		current_read_number_ = 0;
	}

	void close() {
		closed_ = true;
	}

	void set_rand_seed(unsigned int rand_seed) {
		srand(rand_seed);
	}

	void set_error_probability(double probability) {
		inner_set_error_probability((int) (probability * MAX_PROBABILITY));
	}

	void set_max_insert_length_error(int insert_error) {
		if (reading_started_) {
			cerr << "can not change generator parameters while reading" << endl;
			assert(1);
		}
		insert_error_ = insert_error / 2;
		max_position_ = genome_.size() - size_ * cnt_ - insertLength_ * (cnt_ - 1)
				- insert_error_;
		min_position_ = insert_error_;
	}

	ReadGenerator& operator>>(vector<Read>& reads) {
		if (eof()) {
			return *this;
		}
		if (!reading_started_) {
			reading_started_ = true;
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
		if (!reading_started_) {
			reading_started_ = true;
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
		if (!reading_started_) {
			reading_started_ = true;
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
		return current_read_number_ >= read_count_;
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
}
#endif /* READGENERATOR_HPP_ */
