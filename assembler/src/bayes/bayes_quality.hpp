/**
 * bayes_quality.h
 *
 *  Created on: Mar 23, 2011
 *      Author: snikolenko
 */
#ifndef BAYES_QUALITY_H_
#define BAYES_QUALITY_H_

#include <vector>
#include "sequence.hpp"
#include "logging.hpp"

using namespace std;

namespace bayes_quality {


// Vector of quality values
typedef vector<float> QVector;
// A read with quality values
typedef pair<Sequence, QVector> QRead;


// quality value of an insertion in the read
#define INSERT_Q 32

// quality value of a deletion in the read
#define DELETE_Q 32

// number of insertions allowed
#define INS 2

// number of deletions allowed
#define DEL 2

// magical constants
#define TENLOGELEVENTENTH 0.413926852
#define TENLOGTWO 3.010299957
#define VERYLARGEQ 999999


/**
 * @brief Vertex of a condensed graph.
 *
 * Class representing vertex of condensed graph.
 * Contains the Sequence of nucleotides, average coverage of this sequence,
 * pointers to children and coverage of corresponding edges.
 *
 */
class BayesQualityGenome {
private:
	Sequence genome_;
	vector<int> match_;
	string matchstr_, matchprettystr_, matchreadstr_;
	size_t inserts_, deletes_;
	size_t index_;
	int bestq_;
	int totalq_;
	
	vector<QVector*> qv;
	size_t qvsize_;
	
	int AddTwoQualityValues(int curq, int prevq);
	int AddThreeQualityValues(int diagq, int leftq, int rightq);
	float AddTwoQualityValues(float curq, float prevq);
	float AddThreeQualityValues(float diagq, float leftq, float rightq);


public:
	BayesQualityGenome(const char *genome) : genome_(genome), qv(INS+DEL+1) {
		for (size_t i=0; i < INS+DEL+1; ++i) {
			qv[i] = new QVector(genome_.size());
		}
	}
	
	~BayesQualityGenome() {
		for (size_t i=0; i < INS+DEL+1; ++i) {
			qv[i]->clear();  delete qv[i];
		}		
	}

	/**
	 * computes the likelihood of a single read
	 * @return total likelihood
	 */
	double ReadBQInt(QRead);
	
	/**
	 * returns the last set of matches
	 */
	const vector<int> & LastMatch() const { return match_; }
	
	/**
	 * returns the best index of the last match
	 */
	size_t LastMatchIndex() const { return index_; }

	/**
	 * returns the number of inserts in the best last match
	 */
	int LastMatchInserts() const { return inserts_; }

	/**
	 * returns the number of inserts in the best last match
	 */
	int LastMatchDeletes() const { return deletes_; }

	/**
	 * returns the total quality score in the the best last match
	 */
	int LastMatchQ() const { return bestq_; }

	/**
	 * returns the total quality score of the last read
	 */
	int LastTotalQ() const { return totalq_; }
	
	/**
	 * returns the last match genome string
	 */
	const string & LastMatchString() const { return matchstr_; }
	
	/**
	 * returns the last match read string
	 */
	const string & LastMatchReadString() const { return matchreadstr_; }

	/**
	 * returns the last match prettyfying string
	 */
	const string & LastMatchPrettyString() const { return matchprettystr_; }
};

}

#endif /* BAYES_QUALITY_H_ */

