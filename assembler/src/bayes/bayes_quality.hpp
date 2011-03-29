/**
 * bayes_quality.h
 *
 *  Created on: Mar 23, 2011
 *      Author: snikolenko
 */
#ifndef BAYES_QUALITY_H_
#define BAYES_QUALITY_H_

#include <vector>
#include <map>
#include <tr1/unordered_map>
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
#define INS 1

// number of deletions allowed
#define DEL 1

// magical constants
#define TENLOGELEVENTENTH 0.413926852
#define TENLOGTWO 3.010299957
#define VERYLARGEQ 999999

// if the total quality value of the best match is this good,
// we will skip the complementary read and assume this is the one
#define GOOD_ENOUGH_Q 100

// in preprocessing, we construct a map of sequences of this size
// to check which reads are worth pursuing at a given point
#define PREPROCESS_SEQ_LENGTH 10
#define NEED_INTERSECTION 6
typedef Seq<PREPROCESS_SEQ_LENGTH> PSeq;
typedef tr1::unordered_map<PSeq, vector<bool>, PSeq::hash, PSeq::equal_to> PreprocMap;
typedef tr1::unordered_map<PSeq, pair<size_t, size_t>, PSeq::hash, PSeq::equal_to> PreprocReadsMap;

/// @typedef a structure for storing match results in BayesQualityGenome
typedef struct {
	string matchstr_, matchprettystr_, matchreadstr_;
	size_t inserts_, deletes_;
	size_t index_;
	int bestq_;
	int totalq_;
	vector<int> match_;
} MatchResults;

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
	size_t gensize_;
	MatchResults mr_;
	
	// for map-like preprocessing
	PreprocMap availableReadsHead_;
	PreprocMap availableReadsMid_;
	vector<size_t> readMids;
	void fillAvailableReads(size_t vecSize, const string & seq);
	void PreprocessOneMapOneReadWithShift(const PSeq &ps, PreprocReadsMap & rm, int mapno, PreprocMap & map, size_t readno, size_t reads_size, size_t shift);

	vector<QVector*> qv;
	size_t qvsize_;

	// these functions add up two/three logarithms; currently they do it in a very approximate fashion (by taking the min)
	int AddTwoQualityValues(int curq, int prevq);
	int AddThreeQualityValues(int diagq, int leftq, int rightq);
	float AddTwoQualityValues(float curq, float prevq);
	float AddThreeQualityValues(float diagq, float leftq, float rightq);

	/**
	 * process one read, store all results in the corresponding fields, return total likelihood
	 * @param read_no number of this read in the preprocessed maps; if there was no preprocessing, put -1 here
	 */
	double ProcessOneReadBQ(const Sequence &, const QVector &, size_t readno = -1);
	
	// auxiliary function: check whether this read can be applied to this genome part
	bool isAvailable(size_t readno, size_t j, const PSeq & curpseq);

public:
	BayesQualityGenome(const char *genome) : genome_(genome), qv(INS+DEL+1) {
		for (size_t i=0; i < INS+DEL+1; ++i) {
			qv[i] = new QVector(genome_.size());
		}
		gensize_ = genome_.size();
	}
	
	~BayesQualityGenome() {
		for (size_t i=0; i < INS+DEL+1; ++i) {
			qv[i]->clear();  delete qv[i];
		}
	}

	/**
	 * computes the likelihood of a single read
	 * returns the minimum likelihood of a read and its complement
	 * @return total likelihood
	 */
	double ReadBQ(const QRead &);

	/**
	 * computes the likelihood of a single read after preprocessing it in a set of reads
	 * returns the minimum likelihood of a read and its complement
	 * @return total likelihood
	 */
	double ReadBQPreprocessed(const QRead &, size_t readno, size_t readssize);

	/**
	 * preprocess a vector of reads, i.e., be prepared to compute their masks
	 */
	void PreprocessReads(const vector<QRead *> &);
	
	/**
	 * process a vector of reads and log the results
	 */	
	void ProcessReads(const vector<QRead *> &);
	
	/**
	 * returns the last set of matches
	 */
	const vector<int> & LastMatch() const { return mr_.match_; }
	
	/**
	 * returns the best index of the last match
	 */
	size_t LastMatchIndex() const { return mr_.index_; }

	/**
	 * returns the number of inserts in the best last match
	 */
	int LastMatchInserts() const { return mr_.inserts_; }

	/**
	 * returns the number of inserts in the best last match
	 */
	int LastMatchDeletes() const { return mr_.deletes_; }

	/**
	 * returns the total quality score in the the best last match
	 */
	int LastMatchQ() const { return mr_.bestq_; }

	/**
	 * returns the total quality score of the last read
	 */
	int LastTotalQ() const { return mr_.totalq_; }
	
	/**
	 * returns the last match genome string
	 */
	const string & LastMatchString() const { return mr_.matchstr_; }
	
	/**
	 * returns the last match read string
	 */
	const string & LastMatchReadString() const { return mr_.matchreadstr_; }

	/**
	 * returns the last match prettyfying string
	 */
	const string & LastMatchPrettyString() const { return mr_.matchprettystr_; }
};

}

#endif /* BAYES_QUALITY_H_ */

