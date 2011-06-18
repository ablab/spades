/**
 * bayes_quality.h
 *
 *  Created on: Mar 23, 2011
 *      Author: snikolenko
 *
 */
#ifndef BAYES_QUALITY_H_
#define BAYES_QUALITY_H_

#include <vector>
#include <map>
#include <tr1/unordered_map>
#include "sequence.hpp"
#include "logging.hpp"
#include "read.hpp"

using namespace std;

namespace bayes_quality {


// Vector of quality values
typedef vector<float> QVector;
// A read with quality values
typedef pair<Sequence, QVector> QRead;
// Bowtie search results
typedef vector< vector< vector<size_t> > > BowtieResults;

// parameters for bowtie interaction
#define USE_BOWTIE
#define BOWTIE_SEED_LENGTH 15
#define BOWTIE_SEED_NUM 5
#define BWMAP_IND_REVC 0
#define BWMAP_IND_RDNO 1
#define BWMAP_IND_RIND 2
#define BWMAP_IND_GIND 3


// size of read batch that we read from file
#define READ_BATCH 10000

// number of threads
#define THREADS_NUM 3

// number of reads to skip from the beginning of the file
#define SKIP_READS 0

// write a file with statistics
#define WRITE_STATSFILE
#define STATSFILENAME "readstats.txt"

// minimum size for a read
// if a read is shorter (has more N's), it will be skipped
#define MIN_READ_SIZE 60

#define BIGREADNO 100000

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
#define GOOD_ENOUGH_Q 50

// in preprocessing, we construct a map of sequences of this size
// to check which reads are worth pursuing at a given point
// #define USE_PREPROCESSING
#define PREPROCESS_SEQ_LENGTH 10
#define NEED_INTERSECTION 6
typedef Seq<PREPROCESS_SEQ_LENGTH> PSeq;
typedef tr1::unordered_map<PSeq, vector<bool>, PSeq::hash, PSeq::equal_to> PreprocMap;
typedef tr1::unordered_map<PSeq, pair<size_t, size_t>, PSeq::hash, PSeq::equal_to> PreprocReadsMap;


// in look-ahead, we test genome positions in the following way:
//   -- if nucleotides [0, LA_SIZE-1] of the read match the current genome position j with more than LA_ERRORS errors, it's good
//   -- if nucleotides [mid, mid+LA_SIZE-1] of the read match one of positions j+mid-DEL, ..., j+mid+INS-1 with more than LA_ERRORS errors, it's good
//   -- otherwise, it's bad, skip it
// this works well if INS and DEL are small (e.g., INS=DEL=1)
// if INS>1 or DEL>1, this approach may lead to errors
// #define USE_LOOKAHEAD
#define LA_SIZE 25
#define LA_ERRORS 7
typedef Seq<LA_SIZE> LASeq;


/// @typedef a structure for storing match results in BayesQualityGenome
typedef struct {
	string matchstr_, matchprettystr_, matchreadstr_;
	size_t inserts_, deletes_;
	size_t index_;
	int bestq_;
	int totalq_;
	vector<int> match_;
	double prob_;
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

	// do we use external bowtie
	bool bowtie_;
	string bowtieCmd_;
	string bowtieIndex_;
	BowtieResults bwres_;

	// for run statistics
	size_t totalPos_, totalGood_;
	
	// for map-like preprocessing
	PreprocMap availableReadsHead_;
	PreprocMap availableReadsMid_;
	vector<size_t> readMids;
	void fillAvailableReads(size_t vecSize, const string & seq);
	void PreprocessOneMapOneReadWithShift(const PSeq &ps, PreprocReadsMap & rm, int mapno, PreprocMap & map, size_t readno, size_t reads_size, size_t shift);

	//vector<QVector*> qv;
	//size_t qvsize_;

	// these functions add up two/three logarithms; currently they do it in a very approximate fashion (by taking the min)
	int AddTwoQualityValues(int curq, int prevq);
	int AddThreeQualityValues(int diagq, int leftq, int rightq);
	float AddTwoQualityValues(float curq, float prevq);
	float AddThreeQualityValues(float diagq, float leftq, float rightq);

	// simple match for the lookahead
	int simpleMatch(const LASeq & seq, const LASeq & gen);

	/**
	 * process one read, store all results in the corresponding fields, return total likelihood
	 * @param read_no number of this read in the preprocessed maps; if there was no preprocessing, put -1 here
	 * @param indices if we already have a list of indices to check, insert them here; if the input has size 0, we do full search
	 */
	MatchResults ProcessOneReadBQ(const Sequence &, const QVector &, size_t readno, const vector<size_t> & indices);
	
	// auxiliary function: check whether this read can be applied to this genome part
	bool isAvailable(size_t readno, size_t j, const PSeq & curpseq);

public:
	/*BayesQualityGenome(const char *genome) : genome_(genome), totalPos_(0), totalGood_(0), qv(INS+DEL+1) {
		for (size_t i=0; i < INS+DEL+1; ++i) {
			qv[i] = new QVector(genome_.size());
		}
		gensize_ = genome_.size();
	}*/

	BayesQualityGenome(const char *genome) : genome_(genome), totalPos_(0), totalGood_(0) {
		gensize_ = genome_.size();
	}
	
	~BayesQualityGenome() {
		/*for (size_t i=0; i < INS+DEL+1; ++i) {
			qv[i]->clear();  delete qv[i];
		}*/
	}

	/**
	 * computes the likelihood of a single read
	 * returns the minimum likelihood of a read and its complement
	 * @return total likelihood
	 */
	double ReadBQ(const Read &);

	/**
	 * computes the likelihood of a single read after preprocessing it in a set of reads
	 * returns the minimum likelihood of a read and its complement
	 * @return complete match results
	 */
	MatchResults ReadBQPreprocessed(const Read &, size_t readno, size_t readssize, const BowtieResults & bwres);

	/**
	 * preprocess a vector of reads, i.e., be prepared to compute their masks
	 */
	void PreprocessReads(const vector<Read> &);
	
	/**
	 * process a vector of reads and log the results
	 */	
	void ProcessReads(const vector<Read> &);
	
	/**
	 * process all reads from a file and log the results
	 */	
	void ProcessReads(const char * filename);


	/**
	 * returns the last set of matches
	 */
//	const vector<int> & LastMatch() const { return mr_.match_; }
	
	/**
	 * returns the best index of the last match
	 */
//	size_t LastMatchIndex() const { return mr_.index_; }

	/**
	 * returns the number of inserts in the best last match
	 */
//	int LastMatchInserts() const { return mr_.inserts_; }

	/**
	 * returns the number of inserts in the best last match
	 */
//	int LastMatchDeletes() const { return mr_.deletes_; }

	/**
	 * returns the total quality score in the the best last match
	 */
//	int LastMatchQ() const { return mr_.bestq_; }

	/**
	 * returns the total quality score of the last read
	 */
//	int LastTotalQ() const { return mr_.totalq_; }
	
	/**
	 * returns the last match genome string
	 */
//	const string & LastMatchString() const { return mr_.matchstr_; }
	
	/**
	 * returns the last match read string
	 */
//	const string & LastMatchReadString() const { return mr_.matchreadstr_; }

	/**
	 * returns the last match prettyfying string
	 */
//	const string & LastMatchPrettyString() const { return mr_.matchprettystr_; }



	/**
	 * returns the total number of possible matching positions processed
	 */
	size_t TotalPositions() const { return totalPos_; }
	/**
	 * returns the total number of ``good'' positions
	 */
	size_t TotalGoodPositions() const { return totalGood_; }


	/// set bowtie parameters
	void setBowtie(string cmd, string ind) {
		bowtie_ = true;
		bowtieCmd_ = cmd;
		bowtieIndex_ = ind;
	}
	
	
	// Bowtie-related routines
	
	/// write down a batch of reads for Bowtie processing
	void writeReadPartsForBowtie(const vector<Read> & v, int fd, unsigned int readno);
	/// write down a batch of reads for Bowtie processing to a temporary file
	/// return its name
	string writeReadPartsForBowtie(const vector<Read> & v, unsigned int readno) {
		char *tmpname = strdup("/tmp/tmpfileXXXXXX");
		int fd = mkstemp(tmpname);
		assert(fd >= 0);
		writeReadPartsForBowtie(v, fd, readno);
		return tmpname;
	}
	/// read Bowtie search results
	BowtieResults readBowtieResults(int fd, size_t noofreads, unsigned int readno);

private:
	DECL_LOGGER("BayesQualityGenome")
};

}

#endif /* BAYES_QUALITY_H_ */

