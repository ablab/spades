/*
 * hammer_tools.hpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#ifndef HAMMER_TOOLS_HPP
#define HAMMER_TOOLS_HPP

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include "read/read.hpp"
#include "read/ireadstream.hpp"
#include "union.hpp"
#include "hammer_config.hpp"
#include "sequence/seq.hpp"
#include "kmer_stat.hpp"


using namespace std;

#define MAX_INT_64 1000000000000000000
#define READ_BATCH_SIZE 10000000
#define ERROR_RATE 0.01

double oct2phred(string qoct, int qvoffset);
string encode3toabyte (const string & s);

/// join two maps
void join_maps(KMerStatMap & v1, const KMerStatMap & v2);

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const Read &r, uint64_t readno, KMerStatMap *v);

class ReadStatMapContainer {
public:
	ReadStatMapContainer(const vector<KMerStatMap> & vv) : v_(vv) { init(); }

	void init();
	KMerCount next();
	size_t size();
	
private:
	const vector<KMerStatMap> & v_;
	vector<KMerStatMap::const_iterator> i_;

	const KMerStatMap::const_iterator & cur_min();
};

void DoPreprocessing(int tau, int qvoffset, string readsFilename, int nthreads, vector<KMerStatMap> * vv, vector<ReadStat> * rv);
void DoSplitAndSort(int tau, int nthreads, ReadStatMapContainer & rmsc, vector<StringKMerVector> * vs, vector<KMerCount> * kmers, vector<ReadStat> * rv);

/**
  * correct a read in place
  * @return whether the read has changed at all
  */
bool CorrectRead(const vector<KMerCount> & kmers, ReadStat * r, ofstream * ofs = NULL);

#endif

