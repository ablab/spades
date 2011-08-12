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
#include "sequence/seq.hpp"
#include "kmer_stat.hpp"
#include "position_kmer.hpp"
#include "subkmers.hpp"

using namespace std;

#define MAX_INT_64 1000000000000000000
#define READ_BATCH_SIZE 10000000
#define ERROR_RATE 0.01
#define BLOCKSIZE_QUADRATIC_THRESHOLD 100
#define GOOD_SINGLETON_THRESHOLD 1 
#define CONSENSUS_BLOB_MARGIN 0.25

#define TIMEDLN(a) print_time(); cout << a << endl

double oct2phred(string qoct, int qvoffset);
string encode3toabyte (const string & s);
void print_time();

/// join two maps
void join_maps(KMerStatMap & v1, const KMerStatMap & v2);

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const PositionRead &r, hint_t readno, KMerStatMap *v);

void AddKMerNos(const PositionRead &r, hint_t readno, vector<KMerNo> *v);

void DoPreprocessing(int tau, int qvoffset, string readsFilename, int nthreads, vector<KMerNo> * vv);
void ParallelSortKMerNos(vector<KMerNo> * v, vector<KMerCount> * kmers, int nthreads);
void DoSplitAndSort(int tau, int nthreads, vector< vector<hint_t> > * vs, vector<KMerCount> * kmers, vector<SubKMerPQ> * vskpq);

/**
  * correct a read in place
  * @return whether the read has changed at all
  */
bool CorrectRead(const vector<KMerCount> & kmers, hint_t readno, ofstream * ofs = NULL);

#endif

