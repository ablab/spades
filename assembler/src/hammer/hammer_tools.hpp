/*
 * hammer_tools.hpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#ifndef HAMMER_TOOLS_HPP
#define HAMMER_TOOLS_HPP

#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include "sequence/seq.hpp"
#include "read/read.hpp"
#include "union.hpp"
#include "hammerread.hpp"
#include "hammer_config.hpp"

using namespace std;

#define MAX_INT_64 1000000000000000000
#define READ_BATCH_SIZE 10000000
#define ERROR_RATE 0.01

double oct2phred(string qoct, int qvoffset);
string encode3toabyte (const string & s);

class BadConversion : public std::runtime_error {
  public:
	BadConversion(std::string const& s) : std::runtime_error(s) { }
};

inline double convertToInt(std::string const& s) {
   istringstream i(s);
   int x;
   if (!(i >> x))
     throw BadConversion("convertToInt(\"" + s + "\")");
   return x;
}

/// add k-mers from read to map
void addKMers(const Read & r, KMerStatMap & v);

/// join two maps
void join_maps(KMerStatMap & v1, const KMerStatMap & v2);

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

void DoPreprocessing(int tau, int qvoffset, string readsFilename, int nthreads, vector<KMerStatMap> & vv);
void DoSplitAndSort(int tau, int nthreads, ReadStatMapContainer & rmsc, vector<StringKMerVector> * vs, vector<KMerCount> * kmers);

void ClusterProcessBlock(unionFindClass * uf, StringKMerVector & block, int tau);
void ClusterMerge(vector<unionFindClass *> uf, const vector<KMerCount> & kmers, unionFindClass * ufMaster);
void DoClustering(int tau, int nthreads, string dirprefix, const vector<StringKMerVector> & vs, const vector<KMerCount> & kmers);


/// multinomial coefficient
double calcMultCoef(vector<int> & distances, vector<HammerRead> & kmers);
/// find consensus
string find_consensus(vector<HammerRead> & block);
/**
  * SIN
  * find consensus with mask
  * @param mask is a vector of integers of the same size as the block
  * @param maskVal is the integer that we use
  */
string find_consensus_with_mask(vector<HammerRead> & block, const vector<int> & mask, int maskVal);

/**
  * @return total log-likelihood of this particular clustering
  */
double clusterLogLikelihood(const vector<HammerRead> & block, const vector<HammerRead> & centers, const vector<int> & indices);

/**
  * perform l-means clustering on the set of k-mers with initial centers being the l most frequent k-mers here
  * @param indices fill array centers with cluster centers; centers[k].count shows how many different kmers are in this cluster (used later)
  * @param centers fill array indices with ints from 0 to l that denote which kmers belong where
  * @return the resulting likelihood of this clustering
  */
double lMeansClustering(int l, vector< vector<int> > & distances, vector<HammerRead> & kmers, vector<int> & indices, vector<HammerRead> & centers);

/**
  * SIN
  * new version of process_block
  * @param newBlockNum current number of a new cluster; incremented inside
  * @return new value of newBlockNum
  */
void process_block_SIN(vector<HammerRead> & block, vector< vector<HammerRead> > & vec);

#endif



