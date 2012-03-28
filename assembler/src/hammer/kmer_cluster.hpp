/*
 * kmer_cluster.hpp
 *
 *  Created on: 16.07.2011
 *      Author: snikolenko
 */

#ifndef KMER_CLUSTER_HPP
#define KMER_CLUSTER_HPP

#include <string>
#include <vector>
#include "subkmers.hpp"
#include "position_kmer.hpp"

class unionFindClass;

class KMerClustering {
public:
	KMerClustering(std::vector<KMerCount> & kmers, int nthreads, int tau) : k_(kmers), v_(NULL), nthreads_(nthreads), tau_(tau) { hintVector = false; }
	KMerClustering(std::vector<KMerCount> * kmers_empty, std::vector<hint_t> * kmers, int nthreads, int tau) : k_(*kmers_empty), v_(kmers), nthreads_(nthreads), tau_(tau) { hintVector = true; }

	/**
	  * perform k-mer clustering and store the results in the map and the set
	  */
	void process(bool doHamming, std::string dirprefix, SubKMerSorter * skmsorter, boost::shared_ptr<FOStream> ofs, boost::shared_ptr<FOStream> ofs_bad);

	/// free up memory
	void clear() {
		k_.clear();
	}
	
private:
	std::vector<KMerCount> & k_;
	std::vector<hint_t> * v_;
	int nthreads_;
	int tau_;
	bool hintVector;

	/// @return total log-likelihood that x was made from center given x's quality values
	double logLikelihoodKMer(const string & center, const KMerCount & x);
	double logLikelihoodSingleton(const KMerCount & x);

	/// @return Hamming distance between x and y with upper bound tau
	int hamdistKMer(const PositionKMer & x, const PositionKMer & y, int tau = K);
	/// @return Hamming distance between x and y with upper bound tau
	int hamdistKMer(const PositionKMer & x, const string & y, int tau = K);
	/// @return Hamming distance between x and y with upper bound tau
	int hamdistKMer(const string & x, const string & y, int tau = K);
	/// @return Hamming distance between x and y with upper bound tau
	int hamdistKMer(const hint_t & x, const hint_t & y, int tau = K);
	/// @return consensus string for a block
	std::string find_consensus(const std::vector<int> & block);

	/**
	  * find consensus with mask
	  * @param mask is a vector of integers of the same size as the block
	  * @param maskVal is the integer that we use
	  */
	std::string find_consensus_with_mask(const std::vector<int> & block, const std::vector<int> & mask, int maskVal);
	
	/**
	  * @return total log-likelihood of this particular clustering with real quality values
	  */
	double trueClusterLogLikelihood(const vector<int> & cl, const vector<StringCount> & centers, const vector<int> & indices);
	double trueSingletonLogLikelihood(const hint_t & kmerind);

	/**
	  * perform l-means clustering on the set of k-mers with initial centers being the l most frequent k-mers here
	  * @param indices fill array centers with cluster centers; centers[k].count shows how many different kmers are in this cluster (used later)
	  * @param centers fill array indices with ints from 0 to l that denote which kmers belong where
	  * @return the resulting likelihood of this clustering
	  */
	double lMeansClustering(uint32_t l, std::vector< std::vector<int> > & distances, const std::vector<int> & kmerinds, std::vector<int> & indices, std::vector<StringCount> & centers);

	/**
	  * SIN
	  * new version of process_block
	  * @param newBlockNum current number of a new cluster; incremented inside
	  * @return new value of newBlockNum
	  */
	void process_block_SIN(const std::vector<int> & block, std::vector< std::vector<int> > & vec);

	void processBlock(unionFindClass * uf, vector<hint_t> & block, int cur_subkmer, bool fileBased);
	void processBlockQuadratic(unionFindClass * uf, vector<hint_t> & block, bool direct);
	void clusterMerge(std::vector<unionFindClass *> uf, unionFindClass * ufMaster);

};


#endif // KMER_CLUSTER_HPP

