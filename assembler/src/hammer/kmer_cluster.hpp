//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * kmer_cluster.hpp
 *
 *  Created on: 16.07.2011
 *      Author: snikolenko
 */

#ifndef KMER_CLUSTER_HPP
#define KMER_CLUSTER_HPP

#include "hamcluster.hpp"
#include "kmer_index.hpp"

#include <string>
#include <vector>

class KMerClustering {
public:
	KMerClustering(KMerData &data, unsigned nthreads, const std::string &workdir) :
      data_(data), nthreads_(nthreads), workdir_(workdir) { }

	/**
	  * perform k-mer clustering and store the results in the map and the set
	  */
	void process(std::vector<std::vector<unsigned> > classes);
	
private:
	KMerData &data_;
	unsigned nthreads_;
  std::string workdir_;

	/// @return consensus string for a block
	KMer find_consensus(const std::vector<unsigned> & block);

	/**
	  * find consensus with mask
	  * @param mask is a vector of integers of the same size as the block
	  * @param maskVal is the integer that we use
	  */
	KMer find_consensus_with_mask(const std::vector<unsigned> & block, const std::vector<unsigned> & mask, unsigned maskVal);
	
	double ClusterBIC(const vector<unsigned> & cl, const vector<StringCount> & centers, const vector<unsigned> & indices);

	/**
	  * perform l-means clustering on the set of k-mers with initial centers being the l most frequent k-mers here
	  * @param indices fill array centers with cluster centers; centers[k].count shows how many different kmers are in this cluster (used later)
	  * @param centers fill array indices with ints from 0 to l that denote which kmers belong where
	  * @return the resulting likelihood of this clustering
	  */
	double lMeansClustering(unsigned l, const std::vector<unsigned> & kmerinds, std::vector<unsigned> & indices, std::vector<StringCount> & centers);

	/**
	  * SIN
	  * new version of process_block
	  * @param newBlockNum current number of a new cluster; incremented inside
	  * @return new value of newBlockNum
	  */
	size_t process_block_SIN(const std::vector<unsigned> & block, std::vector< std::vector<unsigned> > & vec);

  std::string GetGoodKMersFname() const;
  std::string GetBadKMersFname() const;
  
private:
  DECL_LOGGER("Hamming Subclustering");
};


#endif // KMER_CLUSTER_HPP

