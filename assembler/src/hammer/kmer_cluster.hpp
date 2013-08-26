//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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
#include "kmer_data.hpp"

#include <string>
#include <vector>

class KMerClustering {
public:
  KMerClustering(KMerData &data, unsigned nthreads, const std::string &workdir) :
      data_(data), nthreads_(nthreads), workdir_(workdir) { }

  /**
    * perform k-mer clustering and store the results in the map and the set
    */
  void process(std::vector<std::vector<size_t> > &classes);

private:
  KMerData &data_;
  unsigned nthreads_;
  std::string workdir_;

  /// @return consensus string for a block
  hammer::KMer Consensus(const std::vector<size_t> & block) const;

  hammer::KMer ConsensusWithMask(const std::vector<size_t> & block, const std::vector<size_t> & mask, size_t maskVal) const;

  double ClusterBIC(const std::vector<size_t> & cl, const std::vector<StringCount> & centers, const std::vector<size_t> & indices) const;

  /**
    * perform l-means clustering on the set of k-mers with initial centers being the l most frequent k-mers here
    * @param indices fill array centers with cluster centers; centers[k].count shows how many different kmers are in this cluster (used later)
    * @param centers fill array indices with ints from 0 to l that denote which kmers belong where
    * @return the resulting likelihood of this clustering
    */
  double lMeansClustering(unsigned l, const std::vector<size_t> & kmerinds, std::vector<size_t> & indices, std::vector<StringCount> & centers);

  size_t SubClusterSingle(const std::vector<size_t> & block, std::vector< std::vector<size_t> > & vec);

  std::string GetGoodKMersFname() const;
  std::string GetBadKMersFname() const;

private:
  DECL_LOGGER("Hamming Subclustering");
};


#endif // KMER_CLUSTER_HPP

