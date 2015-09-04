//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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

#include <boost/numeric/ublas/fwd.hpp>

class KMerClustering {
public:
  KMerClustering(KMerData &data, unsigned nthreads, const std::string &workdir, bool debug) :
      data_(data), nthreads_(nthreads), workdir_(workdir), debug_(debug) { }

  void process(const std::string &Prefix);

private:
  KMerData &data_;
  unsigned nthreads_;
  std::string workdir_;
  bool debug_;

  struct Center {
    hammer::ExpandedSeq center_;
    size_t count_;
  };
    
  double ClusterBIC(const std::vector<Center> &centers,
                    const std::vector<size_t> &indices, const std::vector<hammer::ExpandedKMer> &kmers) const;

  /**
    * perform l-means clustering on the set of k-mers with initial centers being the l most frequent k-mers here
    * @param indices fill array centers with cluster centers; centers[k].count shows how many different kmers are in this cluster (used later)
    * @param centers fill array indices with ints from 0 to l that denote which kmers belong where
    * @return the resulting likelihood of this clustering
    */
  double lMeansClustering(unsigned l, const std::vector<hammer::ExpandedKMer> &kmers,
                          std::vector<size_t> & indices, std::vector<Center> & centers);

  size_t SubClusterSingle(const std::vector<size_t> & block, std::vector< std::vector<size_t> > & vec);

  std::string GetGoodKMersFname() const;
  std::string GetBadKMersFname() const;

  size_t ProcessCluster(const std::vector<size_t> &cur_class,
                        boost::numeric::ublas::matrix<uint64_t> &errs,
                        std::ofstream &ofs, std::ofstream &ofs_bad,
                        size_t &gsingl, size_t &tsingl, size_t &tcsingl, size_t &gcsingl,
                        size_t &tcls, size_t &gcls, size_t &tkmers, size_t &tncls);

private:
  DECL_LOGGER("Hamming Subclustering");
};


#endif // KMER_CLUSTER_HPP

