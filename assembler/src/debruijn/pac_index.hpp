/*
 * pac_index.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: lab42
 */

#ifndef PAC_INDEX_HPP_
#define PAC_INDEX_HPP_
#include "debruijn_kmer_index.hpp"
#include "graph_pack.hpp"


template<class Graph>
struct position{
	typename Graph::EdgeId edge;
	int offset_;
};
template<class Graph>
class PacBioMappingIndex{
  protected :
	Graph &g_;
	size_t K_;
	DeBruijnKMerIndex<typename Graph::EdgeId> tmp_index;

  public :
	typedef position<typename Graph::EdgeId> MappedPosition;
	typedef std::vector<MappedPosition> MappedPositionsVector;

	PacBioMappingIndex(Graph &g, size_t k):g_(g), K_(k), tmp_index(K_, "tmp") {
		DeBruijnKMerIndexBuilder<runtime_k::RtSeq> builder;
		builder.BuildIndexFromGraph<typename Graph::EdgeId, Graph>(tmp_index, g_);

	}

	int Count(Sequence s);
	MappedPositionsVector Locate(Sequence s);
};

template<class Graph>
int PacBioMappingIndex<Graph>::Count(Sequence s){
	runtime_k::RtSeq kmer = s.start<runtime_k::RtSeq>(K_);
	buffer[this->GetFileNumForSeq(kmer, num_files)].push_back(kmer);
	for (size_t j = K_; j < s.size(); ++j) {
	kmer <<= seq[j];
	buffer[this->GetFileNumForSeq(kmer, num_files)].push_back(kmer);
	}

	kmers += ((seq.size() - K_ + 1) + 1);

	tmp_index.contains()
	return 0;
}

template<class Graph>
typename PacBioMappingIndex<Graph>::MappedPositionsVector PacBioMappingIndex<Graph>::Locate(Sequence s){
	vector<MappedPosition> res;
	return res;
}

#endif /* PAC_INDEX_HPP_ */
