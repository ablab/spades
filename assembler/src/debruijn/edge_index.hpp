//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "openmp_wrapper.h"

#include "omni/omni_utils.hpp"

#include "kmer_map.hpp"
#include "new_debruijn.hpp"
#include "standard.hpp"
#include "debruijn_kmer_index.hpp"

namespace debruijn_graph {
/**
 * DataHashRenewer listens to add/delete events and updates index according to those events. This class
 * can be used both with vertices and edges of graph.
 * todo EdgeNucls are hardcoded!
 */
template<typename Graph, typename ElementId>
class DataHashRenewer {
	typedef runtime_k::RtSeq Kmer;
	typedef DeBruijnKMerIndex<EdgeId> Index;

	const Graph &g_;
	Index &index_;

	const bool ignore_new_kmers_;

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewKmersHash(ElementId id) {
		Sequence nucls = g_.EdgeNucls(id);
		//		DEBUG("Renewing hashes for k-mers of sequence " << nucls);
		index_.RenewKMers(nucls, id, ignore_new_kmers_);
	}

	void DeleteKmersHash(ElementId id) {
		Sequence nucls = g_.EdgeNucls(id);
		//		DEBUG("Deleting hashes for k-mers of sequence " << nucls);
		index_.DeleteKMers(nucls, id);
	}

public:
	/**
	 * Creates DataHashRenewer for specified graph and index
	 * @param g graph to be indexed
	 * @param index index to be synchronized with graph
	 */
	DataHashRenewer(const Graph& g, Index& index, bool ignore_new_kmers) :
			g_(g), index_(index), ignore_new_kmers_(ignore_new_kmers) {
	}

	virtual ~DataHashRenewer() {

	}

	void HandleAdd(ElementId id) {
		RenewKmersHash(id);
	}

	virtual void HandleDelete(ElementId id) {
		DeleteKmersHash(id);
	}

private:
	DECL_LOGGER("DataHashRenewer")
};

/**
 * EdgeIndex is a structure to store info about location of certain k-mers in graph. It delegates all
 * container procedures to inner_index_ which is DeBruijnKMerIndex and all handling procedures to
 * renewer_ which is DataHashRenewer.
 * @see DeBruijnKMerIndex
 * @see DataHashRenewer
 */
template<class Graph>
class EdgeIndex: public GraphActionHandler<Graph> {

public:
	typedef typename Graph::EdgeId EdgeId;
  typedef DeBruijnKMerIndex<EdgeId> InnerIndex;
  typedef runtime_k::RtSeq Kmer;

private:
	InnerIndex inner_index_;
	DataHashRenewer<Graph, EdgeId> renewer_;
	bool delete_index_;

public:

  EdgeIndex(const Graph& g, size_t k, const std::string &workdir) :
      GraphActionHandler<Graph>(g, "EdgeIndex"), inner_index_(k, workdir),
      renewer_(g, inner_index_, true), delete_index_(true) {
	}

	virtual ~EdgeIndex() {
		TRACE("~EdgeIndex OK")
	}

	InnerIndex &inner_index() {
		VERIFY(this->IsAttached());
		return inner_index_;
	}

	virtual void HandleAdd(EdgeId e) {
		renewer_.HandleAdd(e);
	}

	virtual void HandleDelete(EdgeId e) {
		renewer_.HandleDelete(e);
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
	}

	bool contains(const Kmer& kmer) const {
		VERIFY(this->IsAttached());
		return inner_index_.ContainsInIndex(kmer);
	}

	const pair<EdgeId, size_t> get(const Kmer& kmer) const {
		VERIFY(this->IsAttached());
		return inner_index_.get(kmer);
	}

	void Refill() {
		clear();
    DeBruijnKMerIndexBuilder().BuildIndexFromGraph(inner_index_, this->g());
	}

	void clear() {
		inner_index_.clear();
	}

};
}
