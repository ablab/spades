//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "io/paired_read.hpp"
#include "seq_map.hpp"
#include "omni/omni_utils.hpp"
#include "omni/id_track_handler.hpp"
#include "logger/logger.hpp"
#include "omni/paired_info.hpp"
#include "xmath.h"
#include <boost/optional.hpp>
#include <iostream>
#include "sequence/sequence_tools.hpp"
#include "omni/splitters.hpp"
#include "openmp_wrapper.h"
#include "omni/path_processor.hpp"
#include "omni/distance_estimation.hpp"

#include "runtime_k.hpp"

#include "kmer_map.hpp"
#include "new_debruijn.hpp"
#include "standard.hpp"
#include "debruijn_kmer_index.hpp"

#include "path_helper.hpp"

using namespace std;

namespace debruijn_graph {

using omnigraph::Path;
using omnigraph::MappingPath;
using omnigraph::Range;
using omnigraph::MappingRange;
using omnigraph::PairInfo;
using omnigraph::GraphActionHandler;

template<typename Graph>
class ReadThreaderResult {
	typedef typename Graph::EdgeId EdgeId;

	Path<EdgeId> left_read_, right_read_;
	int gap_;
public:
	ReadThreaderResult(Path<EdgeId> left_read, Path<EdgeId> right_read, int gap) :
			gap_(gap), left_read_(left_read), right_read_(right_read) {
	}
};

template<typename Graph>
class SingleReadThreaderResult {
	typedef typename Graph::EdgeId EdgeId;
public:
	EdgeId edge_;
	int read_position_;
	int edge_position_;
	SingleReadThreaderResult(EdgeId edge, int read_position, int edge_position) :
			edge_(edge), read_position_(read_position), edge_position_(
					edge_position) {
	}
};

template<typename Graph>
class ReadMappingResult {
public:
	Sequence read_;
	vector<SingleReadThreaderResult<Graph> > res_;
	ReadMappingResult(Sequence read,
			vector<SingleReadThreaderResult<Graph> > res) :
			read_(read), res_(res) {

	}
	ReadMappingResult() {

	}
};

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

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewKmersHash(ElementId id) {
		Sequence nucls = g_.EdgeNucls(id);
		//		DEBUG("Renewing hashes for k-mers of sequence " << nucls);
		index_.RenewKMers(nucls, id);
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
	DataHashRenewer(const Graph& g, Index& index) :
			g_(g), index_(index) {
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
 * container procedures to inner_index_ which is SeqMap and all handling procedures to
 * renewer_ which is DataHashRenewer.
 * @see SeqMap
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

	EdgeIndex(const Graph& g, size_t k) :
			GraphActionHandler<Graph>(g, "EdgeIndex"), inner_index_(k), renewer_(
					g, inner_index_), delete_index_(true) {
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
    DeBruijnKMerIndexBuilder().BuildIndexFromGraph(inner_index_,
                                                   this->g());
	}

	void clear() {
		inner_index_.clear();
	}

};

template<class Graph>
class KmerMapper: public omnigraph::GraphActionHandler<Graph> {
	typedef omnigraph::GraphActionHandler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef runtime_k::RtSeq Kmer;
	typedef typename runtime_k::KmerMap<Kmer> MapType;

	MapType mapping_;

	size_t k_;

public:

	KmerMapper(const Graph& g, size_t k): base(g, "KmerMapper"), mapping_(k), k_(k) {
	}

	virtual ~KmerMapper() {

	}

	size_t get_k() const {
	    return k_;
	}

	typename MapType::const_iterator begin() const {
		return mapping_.begin();
	}

	typename MapType::const_iterator end() const {
		return mapping_.end();
	}

	void Normalize() {
		vector<Kmer> all;
		for(auto it = begin(); it != end(); ++it) {
			all.push_back(it->first);
		}
		for(auto it = all.begin(); it != all.end(); ++it) {
			Normalize(*it);
		}
	}

	void Revert(const Kmer &kmer) {
		Kmer old_value = Substitute(kmer);
		if(old_value != kmer) {
			mapping_.erase(kmer);
			mapping_[old_value] = kmer;
		}
	}

	void Normalize(const Kmer &kmer) {
		mapping_[kmer] = Substitute(kmer);
	}

	bool CheckCanRemap(const Sequence& old_s, const Sequence& new_s) const {
		size_t old_length = old_s.size() - k_ + 1;
		size_t new_length = new_s.size() - k_ + 1;
		UniformPositionAligner aligner(old_s.size() - k_ + 1,
				new_s.size() - k_ + 1);
		Kmer old_kmer = old_s.start<Kmer::max_size>(k_);
		old_kmer >>= 0;

		for (size_t i = k_ - 1; i < old_s.size(); ++i) {
			old_kmer <<= old_s[i];
			size_t old_kmer_offset = i - k_ + 1;
			size_t new_kmer_offest = aligner.GetPosition(old_kmer_offset);
			if(old_kmer_offset * 2 + 1 == old_length && new_length % 2 == 0) {
				Kmer middle(k_ - 1, new_s, new_length / 2);
				if(typename Kmer::less2()(middle, !middle)) {
					new_kmer_offest = new_length - 1 - new_kmer_offest;
				}
			}
			Kmer new_kmer(k_, new_s, new_kmer_offest);
			auto it = mapping_.find(new_kmer);
			if(it != mapping_.end()) {
				if(Substitute(new_kmer) != old_kmer) {
					return false;
				}
			}
		}
		return true;
	}

	void RemapKmers(const Sequence& old_s, const Sequence& new_s) {
		VERIFY(this->IsAttached());
		//		cout << endl << "Mapping " << old_s << " to " << new_s << endl;
		size_t old_length = old_s.size() - k_ + 1;
		size_t new_length = new_s.size() - k_ + 1;
		UniformPositionAligner aligner(old_s.size() - k_ + 1,
				new_s.size() - k_ + 1);
		Kmer old_kmer = old_s.start<Kmer::max_size>(k_);
		old_kmer >>= 0;

		for (size_t i = k_ - 1; i < old_s.size(); ++i) {
			old_kmer <<= old_s[i];
			size_t old_kmer_offset = i - k_ + 1;
			size_t new_kmer_offest = aligner.GetPosition(old_kmer_offset);
			if(old_kmer_offset * 2 + 1 == old_length && new_length % 2 == 0) {
				Kmer middle(k_ - 1, new_s, new_length / 2);
				if(typename Kmer::less2()(middle, !middle)) {
					new_kmer_offest = new_length - 1 - new_kmer_offest;
				}
			}
			Kmer new_kmer(k_, new_s, new_kmer_offest);
			auto it = mapping_.find(new_kmer);
			if(it != mapping_.end()) {
				VERIFY(Substitute(new_kmer) == old_kmer);
				mapping_.erase(it);
			}
			mapping_[old_kmer] = new_kmer;
		}
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		VERIFY(this->g().EdgeNucls(new_edge) == this->g().EdgeNucls(edge2));
		RemapKmers(this->g().EdgeNucls(edge1), this->g().EdgeNucls(edge2));
	}

	Kmer Substitute(const Kmer& kmer) const {
		VERIFY(this->IsAttached());
		Kmer answer = kmer;
		auto it = mapping_.find(answer);
		while (it != mapping_.end()) {
			VERIFY(it.first() != it.second());
			answer = it.second();
			it = mapping_.find(answer);
		}
		return answer;
	}

	void BinWrite(std::ostream& file) const {
		u_int32_t size = mapping_.size();
		file.write((const char *) &size, sizeof(u_int32_t));

		for (auto iter = mapping_.begin(); iter != mapping_.end(); ++iter) {
			Kmer::BinWrite(file, iter.first());
			Kmer::BinWrite(file, iter.second());
		}
	}

	void BinRead(std::istream& file) {
		mapping_.clear();
		u_int32_t size;
		file.read((char *) &size, sizeof(u_int32_t));

		for (u_int32_t i = 0; i < size; ++i) {
			Kmer key(k_);
			Kmer value(k_);
			Kmer::BinRead(file, &key);
			Kmer::BinRead(file, &value);
			mapping_[key] = value;
		}
	}

	bool CompareTo(KmerMapper<Graph> const& m) {
		if (mapping_.size() != m.mapping_.size()) {
			INFO("Unequal sizes");
		}
		for (auto iter = mapping_.begin(); iter != mapping_.end(); ++iter) {
			auto cmp = m.mapping_.find(iter.first());
			if (cmp == m.mapping_.end() || cmp.second() != iter.second()) {
				return false;
			}
		}
		return true;
	}

	void clear() {
		mapping_.clear();
	}

	size_t size() const {
		return mapping_.size();
	}
};

/**
 * This class finds how certain sequence is mapped to genome. As it is now it works correct only if sequence
 * is mapped to graph ideally and in unique way.
 */
template<class Graph>
class SimpleSequenceMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef runtime_k::RtSeq Kmer;
	typedef EdgeIndex<Graph> Index;

private:
	const Graph& g_;
	const Index &index_;
	size_t k_;

	bool TryThread(Kmer &kmer, vector<EdgeId> &passed,
			size_t& endPosition) const {
        VERIFY(passed.size() > 0);
		EdgeId last = passed[passed.size() - 1];
        if (endPosition + 1 < g_.length(last)) {
            if (g_.EdgeNucls(last)[endPosition + k_] == kmer[k_ - 1]) {
                endPosition++;
                return true;
            }
        } else {
            vector<EdgeId> edges = g_.OutgoingEdges(g_.EdgeEnd(last));
            for (size_t i = 0; i < edges.size(); i++) {
                if (g_.EdgeNucls(edges[i])[k_ - 1] == kmer[k_ - 1]) {
                    passed.push_back(edges[i]);
                    endPosition = 0;
                    return true;
                }
            }
        }
		return false;
	}

	bool FindKmer(Kmer& kmer, vector<EdgeId> &passed, size_t& startPosition,
			size_t& endPosition) const {
        //TRACE("CONTAINS kmer " << " " << omp_get_thread_num() );
        if (index_.contains(kmer)) {
            //TRACE("YES CONTAINS " << omp_get_thread_num());
            pair<EdgeId, size_t> position = index_.get(kmer);
            //DEBUG("LENGTH " << g_.length(position.first));
            endPosition = position.second;
            if (passed.empty()) {
                startPosition = position.second;
            }
            if (passed.empty() || passed.back() != position.first) {
                passed.push_back(position.first);
            }
            
			return true;
		}
		return false;
	}

	bool ProcessKmer(Kmer &kmer, vector<EdgeId> &passed, size_t &startPosition,
			size_t &endPosition, bool valid) const {
		//DEBUG("process kmer started " << omp_get_thread_num() << " valid " << valid);
		if (valid) {
			return TryThread(kmer, passed, endPosition);
		} else {
            return FindKmer(kmer, passed, startPosition, endPosition);
		}
        return false;
        //DEBUG("process kmer finished " << omp_get_thread_num());
	}

public:
	/**
	 * Creates SimpleSequenceMapper for given graph. Also requires index_ which should be synchronized
	 * with graph.
	 * @param g graph sequences should be mapped to
	 * @param index index synchronized with graph
	 */
	SimpleSequenceMapper(const Graph& g, const Index& index, size_t k) :
			g_(g), index_(index), k_(k) {
	}

	/**
	 * Finds a path in graph which corresponds to given sequence.
	 * @read sequence to be mapped
	 */

	Path<EdgeId> MapSequence(const Sequence &read) const {
		vector<EdgeId> passed;
        //TRACE("Mapping sequence");
		if (read.size() <= k_ - 1) {
			return Path<EdgeId>();
		}

		Kmer kmer = read.start<Kmer::max_size>(k_);
        //DEBUG("started " << kmer.str() << omp_get_thread_num() );
		size_t startPosition = -1;
		size_t endPosition = -1;
		bool valid = ProcessKmer(kmer, passed, startPosition, endPosition,
				false);
        for (size_t i = k_; i < read.size(); ++i) {
            kmer <<= read[i];
            //DEBUG("shifted " << kmer.str() << omp_get_thread_num());
            valid = ProcessKmer(kmer, passed, startPosition, endPosition,
                    valid);
        }
        //DEBUG("Path got " << omp_get_thread_num());
        Path<EdgeId> ans(passed, startPosition, endPosition + 1);
		return ans;
	}

};

template<class Graph>
class ExtendedSequenceMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef vector<MappingRange> RangeMappings;
	typedef runtime_k::RtSeq Kmer;
	typedef EdgeIndex<Graph> Index;
	typedef KmerMapper<Graph> KmerSubs;

private:
	const Graph& g_;
	const Index& index_;
	const KmerSubs& kmer_mapper_;
	size_t k_;

	void FindKmer(Kmer kmer, size_t kmer_pos, vector<EdgeId> &passed,
			RangeMappings& range_mappings) const {

		if (index_.contains(kmer)) {
			pair<EdgeId, size_t> position = index_.get(kmer);
			if (passed.empty() || passed.back() != position.first
					|| kmer_pos != range_mappings.back().initial_range.end_pos
					|| position.second + 1
							< range_mappings.back().mapped_range.end_pos) {
				passed.push_back(position.first);
				MappingRange mapping_range(Range(kmer_pos, kmer_pos + 1),
						Range(position.second, position.second + 1));
				range_mappings.push_back(mapping_range);
			} else {
				range_mappings.back().initial_range.end_pos = kmer_pos + 1;
				range_mappings.back().mapped_range.end_pos = position.second
						+ 1;
			}
		}
	}

	void ProcessKmer(Kmer kmer, size_t kmer_pos, vector<EdgeId> &passed,
			RangeMappings& interval_mapping) const {
		kmer = kmer_mapper_.Substitute(kmer);
		FindKmer(kmer, kmer_pos, passed, interval_mapping);
	}

public:
	ExtendedSequenceMapper(const Graph& g,
			const Index& index,
			const KmerSubs& kmer_mapper,
			size_t k) :
			g_(g), index_(index), kmer_mapper_(kmer_mapper), k_(k) {
	}

	MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
		vector<EdgeId> passed_edges;
		RangeMappings range_mapping;

		if (sequence.size() < k_) {
			return MappingPath<EdgeId>();
		}
		Kmer kmer = sequence.start<Kmer::max_size>(k_);
		kmer >>= 0;
		for (size_t i = k_ - 1; i < sequence.size(); ++i) {
			kmer <<= sequence[i];
			ProcessKmer(kmer, i - k_ + 1, passed_edges, range_mapping);
		}

		//DEBUG
//		for (size_t i = 0; i < passed_edges.size(); ++i) {
//			cerr << int_ids_.ReturnIntId(passed_edges[i]) << " (" << range_mapping[i] << ")"<< "; ";
//		}
//		cerr << endl;
		//DEBUG

		return MappingPath<EdgeId>(passed_edges, range_mapping);
	}
};

//todo compare performance
template<class Graph>
class NewExtendedSequenceMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef vector<MappingRange> RangeMappings;
	typedef runtime_k::RtSeq Kmer;
	typedef EdgeIndex<Graph> Index;
	typedef KmerMapper<Graph> KmerSubs;

private:
	const Graph& g_;
	const Index& index_;
	const KmerSubs& kmer_mapper_;
	size_t k_;
//	mutable size_t mapped_;
//	mutable size_t unmapped_;

	bool FindKmer(Kmer kmer, size_t kmer_pos, vector<EdgeId> &passed,
			RangeMappings& range_mappings) const {
		if (index_.contains(kmer)) {
			pair<EdgeId, size_t> position = index_.get(kmer);
			if (passed.empty() || passed.back() != position.first
					|| kmer_pos != range_mappings.back().initial_range.end_pos
					|| position.second + 1
							< range_mappings.back().mapped_range.end_pos) {
				passed.push_back(position.first);
				range_mappings.push_back(
						MappingRange(Range(kmer_pos, kmer_pos + 1),
								Range(position.second, position.second + 1)));
			} else {
				range_mappings.back().initial_range.end_pos = kmer_pos + 1;
				range_mappings.back().mapped_range.end_pos = position.second
						+ 1;
			}
			return true;
		}
		return false;
	}

	bool TryThread(const Kmer& kmer, size_t kmer_pos, vector<EdgeId> &passed,
			RangeMappings& range_mappings) const {
		EdgeId last_edge = passed.back();
		size_t end_pos = range_mappings.back().mapped_range.end_pos;
		if (end_pos < g_.length(last_edge)) {
			if (g_.EdgeNucls(last_edge)[end_pos + k_ - 1] == kmer[k_ - 1]) {
				range_mappings.back().initial_range.end_pos++;
				range_mappings.back().mapped_range.end_pos++;
				return true;
			}
		} else {
			vector<EdgeId> edges = g_.OutgoingEdges(g_.EdgeEnd(last_edge));
			for (size_t i = 0; i < edges.size(); i++) {
				if (g_.EdgeNucls(edges[i])[k_ - 1] == kmer[k_ - 1]) {
					passed.push_back(edges[i]);
					range_mappings.push_back(
							MappingRange(Range(kmer_pos, kmer_pos + 1),
									Range(0, 1)));
					return true;
				}
			}
		}
		return false;
	}

	bool Substitute(Kmer& kmer) const {
		Kmer subs = kmer_mapper_.Substitute(kmer);
		if (subs != kmer) {
			kmer = subs;
			return true;
		}
		return false;
	}

	bool ProcessKmer(Kmer kmer, size_t kmer_pos, vector<EdgeId> &passed_edges,
			RangeMappings& range_mapping, bool try_thread) const {
		if (try_thread) {
			if(!TryThread(kmer, kmer_pos, passed_edges, range_mapping)) {
				Substitute(kmer);
				FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
				return false;
			} else {
				return true;
			}
		} else {
			if(!Substitute(kmer)) {
				return FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
			} else {
				FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
				return false;
			}
		}
//		if (!Substitute(kmer)) {
//			if (try_thread) {
//				return TryThread(kmer, kmer_pos, passed_edges, range_mapping);
//			} else {
//				return FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
//			}
//		} else {
//			FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
//			return false;
//		}
	}

public:
	NewExtendedSequenceMapper(const Graph& g,
			const Index& index,
			const KmerSubs& kmer_mapper,
			size_t k) :
			g_(g), index_(index),
			kmer_mapper_(kmer_mapper), k_(k) {
	}

	~NewExtendedSequenceMapper() {
//		TRACE("In destructor of sequence mapper");
//		TRACE(mapped_ << " sequences were mapped");
//		TRACE(unmapped_ << " sequences couldn't be mapped");
	}

	MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
		vector<EdgeId> passed_edges;
		RangeMappings range_mapping;

		if (sequence.size() < k_) {
			return MappingPath<EdgeId>();
		}

		Kmer kmer = sequence.start<Kmer::max_size>(k_);
		kmer >>= 0;
		bool try_thread = false;
		for (size_t i = k_ - 1; i < sequence.size(); ++i) {
			kmer <<= sequence[i];
			try_thread = ProcessKmer(kmer, i - k_ + 1, passed_edges,
					range_mapping, try_thread);
		}

//		if (passed_edges.empty()) {
////			TRACE("Sequence " << sequence << "couldn't be mapped");
//			unmapped_++;
//			//todo maybe check path consistency?
//		} else {
//			mapped_++;
//		}

		return MappingPath<EdgeId>(passed_edges, range_mapping);
	}

	MappingPath<EdgeId> MapRead(const io::SingleRead &read) const {
		VERIFY(read.IsValid());
		return MapSequence(read.sequence());
	}

private:
	DECL_LOGGER("NewExtendedSequenceMapper");
};

template<class gp_t>
shared_ptr<const NewExtendedSequenceMapper<typename gp_t::graph_t> > MapperInstance(const gp_t& gp) {
    size_t k_plus_1 = gp.k_value + 1;
	return make_shared<NewExtendedSequenceMapper<typename gp_t::graph_t> >(gp.g, gp.index, gp.kmer_mapper, k_plus_1);
}

template<class Graph>
class OldEtalonPairedInfoCounter {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
	const EdgeIndex<Graph>& index_;
	size_t k_;

	size_t insert_size_;
	size_t read_length_;
	size_t gap_;
	size_t delta_;

	void AddEtalonInfo(omnigraph::PairedInfoIndex<Graph>& paired_info,
			EdgeId e1, EdgeId e2, double d) {
		PairInfo<EdgeId> pair_info(e1, e2, d, 1000.0, 0.);
		paired_info.AddPairInfo(pair_info);
	}

	void ProcessSequence(const Sequence& sequence,
			omnigraph::PairedInfoIndex<Graph>& paired_info) {
		SimpleSequenceMapper<Graph> sequence_mapper(g_, index_, k_ + 1);
		Path<EdgeId> path = sequence_mapper.MapSequence(sequence);

		for (size_t i = 0; i < path.size(); ++i) {
			EdgeId e = path[i];
			if (g_.length(e) + delta_ > gap_ + k_ + 1) {
				AddEtalonInfo(paired_info, e, e, 0);
			}
			size_t j = i + 1;
			size_t length = 0;

			while (j < path.size()
					&& length
							<= omnigraph::PairInfoPathLengthUpperBound(k_,
									insert_size_, delta_)) {
				if (length
						>= omnigraph::PairInfoPathLengthLowerBound(k_,
								g_.length(e), g_.length(path[j]), gap_, delta_)) {
					AddEtalonInfo(paired_info, e, path[j],
							g_.length(e) + length);
				}
				length += g_.length(path[j++]);
			}
		}

	}

	/* DEBUG method
	 void CheckPairInfo(const Sequence& genome, omnigraph::PairedInfoIndex<Graph>& paired_info) {
	 SimpleSequenceMapper<k + 1, Graph> mapper(g_, index_);
	 Path<EdgeId> path = mapper.MapSequence(genome);
	 vector<EdgeId> sequence = path.sequence();
	 EdgeId prev = 0;
	 for (auto it = sequence.begin(); it != sequence.end(); ++it) {
	 if (prev != 0) {
	 vector<PairInfo<EdgeId> > infos = paired_info.GetEdgePairInfo(prev, *it);
	 bool imperfect_flag = false;
	 bool perfect_flag = false;
	 for (auto info_it = infos.begin(); info_it != infos.end(); info_it++) {
	 if (abs((*info_it).d - g_.length(prev)) < 2) {
	 imperfect_flag = true;
	 }
	 if (math::eq((*info_it).d, 0. + g_.length(prev))) {
	 perfect_flag = true;
	 break;
	 }
	 }
	 if (!perfect_flag && imperfect_flag) {
	 cerr<< "AAAAAAAAAAAAAAA" <<endl;
	 }
	 }
	 prev = *it;
	 }
	 }*/

public:

	OldEtalonPairedInfoCounter(const Graph& g,
			const EdgeIndex<Graph>& index, size_t insert_size,
			size_t read_length, size_t delta, size_t k) :
			g_(g), index_(index), k_(k), insert_size_(insert_size), read_length_(
					read_length), gap_(insert_size_ - 2 * read_length_), delta_(
					delta) {
		VERIFY(insert_size_ >= 2 * read_length_);
	}

	void FillEtalonPairedInfo(const Sequence& genome,
			omnigraph::PairedInfoIndex<Graph>& paired_info) {
		ProcessSequence(genome, paired_info);
		ProcessSequence(!genome, paired_info);
		//DEBUG
		//		CheckPairInfo(genome, paired_info);
	}
};

//todo rewrite with extended sequence mapper!
template<class Graph>
class EtalonPairedInfoCounter {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
	const EdgeIndex<Graph>& index_;
	const KmerMapper<Graph>& kmer_mapper_;
	size_t k_;

	size_t insert_size_;
	size_t read_length_;
	int gap_;
	size_t delta_;

	void AddEtalonInfo(set<PairInfo<EdgeId>> & paired_info, EdgeId e1, EdgeId e2,
			double d) {
		PairInfo<EdgeId> pair_info(e1, e2, d, 1000.0, 0.);
		paired_info.insert(pair_info);
	}

	void ProcessSequence(const Sequence& sequence,
			set<PairInfo<EdgeId>>& temporary_info) {
		int mod_gap = (gap_ + (int) k_ > (int) delta_ ) ? gap_ - (int) delta_ :int(0) -k_;
		runtime_k::RtSeq left(k_ +1, sequence);
		left >>= 0;
		for (size_t left_idx = 0;
				left_idx + 2 * (k_ + 1) + mod_gap <= sequence.size();
				++left_idx) {
			left <<= sequence[left_idx + k_];
			runtime_k::RtSeq left_upd = kmer_mapper_.Substitute(left);
			if (!index_.contains(left_upd)) {
				continue;
			}
			pair<EdgeId, size_t> left_pos = index_.get(left_upd);

			size_t right_idx = left_idx + k_ + 1 + mod_gap;
			runtime_k::RtSeq right(k_ + 1, sequence, right_idx);
			right >>= 0;
			for (;
					right_idx + k_ + 1 <= left_idx + insert_size_ + delta_
							&& right_idx + k_ + 1 <= sequence.size();
					++right_idx) {
				right <<= sequence[right_idx + k_];
				runtime_k::RtSeq right_upd = kmer_mapper_.Substitute(right);
				if (!index_.contains(right_upd)) {
					continue;
				}
				pair<EdgeId, size_t> right_pos = index_.get(right_upd);

				AddEtalonInfo(
						temporary_info,
						left_pos.first,
						right_pos.first,
						0. + right_idx - left_idx + left_pos.second
								- right_pos.second);
			}
		}
	}

public:

	EtalonPairedInfoCounter(const Graph& g, const EdgeIndex<Graph>& index
			,
			const KmerMapper<Graph>& kmer_mapper,
			size_t insert_size,
			size_t read_length, size_t delta, size_t k) :
			g_(g), index_(index), kmer_mapper_(kmer_mapper), k_(k), insert_size_(
					insert_size), read_length_(read_length), gap_(
					insert_size_ - 2 * read_length_), delta_(delta) {
//		VERIFY(insert_size_ >= 2 * read_length_);
	}

	void FillEtalonPairedInfo(const Sequence& genome,
			omnigraph::PairedInfoIndex<Graph>& paired_info) {
		set<PairInfo<EdgeId>> temporary_info;
		ProcessSequence(genome, temporary_info);
		ProcessSequence(!genome, temporary_info);
		for (auto it = temporary_info.begin(); it != temporary_info.end();
				++it) {
			paired_info.AddPairInfo(*it);
		}
	}
};

double PairedReadCountWeight(const MappingRange&, const MappingRange&) {
	return 1.;
}

double KmerCountProductWeight(const MappingRange& mr1,
		const MappingRange& mr2) {
	return mr1.initial_range.size() * mr2.initial_range.size();
}

ConjugateDeBruijnGraph::EdgeId conj_wrap(ConjugateDeBruijnGraph& g,
		ConjugateDeBruijnGraph::EdgeId e) {
	return g.conjugate(e);
}

NonconjugateDeBruijnGraph::EdgeId conj_wrap(NonconjugateDeBruijnGraph& g,
		NonconjugateDeBruijnGraph::EdgeId e) {
	VERIFY(0);
	return e;
}

void WrappedSetCoverage(ConjugateDeBruijnGraph& g,
		ConjugateDeBruijnGraph::EdgeId e, int cov) {
		g.coverage_index().SetCoverage(e, cov);
		g.coverage_index().SetCoverage(g.conjugate(e), cov);
}

void WrappedSetCoverage(NonconjugateDeBruijnGraph& g,
		NonconjugateDeBruijnGraph::EdgeId e, int cov) {
	g.coverage_index().SetCoverage(e, cov);
}

/**
 * As for now it ignores sophisticated case of repeated consecutive
 * occurrence of edge in path due to gaps in mapping
 *
 * todo talk with Anton about simplification and speed-up of procedure with little quality loss
 */
template<class Graph, class SequenceMapper, class PairedStream>
class LatePairedIndexFiller {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef runtime_k::RtSeq Kmer;
	typedef boost::function<double(MappingRange, MappingRange)> WeightF;
	const Graph& graph_;
	const SequenceMapper& mapper_;
	io::ReadStreamVector< PairedStream > streams_;
	WeightF weight_f_;

	template<class PairedRead>
	void ProcessPairedRead(omnigraph::PairedInfoIndex<Graph>& paired_index,
			const PairedRead& p_r) {
		Sequence read1 = p_r.first().sequence();
		Sequence read2 = p_r.second().sequence();

		MappingPath<EdgeId> path1 = mapper_.MapSequence(read1);
		MappingPath<EdgeId> path2 = mapper_.MapSequence(read2);
		size_t read_distance = p_r.distance();
		for (size_t i = 0; i < path1.size(); ++i) {
			pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
			for (size_t j = 0; j < path2.size(); ++j) {
				pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
				double weight = weight_f_(mapping_edge_1.second,
						mapping_edge_2.second);
				size_t kmer_distance = read_distance
						+ mapping_edge_2.second.initial_range.start_pos
						- mapping_edge_1.second.initial_range.start_pos;
				int edge_distance = kmer_distance
						+ mapping_edge_1.second.mapped_range.start_pos
						- mapping_edge_2.second.mapped_range.start_pos;

				paired_index.AddPairInfo(
						PairInfo<EdgeId>(mapping_edge_1.first,
								mapping_edge_2.first, (double) edge_distance,
								weight, 0.));
			}
		}
	}

    /**
     * Method reads paired data from stream, maps it to genome and stores it in this PairInfoIndex.
     */
    void FillUsualIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
        for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            paired_index.AddPairInfo(PairInfo<EdgeId>(*it, *it, 0, 0.0, 0.));
        }

        INFO("Processing paired reads (takes a while)");

        PairedStream& stream = streams_.back();
        stream.reset();
        size_t n = 0;
        while (!stream.eof()) {
            typename PairedStream::read_type p_r;
            stream >> p_r;
            ProcessPairedRead(paired_index, p_r);
            VERBOSE_POWER(++n, " paired reads processed");
        }
    }


    void FillParallelIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
        for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            paired_index.AddPairInfo(PairInfo<EdgeId>(*it, *it, 0, 0.0, 0.));
        }

        INFO("Processing paired reads (takes a while)");

        size_t nthreads = streams_.size();
        std::vector< omnigraph::PairedInfoIndex<Graph>* > buffer_pi(nthreads);
        buffer_pi[0] = &paired_index;

        for (size_t i = 1; i < nthreads; ++i) {
            buffer_pi[i] = new omnigraph::PairedInfoIndex<Graph>(graph_, paired_index.GetMaxDifference());
        }

        size_t counter = 0;
        #pragma omp parallel num_threads(nthreads)
        {
            #pragma omp for reduction(+ : counter)
            for (size_t i = 0; i < nthreads; ++i) {

                typename PairedStream::read_type r;
                PairedStream& stream = streams_[i];
                stream.reset();

                while (!stream.eof()) {
                    stream >> r;
                    ++counter;

                    ProcessPairedRead(*buffer_pi[i], r);
                }
            }
        }
        INFO("Used " << counter << " paired reads");

        INFO("Merging paired indices");
        for (size_t i = 1; i < nthreads; ++i) {
            buffer_pi[0]->AddAll(*(buffer_pi[i]));
            delete buffer_pi[i];
        }
    }

public:

	LatePairedIndexFiller(const Graph &graph, const SequenceMapper& mapper, PairedStream& stream, WeightF weight_f) :
			graph_(graph), mapper_(mapper), streams_(1, &stream), weight_f_(weight_f)
	{
	}

    LatePairedIndexFiller(const Graph &graph, const SequenceMapper& mapper, const io::ReadStreamVector< PairedStream >& streams, WeightF weight_f) :
            graph_(graph), mapper_(mapper), streams_(streams), weight_f_(weight_f)
    {
    }

    void FillIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
        if (streams_.size() == 1) {
            FillUsualIndex(paired_index);
        } else {
            FillParallelIndex(paired_index);
        }
    }

private:
	DECL_LOGGER("LatePairedIndexFiller");
};

/**
 * This class finds how certain _paired_ read is mapped to genome. As it is now it is hoped to work correctly only if read
 * is mapped to graph ideally and in unique way.
 */
template<class Graph>
class TemplateReadMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef EdgeIndex<Graph> Index;
private:
	SimpleSequenceMapper<Graph> read_seq_mapper;
	io::IReader<io::SingleRead>& stream_;
	size_t k_;
public:
	/**
	 * Creates TemplateReadMapper for given graph. Also requires index_ which should be synchronized
	 * with graph.
	 * @param g graph sequences should be mapped to
	 * @param index index syncronized with graph
	 */
	TemplateReadMapper(const Graph& g, const Index& index, io::IReader<io::SingleRead>& stream, size_t k) :
			read_seq_mapper(g, index, k), stream_(stream), k_(k) {
		stream_.reset();
	}

	/*	ReadThreaderResult<k + 1, Graph> ThreadNext() {
	 if (!stream_.eof()) {
	 io::PairedRead p_r;
	 stream_ >> p_r;
	 Sequence read1 = p_r.first().sequence();
	 Sequence read2 = p_r.second().sequence();
	 Path<EdgeId> aligned_read[2];
	 aligned_read[0] = read_seq_mapper.MapSequence(read1);
	 aligned_read[1] = read_seq_mapper.MapSequence(read2);
	 size_t insert_size = p_r.insert_size();
	 // TODO bug here with insert_size/distance
	 int current_distance1 = insert_size + aligned_read[0].start_pos() - aligned_read[1].start_pos();
	 return ReadThreaderResult<k + 1, Graph> (aligned_read[0],
	 aligned_read[1], current_distance1);
	 }
	 //		else return NULL;
	 }
	 */
};

template<class Graph>
class SingleReadMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef EdgeIndex<Graph> Index;
private:
	SimpleSequenceMapper<Graph> read_seq_mapper;
	const Graph& g_;
	const Index& index_;
	size_t k_;
public:
	/**
	 * Creates SingleReadMapper for given graph. Also requires index_ which should be synchronized
	 * with graph.
	 * @param g graph sequences should be mapped to
	 * @param index index syncronized with graph
	 */
	SingleReadMapper(const Graph& g, const Index& index, size_t k) :
			read_seq_mapper(g, index, k), g_(g), index_(index), k_(k) {
	}

	vector<EdgeId> GetContainingEdges(io::SingleRead& p_r) {
		vector<EdgeId> res;

		Sequence read = p_r.sequence();
		if (k_ + 1 <= read.size()) {
			runtime_k::RtSeq kmer = read.start<runtime_k::RtSeq::max_size>(k_ + 1);
			bool found;
			for (size_t i = k_ + 1; i <= read.size(); ++i) {
				if (index_.contains(kmer)) {
					pair<EdgeId, size_t> position = index_.get(kmer);
					found = false;
					for (size_t j = 0; j < res.size(); j++)
						if (res[j] == position.first) {
							found = true;
							break;
						}
					if (!found)
						res.push_back(position.first);
				}
				if (i != read.size())
					kmer <<= read[i];
			}
		}

		return res;
	}
};

template<class Graph>
class EdgeQuality: public GraphLabeler<Graph>, public GraphActionHandler<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	map<EdgeId, size_t> quality_;
	size_t k_;

public:

	void FillQuality(const EdgeIndex<Graph> &index
			, const KmerMapper<Graph>& kmer_mapper, const Sequence &genome) {
		if (genome.size() < k_)
			return;
		runtime_k::RtSeq cur = genome.start<runtime_k::RtSeq::max_size>(k_);
		cur >>= 0;
		for (size_t i = 0; i + k_ - 1 < genome.size(); i++) {
			cur <<= genome[i + k_ - 1];
			auto corr_cur = kmer_mapper.Substitute(cur);
			if (index.contains(corr_cur)) {
				quality_[index.get(corr_cur).first]++;
			}
		}
	}

	EdgeQuality(const Graph &graph, const EdgeIndex<Graph> &index,
	const KmerMapper<Graph>& kmer_mapper,
	const Sequence &genome) :

			GraphActionHandler<Graph>(graph, "EdgeQualityLabeler"),
			k_(kmer_mapper.get_k()) {
		FillQuality(index, kmer_mapper, genome);
		FillQuality(index, kmer_mapper, !genome);
	}

	virtual ~EdgeQuality() {
	}

	virtual void HandleAdd(EdgeId e) {
	}

	virtual void HandleDelete(EdgeId e) {
		quality_.erase(e);
	}

	virtual void HandleMerge(vector<EdgeId> old_edges, EdgeId new_edge) {
		size_t res = 0;
		for (size_t i = 0; i < old_edges.size(); i++) {
			res += quality_[old_edges[i]];
		}
		quality_[new_edge] += res;
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		quality_[new_edge] += quality_[edge2];
		quality_[new_edge] += quality_[edge1];
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge1,
			EdgeId new_edge2) {
		quality_[new_edge1] = quality_[old_edge] * this->g().length(new_edge1)
				/ (this->g().length(new_edge1) + this->g().length(new_edge2));
		quality_[new_edge2] = quality_[old_edge] * this->g().length(new_edge2)
				/ (this->g().length(new_edge1) + this->g().length(new_edge2));
	}

	double quality(EdgeId edge) const {
		auto it = quality_.find(edge);
		if (it == quality_.end())
			return 0.;
		else
			return 1. * it->second / this->g().length(edge);
	}

	bool IsPositiveQuality(EdgeId edge) const {
		return math::gr(quality(edge), 0.);
	}

	virtual std::string label(VertexId vertexId) const {
		return "";
	}

	virtual std::string label(EdgeId edge) const {
		double q = quality(edge);
		return (q == 0) ? "" : "quality: " + ToString(q);
	}

};

template<class Graph>
class QualityLoggingRemovalHandler {
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	const EdgeQuality<Graph>& quality_handler_;
//	size_t black_removed_;
//	size_t colored_removed_;
public:
	QualityLoggingRemovalHandler(const Graph& g, const EdgeQuality<Graph>& quality_handler) :
			g_(g), quality_handler_(quality_handler)/*, black_removed_(0), colored_removed_(
	 0)*/{

	}

	void HandleDelete(EdgeId edge) {
		if (math::gr(quality_handler_.quality(edge), 0.)) {
			TRACE("Deleting edge " << g_.int_id(edge) << " with quality " << quality_handler_.quality(edge));
		} else {
//			TRACE("Deleting edge " << g_.int_id(edge) << " with zero quality");
		}
//		if (math::gr(quality_handler_.quality(edge), 0.))
//			colored_removed_++;
//		else
//			black_removed_++;
	}

private:
	DECL_LOGGER("QualityLoggingRemovalHandler")
	;
};

template<class Graph>
class QualityLoggingRemovalCountHandler {
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	const EdgeQuality<Graph>& quality_handler_;
	size_t black_removed_;
    size_t total;

public:
	QualityLoggingRemovalCountHandler(const Graph& g, const EdgeQuality<Graph>& quality_handler) :
			g_(g), quality_handler_(quality_handler)/*, black_removed_(0), colored_removed_(
	 0)*/{
        black_removed_ = 0;
        total = 0;
	}

	void HandleDelete(EdgeId edge) {
        total++;
		if (math::gr(quality_handler_.quality(edge), 0.)) {
            TRACE("Deleting good edge " << g_.int_id(edge) << " with quality " << quality_handler_.quality(edge) << " cov " << g_.coverage(edge) << " length " << g_.length(edge));
		}else{
            black_removed_++;
        }
        if ((total % (1<<10)) != 0)   
            TRACE("Removed still " << black_removed_ << " " << total); 
	}

private:
};

template<class Graph>
class EdgeNeighborhoodFinder: public omnigraph::GraphSplitter<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	EdgeId edge_;
	size_t max_size_;
	size_t edge_length_bound_;
	bool finished_;
public:
	EdgeNeighborhoodFinder(const Graph &graph, EdgeId edge, size_t max_size
			, size_t edge_length_bound) :
			GraphSplitter<Graph>(graph), edge_(edge), max_size_(
					max_size), edge_length_bound_(edge_length_bound), finished_(
					false) {
	}

	/*virtual*/ vector<VertexId> NextComponent() {
		CountingDijkstra<Graph> cf(this->graph(), max_size_,
				edge_length_bound_);
		set<VertexId> result_set;
		cf.run(this->graph().EdgeStart(edge_));
		vector<VertexId> result_start = cf.ReachedVertices();
		result_set.insert(result_start.begin(), result_start.end());
		cf.run(this->graph().EdgeEnd(edge_));
		vector<VertexId> result_end = cf.ReachedVertices();
		result_set.insert(result_end.begin(), result_end.end());

		ComponentCloser<Graph> cc(this->graph(), edge_length_bound_);
		cc.CloseComponent(result_set);

		finished_ = true;
		vector<VertexId> result;
		for (auto it = result_set.begin(); it != result_set.end(); ++it)
			result.push_back(*it);
		return result;
	}

	/*virtual*/ bool Finished() {
		return finished_;
	}
};

void EmptyHandleF(EdgeId edge) {

}

template<class Graph>
class QualityEdgeLocalityPrintingRH {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	const Graph& g_;
	const EdgeQuality<Graph>& quality_handler_;
	const GraphLabeler<Graph>& labeler_;
	const string& output_folder_;
//	size_t black_removed_;
//	size_t colored_removed_;
public:
	QualityEdgeLocalityPrintingRH(const Graph& g
			, const EdgeQuality<Graph>& quality_handler
			, const GraphLabeler<Graph>& labeler
			, const string& output_folder) :
			g_(g), quality_handler_(quality_handler),
			labeler_(labeler), output_folder_(output_folder){
	}

	void HandleDelete(EdgeId edge) {
        if (quality_handler_.IsPositiveQuality(edge)) {
			DEBUG("Deleting edge " << g_.str(edge) << " with quality " << quality_handler_.quality(edge));
			string folder = output_folder_ + "colored_edges_deleted/";
            path::make_dir(folder);
			//todo magic constant
//			map<EdgeId, string> empty_coloring;
			EdgeNeighborhoodFinder<Graph> splitter(g_, edge, 50,
					250);
			WriteComponents(g_, splitter/*, "locality_of_edge_" + ToString(g_.int_id(edge))*/
					, folder + "edge_" +  ToString(g_.int_id(edge)) + "_" + ToString(quality_handler_.quality(edge)) + ".dot"
					, *DefaultColorer(g_), labeler_);
		} else {
			TRACE("Deleting edge " << g_.str(edge) << " with quality " << quality_handler_.quality(edge));
		}
	}

private:
	DECL_LOGGER("QualityEdgeLocalityPrintingRH")
	;
};

template<class Graph>
class EdgePairInfoWrapper {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef omnigraph::PairInfo<EdgeId> PairInfo;
	typedef vector<PairInfo> PairInfos;
	const Graph& g_;
    const PairedInfoIndex<ConjugateDeBruijnGraph>& index_;
//	size_t black_removed_;
//	size_t colored_removed_;
public:
	EdgePairInfoWrapper(const Graph& g
            , const PairedInfoIndex<ConjugateDeBruijnGraph>& index) :
			g_(g), index_(index) {
	}

	double GetTotalWeight(EdgeId edge) {
            PairInfos infos = index_.GetEdgeInfo(edge);
            double total_weight = 0.;
            for (size_t i = 0; i<infos.size(); i++){
                total_weight += infos[i].weight;
                INFO("Tip Info " << infos[i].first << " " << infos[i].second << " " << infos[i].d << " " << infos[i].weight << " " << infos[i].variance);
            }
            return total_weight;
	}

private:
};

template<class Graph>
class QualityPairInfoHandler {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef omnigraph::PairInfo<EdgeId> PairInfo;
	typedef vector<PairInfo> PairInfos;
	const Graph& g_;
	const EdgeQuality<Graph>& quality_handler_;
	const GraphLabeler<Graph>& labeler_;
	const string& output_folder_;
    const PairedInfoIndex<ConjugateDeBruijnGraph>& index_;
//	size_t black_removed_;
//	size_t colored_removed_;
public:
	QualityPairInfoHandler(const Graph& g
			, const EdgeQuality<Graph>& quality_handler
			, const GraphLabeler<Graph>& labeler
			, const string& output_folder
            , const PairedInfoIndex<ConjugateDeBruijnGraph>& index) :
			g_(g), quality_handler_(quality_handler),
			labeler_(labeler), output_folder_(output_folder), index_(index) {
	}

	void HandleDelete(EdgeId edge) {
        if (quality_handler_.IsPositiveQuality(edge)) {
            cout << "Deleting edge " << g_.str(edge) << " with quality " << quality_handler_.quality(edge) << endl;
            string folder = output_folder_ + "colored_edges_deleted/";
            path::make_dir(folder);
            //todo magic constant
            PairInfos infos = index_.GetEdgeInfo(edge);
            if (infos.size() > 0){
                for (size_t i = 0; i<infos.size(); i++){
                    cout << "Tip Info " << g_.int_id(infos[i].first) << " " << g_.int_id(infos[i].second) << " " << infos[i].d << " " << infos[i].weight << " " << infos[i].variance << endl;
                }
            }
            map<EdgeId, string> empty_coloring;
            EdgeNeighborhoodFinder<Graph> splitter(g_, edge, 50,
                    250);
            
            WriteComponents(g_, splitter, TrueFilter<vector<VertexId>>(), "locality_of_edge_" + ToString(g_.int_id(edge))
                    , folder + "edge_" +  ToString(g_.int_id(edge)) + "_" + ToString(quality_handler_.quality(edge)) + ".dot"
                    , empty_coloring, labeler_);
        }
	}

private:
};


template<class Graph>
class EdgeLocalityPrintingRH {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	const Graph& g_;
	const GraphLabeler<Graph>& labeler_;
	const string& output_folder_;
    boost::function<double (EdgeId)>& quality_f_;
//	size_t black_removed_;
//	size_t colored_removed_;
public:
	EdgeLocalityPrintingRH(const Graph& g
			, const GraphLabeler<Graph>& labeler
			, const string& output_folder
            , boost::function<double (EdgeId)> quality_f = 0) :
			g_(g), 
			labeler_(labeler), output_folder_(output_folder),
            quality_f_(quality_f){
	}

	void HandleDelete(EdgeId edge) {
            TRACE("Deleting edge " << g_.str(edge));
            if (quality_f_ && math::gr(quality_f_(edge), 0.)) 
                INFO("EdgeLocalityPrintRH handling the edge with positive quality : " << quality_f_(edge) << " " << g_.str(edge));
        
            string folder = output_folder_ + "edges_deleted/";
            path::make_dir(folder);
            //todo magic constant
            map<EdgeId, string> empty_coloring;
            EdgeNeighborhoodFinder<Graph> splitter(g_, edge, 50,
                    250);

            WriteComponents(g_, splitter, TrueFilter<vector<VertexId>>(), "locality_of_edge_" + ToString(g_.int_id(edge))
                    , folder + "edge_" +  ToString(g_.int_id(edge)) + ".dot", empty_coloring, labeler_);
	}

private:
	DECL_LOGGER("QualityEdgeLocalityPrintingRH")
	;
};

class WeightDEWrapper {
private:

    vector<double> new_hist;
    int left_x;
    size_t insert_size;

	void ExtendLinear(const std::map<int, size_t> & hist) {
        size_t sum_weight = 0;

        for (auto iter = hist.begin(); iter != hist.end(); ++iter)
            sum_weight += iter->second;
            

        auto iter = hist.begin();

        left_x = iter->first;

        int prev = iter->first;
        int prev_val = iter->second;

        new_hist.push_back(prev_val * 1./ sum_weight);
        ++iter;

        for (; iter != hist.end(); ++iter) {
            int x = iter->first;
            int y = iter->second;
            double tan = 1. * (y - prev_val) / (x - prev);

            for (int i = prev + 1; i <= x; ++i) {
                new_hist.push_back((prev_val + tan * (i - prev)) * 1. / sum_weight);
            }
            prev = x;
            prev_val = y;
        }
	}

public:
    WeightDEWrapper(const map<int, size_t>& hist, double IS) {
        insert_size = (size_t) IS;
        ExtendLinear(hist);
    }

    ~WeightDEWrapper() {
    }


    double CountWeight(int x) const {
        int xx = insert_size - left_x + x - 1;
        
        if (!(xx >= 0 && xx < (int) new_hist.size())) return 0.;
        //cout << x << " " << xx << " " << new_hist[xx] << endl;
        VERIFY(math::le(new_hist[xx], 1.));
        return 1000. * new_hist[xx];
    }
};



double UnityFunction(int x) {
    return 1.;   
}

void DivideClusters(const Graph& g, PairedInfoIndex<Graph>& clustered_index, const GraphDistanceFinder<Graph>& dist_finder) {

    vector<PairInfo<EdgeId> > vector_to_remove;
    vector<PairInfo<EdgeId> > vector_to_add;
    for (auto iter = clustered_index.begin(); iter != clustered_index.end(); ++iter) {
        const vector<PairInfo<EdgeId> >& data = *iter;
        EdgeId first_edge = data[0].first;
        EdgeId second_edge = data[0].second;
        DEBUG("Analyzing edges " << g.int_id(first_edge) << " and " << g.int_id(second_edge));

        const vector<vector<EdgeId> >& paths = dist_finder.GetGraphDistances(first_edge, second_edge);
        //count the lengths of corresponding paths
        DEBUG(paths.size() << " paths found");
        vector<size_t> paths_lengths;
        for (size_t i = 0; i < paths.size(); ++i) {
            const vector<EdgeId>& path = paths[i];
            size_t len = g.length(first_edge);
            for (size_t j = 0; j + 1 < path.size(); ++j)
                len += g.length(path[j]);
            paths_lengths.push_back(len);
        }
        DEBUG("Lengths counted");
        
        for (size_t iter = 0; iter < data.size(); ++iter) {
            const PairInfo<EdgeId>& pair_info = data[iter];
            DEBUG("New pair info " << pair_info);
            // choose only clusters
            if (math::gr(pair_info.variance, 0.)) {
                // filtering paths with corresponding length
                DEBUG("Variance > 0");
                vector_to_remove.push_back(pair_info);
                double average_weight = 0.;
                vector<double> paths_weights;
                size_t num_of_clustered_paths = 0;
                for (size_t i = 0; i < paths.size(); ++i) {
                    double weight_total = 0.;
                    if (math::le(abs(paths_lengths[i] - pair_info.d), pair_info.variance)) {
                        ++num_of_clustered_paths;
                        const vector<EdgeId>& path = paths[i];
                        double cur_len = g.length(first_edge);
                        for (size_t j = 0; j + 1 < path.size(); ++j) {
                            const vector<PairInfo<EdgeId> >& back_infos = clustered_index.GetEdgePairInfo(first_edge, path[j]);
                            for (auto iter = back_infos.begin(); iter != back_infos.end(); ++iter) {
                                const PairInfo<EdgeId>& info = *iter;
                                if (info.d == cur_len) {
                                    weight_total += info.weight;   
                                }
                            }
                            const vector<PairInfo<EdgeId> >& forward_infos = clustered_index.GetEdgePairInfo(path[j], second_edge);
                            for (auto iter = forward_infos.begin(); iter != forward_infos.end(); ++iter) {
                                const PairInfo<EdgeId>& info = *iter;
                                if (info.d + cur_len == paths_lengths[i]) {
                                    weight_total += info.weight;   
                                }
                            }
                        }
                    }
                    paths_weights.push_back(weight_total);
                    average_weight += weight_total;
                }
                double sum_weight = average_weight;
                average_weight /= 1. * num_of_clustered_paths;
                // filtering bad paths
                for (size_t i = 0; i < paths.size(); ++i) {
                    if (math::le(abs(paths_lengths[i] - pair_info.d), pair_info.variance)) {
                        if (math::gr(paths_weights[i], average_weight)) {
                            PairInfo<EdgeId> new_pair_info(first_edge, second_edge, paths_lengths[i], 
                                                            pair_info.weight / num_of_clustered_paths * paths_weights[i] / sum_weight, 0.);
                            DEBUG("Adding new pair info " << new_pair_info);
                            vector_to_add.push_back(new_pair_info);
                        }
                    }
                }
            } else {
                DEBUG("Variance zero");   
            }
        }
    }
    for (size_t i = 0; i < vector_to_remove.size(); ++i)
        clustered_index.RemovePairInfo(vector_to_remove[i]);

    for (size_t i = 0; i < vector_to_add.size(); ++i)
        clustered_index.AddPairInfo(vector_to_add[i], false);
       
}


}
