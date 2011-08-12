/*
 * utils.hpp
 *
 *  Created on: Apr 5, 2011
 *      Author: sergey
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "io/paired_read.hpp"
#include "seq_map.hpp"
#include "omni_utils.hpp"
#include "logging.hpp"
#include "paired_info.hpp"
#include "statistics.hpp"
#include "xmath.h"
#include <boost/optional.hpp>
//#include "common/io/paired_read.hpp"
namespace debruijn_graph {

using omnigraph::Path;
using omnigraph::MappingPath;
using omnigraph::Range;
using omnigraph::MappingRange;
using omnigraph::PairInfo;
using omnigraph::GraphActionHandler;
//using io::PairedRead;

template<size_t kmer_size_, typename Graph>
class ReadThreaderResult {
	typedef typename Graph::EdgeId EdgeId;

	Path<EdgeId> left_read_, right_read_;
	int gap_;
public:
	ReadThreaderResult(Path<EdgeId> left_read, Path<EdgeId> right_read, int gap) :
			gap_(gap), left_read_(left_read), right_read_(right_read) {
	}
};

template< typename Graph>
class SingleReadThreaderResult {
	typedef typename Graph::EdgeId EdgeId;
public:
	EdgeId edge_;
	int read_position_;
	int edge_position_;
	SingleReadThreaderResult(EdgeId edge, int read_position, int edge_position): edge_(edge), read_position_(read_position), edge_position_(edge_position){
	}
};

template< typename Graph>
class ReadMappingResult{
public:
	Sequence read_;
	vector<SingleReadThreaderResult<Graph> > res_;
	ReadMappingResult(Sequence read, vector<SingleReadThreaderResult<Graph> > res): read_(read), res_(res){

	}
	ReadMappingResult(){

	}
};
/**
 * DataHashRenewer listens to add/delete events and updates index according to those events. This class
 * can be used both with vertices and edges of graph.
 * todo EdgeNucls are hardcoded!
 */
template<size_t kmer_size_, typename Graph, typename ElementId>
class DataHashRenewer {

	typedef Seq<kmer_size_> Kmer;
	typedef SeqMap<kmer_size_, ElementId> Index;
	const Graph &g_;

	Index &index_;

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewKmersHash(ElementId id) {
		Sequence nucls = g_.EdgeNucls(id);
//		DEBUG("Renewing hashes for k-mers of sequence " << nucls);
		index_.RenewKmersHash(nucls, id);
	}

	void DeleteKmersHash(ElementId id) {
		Sequence nucls = g_.EdgeNucls(id);
//		DEBUG("Deleting hashes for k-mers of sequence " << nucls);
		index_.DeleteKmersHash(nucls, id);
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
template<size_t k, class Graph>
class EdgeIndex: public GraphActionHandler<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef SeqMap<k, EdgeId> InnerIndex;
	typedef Seq<k> Kmer;
	InnerIndex inner_index_;
	DataHashRenewer<k, Graph, EdgeId> renewer_;
	bool delete_index_;
public:

	EdgeIndex(const Graph& g) :
			GraphActionHandler<Graph>(g, "EdgeIndex"), inner_index_(), renewer_(
					g, inner_index_), delete_index_(true) {
	}

	virtual ~EdgeIndex() {
		TRACE("~EdgeIndex OK")
	}

	InnerIndex &inner_index() {
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

	bool containsInIndex(const Kmer& kmer) const {
		return inner_index_.containsInIndex(kmer);
	}

	const pair<EdgeId, size_t>& get(const Kmer& kmer) const {
		return inner_index_.get(kmer);
	}

};

template<size_t k, class Graph>
class KmerMapper : public omnigraph::GraphActionHandler<Graph> {
	typedef omnigraph::GraphActionHandler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef Seq<k> Kmer;
	typedef typename std::tr1::unordered_map<Kmer, Kmer, typename Kmer::hash> MapType;

	void RemapKmers(const Sequence& old_s, const Sequence& new_s) {
		Kmer old_kmer = old_s.start<k>() >> 0;
		for (size_t i = k - 1; i < old_s.size(); ++i) {
			old_kmer << old_s[i];
			size_t old_kmer_offset = i - k + 1;
			size_t new_kmer_offest = std::floor(1. * old_kmer_offset / (old_s.size() - k + 1) * (new_s.size() - k + 1) + 0.5);
			Kmer new_kmer(new_s, new_kmer_offest);
			mapping_[old_kmer] = new_kmer;
		}
	}

	MapType mapping_;

public:
	KmerMapper(const Graph& g) : base(g, "KmerMapper") {

	}

	virtual ~KmerMapper() {

	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		assert(this->g().GetNucls(new_edge) == this->g().GetNucls(edge2));
		RemapKmers(this->g().GetNucls(edge1), this->g().GetNucls(edge2));
	}

	Kmer Substitute(const Kmer& kmer) {
		Kmer answer = kmer;
		auto it = mapping_.find(answer);
		while (it != mapping_.end()) {
			answer = (*it).second;
			it = mapping_.find(answer);
		}
		return answer;
	}
};

/**
 * This class finds how certain sequence is mapped to genome. As it is now it works correct only if sequence
 * is mapped to graph ideally and in unique way.
 */
template<size_t k, class Graph>
class SimpleSequenceMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef Seq<k> Kmer;
	typedef EdgeIndex<k, Graph> Index;
private:
	const Graph& g_;
	const Index &index_;

	bool TryThread(Kmer &kmer, vector<EdgeId> &passed,
			size_t &endPosition) const {
		EdgeId last = passed[passed.size() - 1];
		if (endPosition + 1 < g_.length(last)) {
			if (g_.EdgeNucls(last)[endPosition + k] == kmer[k - 1]) {
				endPosition++;
				return true;
			}
		} else {
			vector<EdgeId> edges = g_.OutgoingEdges(g_.EdgeEnd(last));
			for (size_t i = 0; i < edges.size(); i++) {
				if (g_.EdgeNucls(edges[i])[k - 1] == kmer[k - 1]) {
					passed.push_back(edges[i]);
					endPosition = 0;
					return true;
				}
			}
		}
		return false;
	}

	bool FindKmer(Kmer kmer, vector<EdgeId> &passed,
			size_t &startPosition, size_t &endPosition) const {
		if (index_.containsInIndex(kmer)) {
			pair<EdgeId, size_t> position = index_.get(kmer);
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

	bool ProcessKmer(Kmer &kmer, vector<EdgeId> &passed,
			size_t &startPosition, size_t &endPosition, bool valid) const {
		if (valid) {
			return TryThread(kmer, passed, endPosition);
		} else {
			return FindKmer(kmer, passed, startPosition, endPosition);
		}
	}

public:
	/**
	 * Creates SimpleSequenceMapper for given graph. Also requires index_ which should be synchronized
	 * with graph.
	 * @param g graph sequences should be mapped to
	 * @param index index synchronized with graph
	 */
	SimpleSequenceMapper(const Graph& g, const Index& index) :
			g_(g), index_(index) {
	}

	/**
	 * Finds a path in graph which corresponds to given sequence.
	 * @read sequence to be mapped
	 */

	Path<EdgeId> MapSequence(const Sequence &read) const {
		vector<EdgeId> passed;
		if (read.size() <= k - 1) {
			return Path<EdgeId>();
		}
		Kmer kmer = read.start<k>();
		size_t startPosition = -1;
		size_t endPosition = -1;
		bool valid = ProcessKmer(kmer, passed, startPosition, endPosition,
				false);
		for (size_t i = k; i < read.size(); ++i) {
			kmer = kmer << read[i];
			valid = ProcessKmer(kmer, passed, startPosition, endPosition,
					valid);
		}
		return Path<EdgeId>(passed, startPosition, endPosition + 1);
	}

};

//todo optimize if needed
template<size_t k, class Graph>
class ExtendedSequenceMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef vector<MappingRange> RangeMappings;
	typedef Seq<k> Kmer;
	typedef EdgeIndex<k, Graph> Index;
	typedef KmerMapper<k, Graph> KmerSubs;

private:
	const Graph& g_;
	const Index& index_;
	const KmerSubs& kmer_mapper_;
	bool glue_across_gap_;

	void FindKmer(Kmer kmer, size_t kmer_pos, vector<EdgeId> &passed,
			RangeMappings& range_mappings) const {

		if (index_.containsInIndex(kmer)) {
			pair<EdgeId, size_t> position = index_.get(kmer);
			if (passed.empty() || passed.back() != position.first
					|| kmer_pos != range_mappings.back().initial_range.end_pos
					|| position.second + 1 < range_mappings.back().mapped_range.end_pos) {
				passed.push_back(position.first);
				MappingRange mapping_range(Range(kmer_pos, kmer_pos + 1), Range(position.second, position.second + 1));
				range_mappings.push_back(mapping_range);
			} else {
				range_mappings.back().initial_range.end_pos = kmer_pos + 1;
				range_mappings.back().mapped_range.end_pos = position.second + 1;
			}
		}
	}

	void ProcessKmer(Kmer kmer, size_t kmer_pos, vector<EdgeId> &passed, RangeMappings& interval_mapping) const {
		kmer = kmer_mapper_.Substitute(kmer);
		FindKmer(kmer, kmer_pos, passed, interval_mapping);
	}

public:
	ExtendedSequenceMapper(const Graph& g, const Index& index, const KmerSubs& kmer_mapper) :
			g_(g), index_(index), kmer_mapper_(kmer_mapper) {
	}

	MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
		vector<EdgeId> passed_edges;
		RangeMappings range_mapping;

		assert(sequence.size() >= k);
		Kmer kmer = sequence.start<k>() >> 0;
		for (size_t i = k - 1; i < sequence.size(); ++i) {
			kmer = kmer << sequence[i];
			ProcessKmer(kmer, i - k + 1, passed_edges, range_mapping);
		}
		return MappingPath<EdgeId>(passed_edges, range_mapping);
	}
};

template<size_t k, class Graph>
class EtalonPairedInfoCounter {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
	const EdgeIndex<k + 1, Graph>& index_;
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
		SimpleSequenceMapper<k + 1, Graph> sequence_mapper(g_, index_);
		Path<EdgeId> path = sequence_mapper.MapSequence(sequence);

		for (size_t i = 0; i < path.size(); ++i) {
			EdgeId e = path[i];
			if (g_.length(e) + delta_ > gap_ + k + 1) {
				AddEtalonInfo(paired_info, e, e, 0);
			}
			size_t j = i + 1;
			size_t length = 0;

			while (j < path.size() && length <= omnigraph::PairInfoPathLengthUpperBound(k, insert_size_, delta_)) {
				if (length >= omnigraph::PairInfoPathLengthLowerBound(k, g_.length(e), g_.length(path[j]), gap_, delta_)) {
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

	EtalonPairedInfoCounter(const Graph& g, const EdgeIndex<k + 1, Graph>& index
			, size_t insert_size, size_t read_length, size_t delta) :
			g_(g), index_(index), insert_size_(insert_size), read_length_(
					read_length), gap_(insert_size_ - 2 * read_length_), delta_(
					delta) {
		assert(insert_size_ >= 2 * read_length_);
	}

	void FillEtalonPairedInfo(const Sequence& genome,
			omnigraph::PairedInfoIndex<Graph>& paired_info) {
		ProcessSequence(genome, paired_info);
		ProcessSequence(!genome, paired_info);
		//DEBUG
//		CheckPairInfo(genome, paired_info);
	}
};

template<size_t k, class Graph>
class NewEtalonPairedInfoCounter {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
	const EdgeIndex<k + 1, Graph>& index_;
	size_t insert_size_;
	size_t read_length_;
	size_t gap_;
	size_t delta_;

	void AddEtalonInfo(set<PairInfo<EdgeId>>& paired_info, EdgeId e1, EdgeId e2,
			double d) {
		PairInfo<EdgeId> pair_info(e1, e2, d, 1000.0, 0.);
		paired_info.insert(pair_info);
	}

	void ProcessSequence(const Sequence& sequence,
			set<PairInfo<EdgeId>>& temporary_info) {
		int mod_gap = (gap_ > delta_) ? gap_ - delta_ : 0;
		Seq<k + 1> left(sequence);
		left = left >> 0;
		for (size_t left_idx = 0; left_idx + k + 1 + mod_gap <= sequence.size(); ++left_idx) {
			left = left << sequence[left_idx + k];
			if (!index_.containsInIndex(left)) {
				continue;
			}
			pair<EdgeId, size_t> left_pos = index_.get(left);

			size_t right_idx = left_idx + mod_gap;
			Seq<k + 1> right(sequence, right_idx);
			right = right >> 0;
			for (; right_idx + k + 1 <= left_idx + insert_size_ + delta_ && right_idx + k + 1 <= sequence.size(); ++right_idx) {
				right = right << sequence[right_idx + k];
				if (!index_.containsInIndex(right)) {
					continue;
				}
				pair<EdgeId, size_t> right_pos = index_.get(right);

				AddEtalonInfo(temporary_info, left_pos.first, right_pos.first, right_idx - left_idx + left_pos.second - right_pos.second);
			}
		}
	}

public:

	NewEtalonPairedInfoCounter(const Graph& g, const EdgeIndex<k + 1, Graph>& index
			, size_t insert_size, size_t read_length, size_t delta) :
			g_(g), index_(index), insert_size_(insert_size), read_length_(
					read_length), gap_(insert_size_ - 2 * read_length_), delta_(
					delta) {
		assert(insert_size_ >= 2 * read_length_);
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

template<class Graph, size_t k>
class GenomeMappingStat: public omnigraph::AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	const EdgeIndex<k + 1, Graph>& index_;
	Sequence genome_;
public:
	GenomeMappingStat(Graph &graph, const EdgeIndex<k + 1, Graph> &index,
	Sequence genome) :
			graph_(graph), index_(index), genome_(genome) {
	}

	virtual ~GenomeMappingStat() {
	}

	virtual void Count() {
		INFO("Mapping genome");
		size_t break_number = 0;
		size_t covered_kp1mers = 0;
		size_t fail = 0;
		Seq<k + 1> cur = genome_.start<k + 1>() >> 0;
		bool breaked = true;
		pair<EdgeId, size_t> cur_position;
		for (size_t cur_nucl = k; cur_nucl < genome_.size(); cur_nucl++) {
			cur = cur << genome_[cur_nucl];
			if (index_.containsInIndex(cur)) {
				pair<EdgeId, size_t> next = index_.get(cur);
				if (!breaked
						&& cur_position.second + 1
								< graph_.length(cur_position.first)) {
					if (next.first != cur_position.first
							|| cur_position.second + 1 != next.second) {
						fail++;
					}
				}
				cur_position = next;
				covered_kp1mers++;
				breaked = false;
			} else {
				if (!breaked) {
					breaked = true;
					break_number++;
				}
			}
		}INFO("Genome mapped");
		INFO("Genome mapping results:");
		INFO(
				"Covered k+1-mers:" << covered_kp1mers << " of " << (genome_.size() - k) << " which is " << (100.0 * covered_kp1mers / (genome_.size() - k)) << "%");
		INFO(
				"Covered k+1-mers form " << break_number + 1 << " contigious parts");
		INFO("Continuity failtures " << fail);
	}
};

template<class Graph, size_t k>
class StatCounter: public omnigraph::AbstractStatCounter {
private:
	omnigraph::StatList stats_;
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	StatCounter(Graph& graph, const EdgeIndex<k + 1, Graph>& index,
	const Sequence& genome) {
		SimpleSequenceMapper<k + 1, Graph> sequence_mapper(graph, index);
		Path<EdgeId> path1 = sequence_mapper.MapSequence(Sequence(genome));
		Path<EdgeId> path2 = sequence_mapper.MapSequence(!Sequence(genome));
		stats_.AddStat(new omnigraph::VertexEdgeStat<Graph>(graph));
		stats_.AddStat(
				new omnigraph::BlackEdgesStat<Graph>(graph, path1, path2));
		stats_.AddStat(new omnigraph::NStat<Graph>(graph, path1, 50));
		stats_.AddStat(new omnigraph::SelfComplementStat<Graph>(graph));
		stats_.AddStat(
				new GenomeMappingStat<Graph, k>(graph, index,
						Sequence(genome)));
	}

	virtual ~StatCounter() {
		stats_.DeleteStats();
	}

	virtual void Count() {
		stats_.Count();
	}

private:
	DECL_LOGGER("StatCounter")
};

/**
 * As for now it ignores sophisticated case of repeated consecutive
 * occurrence of edge in path due to gaps in mapping
 *
 * todo talk with Anton about simplification and speed-up of procedure with little quality loss
 */
template<size_t k, class Graph, class Stream>
class ReadCountPairedIndexFiller {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef Seq<k> Kmer;
	Graph &graph_;
	const ExtendedSequenceMapper<k, Graph>& mapper_;
	Stream& stream_;

	inline size_t CountDistance(const io::PairedRead& paired_read) {
		return paired_read.distance() - paired_read.second().size();
	}

	void ProcessPairedRead(
			omnigraph::PairedInfoIndex<Graph> &paired_index,
			const io::PairedRead& p_r) {
		Sequence read1 = p_r.first().sequence();
		Sequence read2 = p_r.second().sequence();

		MappingPath<EdgeId> path1 = mapper_.MapSequence(read1);
		MappingPath<EdgeId> path2 = mapper_.MapSequence(read2);
		size_t read_distance = CountDistance(p_r);
		for (size_t i = 0; i < path1.size(); ++i) {
			pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
			for (size_t j = 0; j < path2.size(); ++j) {
				pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
				double weight = 1;
				size_t kmer_distance = read_distance + mapping_edge_2.second.initial_range.start_pos - mapping_edge_1.second.initial_range.start_pos;
				size_t edge_distance = kmer_distance + mapping_edge_1.second.mapped_range.start_pos - mapping_edge_2.second.mapped_range.start_pos;
				PairInfo<EdgeId> new_info(mapping_edge_1.first, mapping_edge_2.second, edge_distance, weight);
				paired_index.AddPairInfo(new_info);
			}
		}
	}

public:

	ReadCountPairedIndexFiller(const Graph &graph, const ExtendedSequenceMapper<k, Graph>& mapper, Stream& stream) :
			graph_(graph), mapper_(mapper), stream_(stream) {

	}

	void FillIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			paired_index.AddPairInfo(PairInfo<EdgeId>(*it, *it, 0, 0.0));
		}
		stream_.reset();
		while (!stream_.eof()) {
			io::PairedRead p_r;
			stream_ >> p_r;
			ProcessPairedRead(paired_index, p_r);
		}
	}

};


/**
 * This class finds how certain _paired_ read is mapped to genome. As it is now it is hoped to work correctly only if read
 * is mapped to graph ideally and in unique way.
 */
template<size_t k, class Graph, class Stream>
class TemplateReadMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef EdgeIndex<k + 1, Graph> Index;
private:
	SimpleSequenceMapper<k, Graph> read_seq_mapper;
	Stream& stream_;
public:
	/**
	 * Creates TemplateReadMapper for given graph. Also requires index_ which should be synchronized
	 * with graph.
	 * @param g graph sequences should be mapped to
	 * @param index index syncronized with graph
	 */
	TemplateReadMapper(const Graph& g, const Index& index, Stream & stream):
		read_seq_mapper(g, index), stream_(stream) {
		stream_.reset();
	}

	ReadThreaderResult<k + 1, Graph> ThreadNext() {
		if (!stream_.eof()) {
      io::PairedRead p_r;
			stream_ >> p_r;
			Sequence read1 = p_r.first().sequence();
			Sequence read2 = p_r.second().sequence();
			Path<EdgeId> aligned_read[2];
			aligned_read[0] = read_seq_mapper.MapSequence(read1);
			aligned_read[1] = read_seq_mapper.MapSequence(read2);
			size_t distance = p_r.distance();
			int current_distance1 = distance + aligned_read[0].start_pos()
					- aligned_read[1].start_pos();
			return ReadThreaderResult<k + 1, Graph>(aligned_read[0], aligned_read[1], current_distance1);
		}
//		else return NULL;
	}

};

template<size_t k, class Graph>
class SingleReadMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef EdgeIndex<k + 1, Graph> Index;
private:
	SimpleSequenceMapper<k + 1, Graph> read_seq_mapper;
	const Graph& g_;
	const Index& index_;
public:
	/**
	 * Creates SingleReadMapper for given graph. Also requires index_ which should be synchronized
	 * with graph.
	 * @param g graph sequences should be mapped to
	 * @param index index syncronized with graph
	 */
	SingleReadMapper(const Graph& g, const Index& index):
		read_seq_mapper(g, index),  g_(g), index_(index) {
	}

	vector<EdgeId> GetContainingEdges(io::SingleRead& p_r){
		vector<EdgeId> res;

		Sequence read = p_r.sequence();
		if (k+1 <= read.size()) {
			Seq<k + 1> kmer = read.start<k + 1>();
			bool found;
			for (size_t i = k + 1; i <= read.size(); ++i) {
				if (index_.containsInIndex(kmer)) {
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
					kmer = kmer << read[i];
			}
		}

		return res;
	}

};

}

#endif /* UTILS_HPP_ */
