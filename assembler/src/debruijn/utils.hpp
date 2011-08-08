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
//#include "common/io/paired_read.hpp"
namespace debruijn_graph {

using omnigraph::Path;
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
	Graph& g_;
	InnerIndex inner_index_;
	DataHashRenewer<k, Graph, EdgeId> renewer_;
	bool delete_index_;
public:

	EdgeIndex(Graph& g) :
			GraphActionHandler<Graph>("EdgeIndex"), g_(g), inner_index_(), renewer_(
					g, inner_index_), delete_index_(true) {
		g_.AddActionHandler(this);
	}

	virtual ~EdgeIndex() {
		TRACE("~EdgeIndex")
		g_.RemoveActionHandler(this);
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

template<size_t kmer_size_, typename Graph>
class VertexHashRenewer: public GraphActionHandler<Graph> {

	typedef typename Graph::VertexId VertexId;

	DataHashRenewer<kmer_size_, Graph, VertexId> renewer_;

public:
	VertexHashRenewer(const Graph& g, SeqMap<kmer_size_, VertexId> *index) :
			renewer_(g, index) {
	}

	virtual void HandleAdd(VertexId e) {
		renewer_.HandleAdd(e);
	}

	virtual void HandleDelete(VertexId e) {
		renewer_.HandleDelete(e);
	}
};

class NoInfo {
};

/**
 * This class finds how certain sequence is mapped to genome. As it is now it works correct only if sequence
 * is mapped to graph ideally and in unique way.
 */
template<size_t k, class Graph>
class SimpleSequenceMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef EdgeIndex<k + 1, Graph> Index;
private:
	const Graph& g_;
	const Index &index_;

	bool TryThread(Seq<k + 1> &kmer, vector<EdgeId> &passed,
			size_t &endPosition) const {
		EdgeId last = passed[passed.size() - 1];
		if (endPosition + 1 < g_.length(last)) {
			if (g_.EdgeNucls(last)[endPosition + k + 1] == kmer[k]) {
				endPosition++;
				return true;
			}
		} else {
			vector<EdgeId> edges = g_.OutgoingEdges(g_.EdgeEnd(last));
			for (size_t i = 0; i < edges.size(); i++) {
				if (g_.EdgeNucls(edges[i])[k] == kmer[k]) {
					passed.push_back(edges[i]);
					endPosition = 0;
					return true;
				}
			}
		}
		return false;
	}

	bool FindKmer(Seq<k + 1> &kmer, vector<EdgeId> &passed,
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

	bool ProcessKmer(Seq<k + 1> &kmer, vector<EdgeId> &passed,
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
	 * @param index index syncronized with graph
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
		if (read.size() <= k) {
			return Path<EdgeId>();
		}
		Seq<k + 1> kmer = read.start<k + 1>();
		size_t startPosition = -1;
		size_t endPosition = -1;
		bool valid = ProcessKmer(kmer, passed, startPosition, endPosition,
				false);
		for (size_t i = k + 1; i < read.size(); ++i) {
			kmer = kmer << read[i];
			valid = ProcessKmer(kmer, passed, startPosition, endPosition,
					valid);
		}
		return Path<EdgeId>(passed, startPosition, endPosition + 1);
	}

};

template<size_t k, class Graph>
class EtalonPairedInfoCounter {
	typedef typename Graph::EdgeId EdgeId;

	Graph& g_;
	const EdgeIndex<k + 1, Graph>& index_;
	size_t insert_size_;
	size_t read_length_;
	size_t gap_;
	size_t delta_;

	void AddEtalonInfo(omnigraph::PairedInfoIndex<Graph>& paired_info,
			EdgeId e1, EdgeId e2, double d) {
		PairInfo<EdgeId> pair_info(e1, e2, d, 1000.0);
		paired_info.AddPairInfo(pair_info);
	}

	void ProcessSequence(const Sequence& sequence,
			omnigraph::PairedInfoIndex<Graph>& paired_info) {
		SimpleSequenceMapper<k, Graph> sequence_mapper(g_, index_);
		Path<EdgeId> path = sequence_mapper.MapSequence(sequence);

//		cout << "PATH SIZE " << path.size() << endl;

		for (size_t i = 0; i < path.size(); ++i) {
			EdgeId e = path[i];
			if (g_.length(e) + delta_ > gap_) {
//				cout << "HERE1 " << endl;
				AddEtalonInfo(paired_info, e, e, 0);
			}
			size_t j = i + 1;
			size_t length = 0;

			while (j < path.size() && length  + k + 2 <= (insert_size_ + delta_)) {
				if (length + g_.length(e) + g_.length(path[j]) + delta_ + k
						>= gap_ + 2 * (k + 1)) {
//					cout << "HERE2 " <<  /*g_.length(e) + */length << endl;
					AddEtalonInfo(paired_info, e, path[j],
							g_.length(e) + length);
				}
				length += g_.length(path[j++]);
			}
		}

	}

public:

	EtalonPairedInfoCounter(Graph& g, const EdgeIndex<k + 1, Graph>& index
			, size_t insert_size, size_t read_length, size_t delta) :
			g_(g), index_(index), insert_size_(insert_size), read_length_(
					read_length), gap_(insert_size_ - 2 * read_length_), delta_(
					delta) {
		assert(insert_size_ >= 2 * read_length_);
//		cout << "IS " << insert_size_ << endl;
//		cout << "RL " << read_length_ << endl;
//		cout << "GAP " << gap_ << endl;
//		cout << "DELTA " << delta_ << endl;
	}

	void FillEtalonPairedInfo(const Sequence& genome,
			omnigraph::PairedInfoIndex<Graph>& paired_info) {
		ProcessSequence(genome, paired_info);
		ProcessSequence(!genome, paired_info);
	}
};

template<size_t k, class Graph>
class NewEtalonPairedInfoCounter {
	typedef typename Graph::EdgeId EdgeId;

	Graph& g_;
	const EdgeIndex<k + 1, Graph>& index_;
	size_t insert_size_;
	size_t read_length_;
	size_t gap_;
	size_t delta_;

	void AddEtalonInfo(set<PairInfo<EdgeId>> paired_info, EdgeId e1, EdgeId e2,
			double d) {
		PairInfo<EdgeId> pair_info(e1, e2, d, 1000.0);
		paired_info.insert(pair_info);
	}

	void ProcessSequence(const Sequence& sequence,
			set<PairInfo<EdgeId>>& temporary_info) {
		int mod_gap = (gap_ > delta_) ? gap_ - delta_ : 0;
		Seq<k + 1> left(sequence);
		left >> 0;
		for (size_t left_idx = 0; left_idx + k + 1 + mod_gap <= sequence.size(); ++left_idx) {
			left = left << sequence[left_idx + k];
			size_t right_idx = left_idx + mod_gap;
			if (!index_.containsInIndex(left)) {
				continue;
			}
			pair<EdgeId, size_t> left_pos = index_.get(left);
			Seq<k + 1> right(sequence, right_idx);
			right >> 0;
			for (; right_idx + k + 1 <= left_idx + insert_size_ + delta_ && right_idx + k + 1 <= sequence.size(); ++right_idx) {
				right << sequence[right_idx + k];
				if (!index_.containsInIndex(right)) {
					continue;
				}
				pair<EdgeId, size_t> right_pos = index_.get(right);
				AddEtalonInfo(temporary_info, left_pos.first, right_pos.first, right_idx - left_idx + left_pos.second - right_pos.second);
			}
		}
	}

public:

	NewEtalonPairedInfoCounter(Graph& g, const EdgeIndex<k + 1, Graph>& index
			, size_t insert_size, size_t read_length, size_t delta) :
			g_(g), index_(index), insert_size_(insert_size), read_length_(
					read_length), gap_(insert_size_ - 2 * read_length_), delta_(
					delta) {
		assert(insert_size_ >= 2 * read_length_);
//		cout << "IS " << insert_size_ << endl;
//		cout << "RL " << read_length_ << endl;
//		cout << "GAP " << gap_ << endl;
//		cout << "DELTA " << delta_ << endl;
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

template<class Graph>
class UniqueDistanceStat: public omnigraph::AbstractStatCounter {
	typedef omnigraph::PairedInfoIndex<Graph> PairedIndex;

	PairedIndex& paired_info_;
	size_t unique_;
	size_t non_unique_;
public:

	UniqueDistanceStat(PairedIndex& paired_info) :
			paired_info_(paired_info), unique_(0), non_unique_(0) {

	}

	virtual ~UniqueDistanceStat() {

	}

	virtual void Count() {
		for (auto it = paired_info_.begin(); it != paired_info_.end(); ++it) {
			assert((*it).size() > 0);
			if ((*it).size() > 1) {
				non_unique_++;
//				for (auto info_it = (*it).begin(); info_it != (*it).end(); ++info_it) {
//					//todo
//				}
			} else {
				unique_++;
			}
		}INFO(unique_ << " unique edge distances");
		INFO(non_unique_ << " non unique edge distances");
	}

	size_t unique() {
		return unique_;
	}

	size_t non_unique() {
		return non_unique_;
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
		SimpleSequenceMapper<k, Graph> sequence_mapper(graph, index);
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

template<class Graph, size_t kmer_size, class Stream>
class PairedIndexFiller {
private:
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	const EdgeIndex<kmer_size + 1, Graph>& index_;
	Stream& stream_;

	inline size_t CountDistance(const io::PairedRead& paired_read) {
		return paired_read.distance() - paired_read.second().size();
	}

	size_t CorrectLength(Path<EdgeId> path, size_t idx) {
		size_t answer = graph_.length(path[idx]);
		if (idx == 0)
			answer -= path.start_pos();
		if (idx == path.size() - 1)
			answer -= graph_.length(path[idx]) - path.end_pos();
		return answer;
	}

	void ProcessPairedRead(
			omnigraph::PairedInfoIndex<Graph> &paired_index,
			const io::PairedRead& p_r,
			debruijn_graph::SimpleSequenceMapper<kmer_size, Graph> &read_threader) {
		Sequence read1 = p_r.first().sequence();
		Sequence read2 = p_r.second().sequence();
		Path<EdgeId> path1 = read_threader.MapSequence(read1);
		Path<EdgeId> path2 = read_threader.MapSequence(read2);
		size_t distance = CountDistance(p_r);
		int current_distance1 = distance + path1.start_pos()
				- path2.start_pos();
		for (size_t i = 0; i < path1.size(); ++i) {
			int current_distance2 = current_distance1;
			for (size_t j = 0; j < path2.size(); ++j) {
				double weight = CorrectLength(path1, i)
						* CorrectLength(path2, j);
				PairInfo<EdgeId> new_info(path1[i], path2[j], current_distance2,
						weight);
				paired_index.AddPairInfo(new_info);
				current_distance2 += graph_.length(path2[j]);
			}
			current_distance1 -= graph_.length(path1[i]);
		}
	}

public:

	PairedIndexFiller(Graph &graph, const EdgeIndex<kmer_size + 1, Graph>& index
			, Stream& stream) :
			graph_(graph), index_(index), stream_(stream) {

	}

	/**
	 * Method reads paired data from stream, maps it to genome and stores it in this PairInfoIndex.
	 */
	void FillIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			paired_index.AddPairInfo(PairInfo<EdgeId>(*it, *it, 0, 0.0));
		}
		typedef Seq<kmer_size + 1> KPOMer;
		debruijn_graph::SimpleSequenceMapper<kmer_size, Graph> read_threader(
				graph_, index_);
		stream_.reset();
		while (!stream_.eof()) {
      io::PairedRead p_r;
			stream_ >> p_r;
			ProcessPairedRead(paired_index, p_r, read_threader);
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
	SimpleSequenceMapper<k, Graph> read_seq_mapper;
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
