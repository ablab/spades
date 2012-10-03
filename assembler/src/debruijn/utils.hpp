//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "openmp_wrapper.h"

#include "io/paired_read.hpp"
#include "omni/omni_utils.hpp"
#include "omni/id_track_handler.hpp"
#include "omni/splitters.hpp"
#include "omni/path_processor.hpp"
#include "omni/distance_estimation.hpp"
#include "omni/paired_info.hpp"

#include "logger/logger.hpp"
#include "xmath.h"
#include "sequence/sequence_tools.hpp"

#include "runtime_k.hpp"

#include <iostream>

#include "kmer_map.hpp"
#include "new_debruijn.hpp"
#include "debruijn_kmer_index.hpp"
#include "edge_index.hpp"
#include "sequence_mapper.hpp"
#include "standard.hpp"

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
        DEBUG("Read distance " << read_distance);
		for (size_t i = 0; i < path1.size(); ++i) {
			pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
			for (size_t j = 0; j < path2.size(); ++j) {
				pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
				double weight = weight_f_(mapping_edge_1.second,
						mapping_edge_2.second);
				size_t kmer_distance = read_distance
						+ mapping_edge_2.second.initial_range.end_pos
						- mapping_edge_1.second.initial_range.start_pos;
				int edge_distance = kmer_distance
						+ mapping_edge_1.second.mapped_range.start_pos
						- mapping_edge_2.second.mapped_range.end_pos;

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

	virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
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
