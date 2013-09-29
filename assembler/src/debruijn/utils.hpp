//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "openmp_wrapper.h"
#include "standard.hpp"

#include "io/paired_read.hpp"
#include "omni/omni_utils.hpp"
#include "omni/id_track_handler.hpp"
#include "omni/splitters.hpp"
#include "omni/path_processor.hpp"

#include "de/distance_estimation.hpp"
#include "de/paired_info.hpp"

#include "adt/kmer_map.hpp"

#include "logger/logger.hpp"
#include "xmath.h"
#include "sequence/sequence_tools.hpp"

#include "runtime_k.hpp"

#include "path_helper.hpp"

#include "debruijn_graph.hpp"
#include "indices/debruijn_kmer_index.hpp"
#include "edge_index.hpp"
#include "sequence_mapper.hpp"
#include "genomic_quality.hpp"

#include <iostream>

namespace debruijn_graph {

using omnigraph::Path;
using omnigraph::MappingPath;
using omnigraph::Range;
using omnigraph::MappingRange;
using namespace omnigraph::de;
using omnigraph::GraphActionHandler;

//todo rewrite with extended sequence mapper!
template<class Graph, class Index>
class EtalonPairedInfoCounter {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
	const Index& index_;
	const KmerMapper<Graph>& kmer_mapper_;
	size_t k_;

	size_t insert_size_;
	size_t read_length_;
	int gap_;
	size_t delta_;

  void AddEtalonInfo(PairedInfoIndexT<Graph>& index, EdgeId e1, EdgeId e2, double d) {
    index.AddPairInfo(e1, e2, d, 1000., 0.);
	}

  void ProcessSequence(const Sequence& sequence, PairedInfoIndexT<Graph>& index)
  {
		int mod_gap = (gap_ + (int) k_ > (int) delta_ ) ? gap_ - (int) delta_ : 0 - (int) k_;
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
			     right_idx + k_ + 1 <= left_idx + insert_size_ + delta_ && right_idx + k_ + 1 <= sequence.size();
			     ++right_idx) {
				right <<= sequence[right_idx + k_];
				runtime_k::RtSeq right_upd = kmer_mapper_.Substitute(right);
				if (!index_.contains(right_upd)) {
					continue;
				}
				pair<EdgeId, size_t> right_pos = index_.get(right_upd);

				AddEtalonInfo(index, left_pos.first, right_pos.first,
				              0. + (double) right_idx - (double) left_idx +
				              (double) left_pos.second - (double) right_pos.second);
			}
		}
	}

public:
    EtalonPairedInfoCounter(const Graph& g, const Index& index,
                            const KmerMapper<Graph>& kmer_mapper,
                            size_t insert_size, size_t read_length,
                            size_t delta, size_t k)
            : g_(g),
              index_(index),
              kmer_mapper_(kmer_mapper),
              k_(k),
              insert_size_(insert_size),
              read_length_(read_length),
              gap_((int) (insert_size_ - 2 * read_length_)),
              delta_(delta) {
//		VERIFY(insert_size_ >= 2 * read_length_);
    }

    void FillEtalonPairedInfo(const Sequence& genome,
                              omnigraph::de::PairedInfoIndexT<Graph>& paired_info) {
        ProcessSequence(genome, paired_info);
        ProcessSequence(!genome, paired_info);
    }
};

inline double PairedReadCountWeight(const MappingRange&, const MappingRange&) {
	return 1.;
}

inline double KmerCountProductWeight(const MappingRange& mr1,
                              const MappingRange& mr2) {
    return (double)(mr1.initial_range.size() * mr2.initial_range.size());
}

inline ConjugateDeBruijnGraph::EdgeId conj_wrap(ConjugateDeBruijnGraph& g,
                                         ConjugateDeBruijnGraph::EdgeId e) {
	return g.conjugate(e);
}

inline NonconjugateDeBruijnGraph::EdgeId conj_wrap(NonconjugateDeBruijnGraph& /*g*/,
                                                   NonconjugateDeBruijnGraph::EdgeId e) {
	VERIFY(0);
	return e;
}

inline void WrappedSetCoverage(ConjugateDeBruijnGraph& g,
                               ConjugateDeBruijnGraph::EdgeId e, int cov) {
    g.coverage_index().SetCoverage(e, cov);
    g.coverage_index().SetCoverage(g.conjugate(e), cov);
}

inline void WrappedSetCoverage(NonconjugateDeBruijnGraph& g,
                               NonconjugateDeBruijnGraph::EdgeId e, int cov) {
	g.coverage_index().SetCoverage(e, cov);
}

//template<class Graph>
//class EdgeNeighborhoodFinder: public omnigraph::GraphSplitter<Graph> {
//private:
//  typedef typename Graph::EdgeId EdgeId;
//  typedef typename Graph::VertexId VertexId;
//  EdgeId edge_;
//  size_t max_size_;
//  size_t edge_length_bound_;
//  bool finished_;
//public:
//  EdgeNeighborhoodFinder(const Graph &graph, EdgeId edge, size_t max_size
//          , size_t edge_length_bound) :
//          GraphSplitter<Graph>(graph), edge_(edge), max_size_(
//                  max_size), edge_length_bound_(edge_length_bound), finished_(
//                  false) {
//  }
//
//  GraphComponent<Graph> NextComponent() {
//      CountingDijkstra<Graph> cf(this->graph(), max_size_,
//              edge_length_bound_);
//      set<VertexId> result_set;
//      cf.run(this->graph().EdgeStart(edge_));
//      vector<VertexId> result_start = cf.ReachedVertices();
//      result_set.insert(result_start.begin(), result_start.end());
//      cf.run(this->graph().EdgeEnd(edge_));
//      vector<VertexId> result_end = cf.ReachedVertices();
//      result_set.insert(result_end.begin(), result_end.end());
//
//      ComponentCloser<Graph> cc(this->graph(), edge_length_bound_);
//      cc.CloseComponent(result_set);
//
//      finished_ = true;
//      return GraphComponent<Graph>(this->graph(), result_set.begin(), result_set.end());
//  }
//
//  /*virtual*/ bool Finished() {
//      return finished_;
//  }
//};
//
//template<class Graph>
//class EdgeLocalityPrintingRH {
//	typedef typename Graph::EdgeId EdgeId;
//	typedef typename Graph::VertexId VertexId;
//	const Graph& g_;
//	const GraphLabeler<Graph>& labeler_;
//	const string& output_folder_;
//    boost::function<double (EdgeId)>& quality_f_;
////	size_t black_removed_;
////	size_t colored_removed_;
//public:
//	EdgeLocalityPrintingRH(const Graph& g
//			, const GraphLabeler<Graph>& labeler
//			, const string& output_folder
//            , boost::function<double (EdgeId)> quality_f = 0) :
//			g_(g),
//			labeler_(labeler), output_folder_(output_folder),
//            quality_f_(quality_f){
//	}
//
//	void HandleDelete(EdgeId edge) {
//            TRACE("Deleting edge " << g_.str(edge));
//            if (quality_f_ && math::gr(quality_f_(edge), 0.))
//                INFO("EdgeLocalityPrintRH handling the edge with positive quality : " << quality_f_(edge) << " " << g_.str(edge));
//
//            string folder = output_folder_ + "edges_deleted/";
//            path::make_dir(folder);
//            //todo magic constant
//            map<EdgeId, string> empty_coloring;
//            omnigraph::visualization::WriteComponent(g_, EdgeNeighborhood<Graph>(g_, edge, 50, 250),
//            		folder + "edge_" +  ToString(g_.int_id(edge)) + ".dot", empty_coloring, labeler_);
//	}
//
//private:
//	DECL_LOGGER("QualityEdgeLocalityPrintingRH")
//	;
//};

class WeightDEWrapper {
private:

    vector<double> new_hist;
    int left_x;
    int insert_size;

	void ExtendLinear(const std::map<int, size_t> & hist) {
        size_t sum_weight = 0;

        for (auto iter = hist.begin(); iter != hist.end(); ++iter)
            sum_weight += iter->second;
        DEBUG(sum_weight);

        VERIFY(hist.size() > 0);
        auto iter = hist.begin();

        left_x = iter->first;

        int prev = iter->first;
        size_t prev_val = iter->second;

        new_hist.push_back((double)prev_val / (double)sum_weight);
        ++iter;

        for (; iter != hist.end(); ++iter) {
            int x = iter->first;
            size_t y = iter->second;
            double tan = ((double)y - (double)prev_val) / (x - prev);

            VERIFY(prev < x);
            for (int i = prev + 1; i <= x; ++i) {
                new_hist.push_back(((double)prev_val + tan * (i - prev)) / (double)sum_weight);
            }
            prev = x;
            prev_val = y;
            DEBUG("hist " << x << " " << y);
        }
	}

public:
    WeightDEWrapper(const map<int, size_t>& hist, double IS) {
        DEBUG("WeightDEWrapper " << IS);
        insert_size = (int) IS;
        DEBUG("Extending linear");
        ExtendLinear(hist);
    }

    ~WeightDEWrapper() {
    }


    double CountWeight(int x) const {
        int xx = insert_size - left_x + x - 1;

        if (!(xx >= 0 && xx < (int) new_hist.size())) return 0.;
        VERIFY(math::le(new_hist[xx], 1.));
        return 1000. * new_hist[xx];
    }
};

template<class graph_pack, class PairedRead, class ConfigType>
bool RefineInsertSize(const graph_pack& gp,
                      io::ReadStreamVector<io::IReader<PairedRead> >& streams,
                      ConfigType& config,
                      size_t edge_length_threshold) {
  size_t rl;
  double mean;
  double delta;
  double median;
  double mad;
  std::map<size_t, size_t> percentiles;
  std::map<int, size_t> hist;
  // calling default method
  refine_insert_size(streams, gp, edge_length_threshold, rl, mean, delta, median, mad, percentiles, hist);

  if (hist.size() == 0) {
    config.paired_mode = false;
    WARN("Failed to estimate the insert size of paired reads, because none of the paired reads aligned to long edges.");
    WARN("Paired reads will not be used.");
    return false;
  }

  config.ds.set_IS(mean);
  config.ds.set_is_var(delta);
  config.ds.set_median(median);
  config.ds.set_mad(mad);
  config.ds.set_hist(hist);
  INFO("Mean Insert Size = " << mean);
  INFO("Insert Size stddev= " << delta);
  INFO("Median Insert Size = " << median);
  INFO("Insert Size MAD = " << mad);
  DEBUG("Delta_Mad = " << 1.4826 * mad);

  return true;
}

template<class graph_pack, class PairedRead, class DataSet>
bool RefineInsertSizeForLib(const graph_pack& gp,
                      io::ReadStreamVector<io::IReader<PairedRead> >& streams,
                      DataSet& data,
                      size_t edge_length_threshold) {

  std::map<size_t, size_t> percentiles;
  // calling default method
  data.read_length = 0;
  refine_insert_size(streams, gp, edge_length_threshold,
          data.read_length,
          data.mean_insert_size,
          data.insert_size_deviation,
          data.median_insert_size,
          data.insert_size_mad,
          percentiles,
          data.insert_size_distribution);

  if (data.insert_size_distribution.size() == 0) {
    return false;
  }

  return true;
}

inline double UnityFunction(int /*x*/) {
    return 1.;
}

template<class Graph>
Sequence MergeSequences(const Graph& g,
		const vector<typename Graph::EdgeId>& continuous_path) {
	vector < Sequence > path_sequences;
	path_sequences.push_back(g.EdgeNucls(continuous_path[0]));
	for (size_t i = 1; i < continuous_path.size(); ++i) {
		VERIFY(
				g.EdgeEnd(continuous_path[i - 1])
						== g.EdgeStart(continuous_path[i]));
		path_sequences.push_back(g.EdgeNucls(continuous_path[i]));
	}
	return MergeOverlappingSequences(path_sequences, g.k());
}

}
