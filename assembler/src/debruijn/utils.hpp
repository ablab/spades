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
#include "sequence_mapper_notifier.hpp"

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

double PairedReadCountWeight(const MappingRange&, const MappingRange&) {
	return 1.;
}

double KmerCountProductWeight(const MappingRange& mr1,
        const MappingRange& mr2) {
    return (double)(mr1.initial_range.size() * mr2.initial_range.size());
}

ConjugateDeBruijnGraph::EdgeId conj_wrap(ConjugateDeBruijnGraph& g,
		ConjugateDeBruijnGraph::EdgeId e) {
	return g.conjugate(e);
}

NonconjugateDeBruijnGraph::EdgeId conj_wrap(NonconjugateDeBruijnGraph& /*g*/,
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
class LatePairedIndexFiller : public SequenceMapperListener {

    typedef boost::function<double(MappingRange, MappingRange)> WeightF;

public:
    LatePairedIndexFiller(const Graph &graph, WeightF weight_f, omnigraph::de::PairedInfoIndexT<Graph>& paired_index)
            : graph_(graph),
              weight_f_(weight_f),
              paired_index_(paired_index) {
    }

    virtual void StartProcessLibrary(size_t threads_count) {
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            paired_index_.AddPairInfo(*it, *it, 0., 0., 0.);
        }
        for (size_t i = 0; i < threads_count; ++i) {
            buffer_pi_.push_back(
                    new omnigraph::de::PairedInfoIndexT<Graph>(graph_));
        }
    }
    virtual void StopProcessLibrary() {
        for (size_t i = 0; i < buffer_pi_.size(); ++i) {
            MergeBuffer(i);
            delete buffer_pi_[i];
        }
        buffer_pi_.clear();
    }
    virtual void ProcessPairedRead(size_t thread_index,
                                   const MappingPath<EdgeId>& read1,
                                   const MappingPath<EdgeId>& read2,
                                   size_t dist) {
        ProcessPairedRead(*(buffer_pi_[thread_index]), read1, read2, dist);
    }
    virtual void ProcessSingleRead(size_t thread_index,
                                   const MappingPath<EdgeId>& read) {
    }
    virtual void MergeBuffer(size_t thread_index) {
        paired_index_.AddAll(*(buffer_pi_[thread_index]));
        buffer_pi_[thread_index]->Clear();
    }
    virtual ~LatePairedIndexFiller() {
    }

private:
    void ProcessPairedRead(omnigraph::de::PairedInfoIndexT<Graph>& paired_index,
                           const MappingPath<EdgeId>& path1,
                           const MappingPath<EdgeId>& path2, size_t read_distance) {
        for (size_t i = 0; i < path1.size(); ++i) {
            pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
            for (size_t j = 0; j < path2.size(); ++j) {
                pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
                double weight = weight_f_(mapping_edge_1.second,
                                          mapping_edge_2.second);
                size_t kmer_distance = read_distance
                        + mapping_edge_2.second.initial_range.end_pos
                        - mapping_edge_1.second.initial_range.start_pos;
                int edge_distance = (int) kmer_distance
                        + (int) mapping_edge_1.second.mapped_range.start_pos
                        - (int) mapping_edge_2.second.mapped_range.end_pos;

                paired_index.AddPairInfo(mapping_edge_1.first,
                                         mapping_edge_2.first,
                                         (double) edge_distance, weight, 0.);
            }
        }
    }

private:
    const Graph& graph_;
    WeightF weight_f_;
    omnigraph::de::PairedInfoIndexT<Graph>& paired_index_;
    vector<omnigraph::de::PairedInfoIndexT<Graph>*> buffer_pi_;

    DECL_LOGGER("LatePairedIndexFiller")
    ;
};

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

double UnityFunction(int /*x*/) {
    return 1.;
}

//postprocessing, checking that clusters do not intersect
template<class Graph>
void RefinePairedInfo(const Graph& graph, PairedInfoIndexT<Graph>& clustered_index) {
  for (auto iter = clustered_index.begin(); iter != clustered_index.end(); ++iter) {
    EdgeId first_edge = iter.first();
    EdgeId second_edge = iter.second();
    const de::Histogram& infos = *iter;
    if (infos.size() == 0)
        continue;

    auto prev_it = infos.begin();
    auto it = prev_it;
    ++it;
    for (auto end_it = infos.end(); it != end_it; ++it) {
      if (math::le(abs(it->d - prev_it->d), it->var + prev_it->var)) {
        WARN("Clusters intersect, edges -- " << graph.int_id(first_edge)
             << " " << graph.int_id(second_edge));
        INFO("Trying to handle this case");
        // seeking the symmetric pair info to [i - 1]
        bool success = false;
        double total_weight = prev_it->weight;
        for (auto inner_it = it; inner_it != end_it; ++inner_it) {
          total_weight += inner_it->weight;
          if (math::eq(inner_it->d + prev_it->d, 0.)) {
            success = true;
            double center = 0.;
            double var = inner_it->d + inner_it->var;
            for (auto inner_it_2 = prev_it; inner_it_2 != inner_it; ++inner_it_2)
            {
              TRACE("Removing pair info " << *inner_it_2);
              clustered_index.RemovePairInfo(first_edge, second_edge, *inner_it_2);
            }
            clustered_index.RemovePairInfo(first_edge, second_edge, *inner_it);
            Point new_point(center, total_weight, var);
            TRACE("Adding new pair info " << first_edge << " " << second_edge << " " << new_point);
            clustered_index.AddPairInfo(first_edge, second_edge, new_point);
            break;
          }
        }
        INFO("Pair information was resolved");

        if (!success)
          WARN("This intersection can not be handled in the right way");
        break;
      }
    }
  }
}

}
