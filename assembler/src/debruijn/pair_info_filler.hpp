/*
 * pair_info_filler.hpp
 *
 *  Created on: Oct 3, 2013
 *      Author: andrey
 */

#ifndef PAIR_INFO_FILLER_HPP_
#define PAIR_INFO_FILLER_HPP_

#include "sequence_mapper_notifier.hpp"

namespace debruijn_graph {

/**
 * As for now it ignores sophisticated case of repeated consecutive
 * occurrence of edge in path due to gaps in mapping
 *
 * todo talk with Anton about simplification and speed-up of procedure with little quality loss
 */
class LatePairedIndexFiller : public SequenceMapperListener {
    typedef boost::function<double(MappingRange, MappingRange)> WeightF;
    typedef std::pair<EdgeId, EdgeId> EdgePair;
public:
    LatePairedIndexFiller(const Graph &graph, WeightF weight_f, omnigraph::de::UnclusteredPairedInfoIndexT<Graph>& paired_index)
            : graph_(graph),
              weight_f_(weight_f),
              paired_index_(paired_index) {
    }

    virtual void StartProcessLibrary(size_t threads_count) {
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it)
            paired_index_.AddPairInfo(*it, *it, { 0., 0. });
        for (size_t i = 0; i < threads_count; ++i)
            buffer_pi_.emplace_back();
    }

    virtual void StopProcessLibrary() {
        for (size_t i = 0; i < buffer_pi_.size(); ++i)
            MergeBuffer(i);

        buffer_pi_.clear();
    }

    virtual void ProcessPairedRead(size_t thread_index,
                                   const MappingPath<EdgeId>& read1,
                                   const MappingPath<EdgeId>& read2,
                                   size_t dist) {
        ProcessPairedRead(buffer_pi_[thread_index], read1, read2, dist);
    }

    virtual void ProcessSingleRead(size_t,
                                   const MappingPath<EdgeId>&) {}

    virtual void MergeBuffer(size_t thread_index) {
        paired_index_.AddAll(buffer_pi_[thread_index]);
        buffer_pi_[thread_index].Clear();
    }

    virtual ~LatePairedIndexFiller() {}

private:
    EdgePair ConjugatePair(EdgePair ep) const {
        return std::make_pair(graph_.conjugate(ep.second), graph_.conjugate(ep.first));
    }

    void ProcessPairedRead(omnigraph::de::PairedInfoBuffer<Graph>& paired_index,
                           const MappingPath<EdgeId>& path1,
                           const MappingPath<EdgeId>& path2, size_t read_distance) const {
        for (size_t i = 0; i < path1.size(); ++i) {
            std::pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
            for (size_t j = 0; j < path2.size(); ++j) {
                std::pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];

                EdgePair ep{mapping_edge_1.first, mapping_edge_2.first};

                if (ep > ConjugatePair(ep))
                    continue;

                double weight = weight_f_(mapping_edge_1.second,
                                          mapping_edge_2.second);
                size_t kmer_distance = read_distance
                        + mapping_edge_2.second.initial_range.end_pos
                        - mapping_edge_1.second.initial_range.start_pos;
                int edge_distance = (int) kmer_distance
                        + (int) mapping_edge_1.second.mapped_range.start_pos
                        - (int) mapping_edge_2.second.mapped_range.end_pos;

                paired_index.AddPairInfo(mapping_edge_1.first, mapping_edge_2.first,
                                         { (double) edge_distance, weight });
            }
        }
    }

private:
    const Graph& graph_;
    WeightF weight_f_;
    omnigraph::de::UnclusteredPairedInfoIndexT<Graph>& paired_index_;
    std::vector<omnigraph::de::PairedInfoBuffer<Graph> > buffer_pi_;

    DECL_LOGGER("LatePairedIndexFiller")
    ;
};


}


#endif /* PAIR_INFO_FILLER_HPP_ */
