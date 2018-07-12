//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <algorithm>
#include <vector>
#include <set>

#include "assembly_graph/index/edge_multi_index.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"

#include "modules/alignment/edge_index_refiller.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "modules/alignment/gap_info.hpp"

#include "modules/alignment/pacbio/pacbio_read_structures.hpp"
#include "modules/alignment/pacbio/gap_filler.hpp"

namespace sensitive_aligner {

//TODO:: invent appropriate name, move code to .cpp
class PacBioMappingIndex {
public:
    typedef std::set<QualityRange> RangeSet;
    typedef std::pair<QualityRange, int> ColoredRange;

    typedef debruijn_graph::Graph Graph;
    typedef debruijn_graph::VertexId VertexId;
    typedef debruijn_graph::EdgeId EdgeId;


private:
    DECL_LOGGER("PacIndex")

    const Graph &g_;

    static const size_t SHORT_SPURIOUS_LENGTH = 500;
    //presumably separate class for this and GetDistance
    mutable std::map<std::pair<VertexId, VertexId>, size_t> distance_cashed_;
    size_t read_count_;
    debruijn_graph::config::pacbio_processor pb_config_;

    alignment::BWAReadMapper<Graph> bwa_mapper_;

    bool similar_in_graph(const MappingInstance &a, const MappingInstance &b,
                          int shift = 0) const {
        if (b.read_position + shift < a.read_position) {
            return similar_in_graph(b, a, -shift);
        } else if (b.read_position == a.read_position) {
            return (abs(int(b.edge_position) + shift - int(a.edge_position)) < 2);
        } else {
//3 to allow small deletion in read on the consecutive edges
            return ((b.edge_position + shift - a.edge_position) * pb_config_.compression_cutoff <= max((b.read_position - a.read_position), 3));
        }
    }

    RangeSet GetBWAClusters(const io::SingleRead &read) const {
        DEBUG("BWA started")
        RangeSet res;
        Sequence s = read.sequence();
        if (s.size() < g_.k()) {
            return res;
        }

        omnigraph::MappingPath<EdgeId> mapped_path = FilterShortAlignments(FilterSpuriousAlignments(bwa_mapper_.MapSequence(s), s.size()));

        TRACE(read_count_ << " read_count_");
        TRACE("BWA ended")
        DEBUG(mapped_path.size() <<"  clusters");
        for (const auto &e_mr : mapped_path) {
            EdgeId e = e_mr.first;
            omnigraph::MappingRange mr = e_mr.second;
            DEBUG("ReadName=" << read.name() << " BWA loading edge=" << g_.int_id(e) << " e_start=" << mr.mapped_range.start_pos << " e_end=" << mr.mapped_range.end_pos
                              << " r_start=" << mr.initial_range.start_pos << " r_end=" << mr.initial_range.end_pos << " qual " << mr.quality <<" len "<< g_.length(e) );
            size_t cut = 0;
            size_t edge_start_pos = mr.mapped_range.start_pos;
            size_t edge_end_pos = mr.mapped_range.end_pos;
            size_t read_start_pos = mr.initial_range.start_pos + cut;
            size_t read_end_pos = mr.initial_range.end_pos;
            if (edge_start_pos >= edge_end_pos || read_start_pos >= read_end_pos) {
                DEBUG ("skipping extra-short alignment");
                continue;
            }
            res.insert(QualityRange(e, edge_start_pos, edge_end_pos, read_start_pos, read_end_pos, mr.quality));
        }
        DEBUG("Ended loading bwa")
        return res;
    }

    typename omnigraph::MappingPath<EdgeId> FilterShortAlignments(typename omnigraph::MappingPath<EdgeId> mapped_path) const {
        omnigraph::MappingPath<EdgeId> res;
        size_t length_cutoff = pb_config_.internal_length_cutoff;

        for (size_t i = 0; i < mapped_path.size(); i++) {
            size_t rlen = mapped_path[i].second.initial_range.size();
//left and right ends of ranges;
//TODO:: think whether it is right condition
            if (rlen > length_cutoff/* || (rlen > g_.k() && (i == mini || i == maxi))*/){
                vector<pair<size_t, int> > range_limits;
                for (size_t j = 0; j < mapped_path.size(); j++) {
                    if (i != j) {
                        if (mapped_path[i].second.initial_range.Intersect(mapped_path[j].second.initial_range)) {
                            size_t pos_start = std::max (mapped_path[i].second.initial_range.start_pos, mapped_path[j].second.initial_range.start_pos)
                                               - mapped_path[i].second.initial_range.start_pos;
                            size_t pos_end = std::min (mapped_path[i].second.initial_range.end_pos, mapped_path[j].second.initial_range.end_pos)
                                             - mapped_path[i].second.initial_range.start_pos;
                            range_limits.push_back(make_pair(pos_start, 1));
                            range_limits.push_back(make_pair(pos_end, -1));
                        }
                    }
                }
                sort(range_limits.begin(), range_limits.end());
                size_t current_cover = 0;
                if (range_limits.size() > 0) {
//may be negative if some ranges are zero-sized
                    int covered_lays = range_limits[0].second;
                    for (size_t j = 1; j < range_limits.size(); j++) {
                        if (covered_lays > 0) {
                            current_cover += range_limits[j].first - range_limits[j - 1].first;
                        }
                        covered_lays += range_limits[j].second;
                    }
                }
                if (current_cover * 2 <  rlen)
                    res.push_back(mapped_path[i].first, mapped_path[i].second);
                else {
                    DEBUG ("Filtering " << g_.int_id(mapped_path[i].first) << " " << mapped_path[i].second)
                }
            }
        }
        return res;
    }

    typename omnigraph::MappingPath<EdgeId> FilterSpuriousAlignments(typename omnigraph::MappingPath<EdgeId> mapped_path, size_t seq_len) const {
        omnigraph::MappingPath<EdgeId> res;

        for (size_t i = 0; i < mapped_path.size(); i++) {
            omnigraph::MappingRange current = mapped_path[i].second;
//Currently everything in kmers
            size_t expected_additional_left = std::min(current.initial_range.start_pos, current.mapped_range.start_pos);
            size_t expected_additional_right = std::min(seq_len - current.initial_range.end_pos - g_.k(),
                                                        g_.length(mapped_path[i].first) - current.mapped_range.end_pos);

            size_t rlen =
                    mapped_path[i].second.initial_range.end_pos - mapped_path[i].second.initial_range.start_pos;

//FIXME more reasonable condition
//g_.k() to compare sequence lengths, not kmers
            if (rlen < SHORT_SPURIOUS_LENGTH && (rlen + g_.k()) * 2 < expected_additional_left + expected_additional_right) {
                DEBUG ("Skipping spurious alignment " << i << " on edge " << mapped_path[i].first);
            } else
                res.push_back(mapped_path[i].first, mapped_path[i].second);
        }

//        if (res.size() != mapped_path.size()) {
        if (true) {
            DEBUG("Seq len " << seq_len);
            for (const auto &e_mr : mapped_path) {
                EdgeId e = e_mr.first;
                omnigraph::MappingRange mr = e_mr.second;
                DEBUG("Alignments were" << g_.int_id(e) << " e_start=" << mr.mapped_range.start_pos << " e_end=" <<
                                        mr.mapped_range.end_pos
                                        << " r_start=" << mr.initial_range.start_pos << " r_end=" << mr.initial_range.end_pos);
            }
        }
        return res;
    }

    std::vector<std::vector<bool>> FillConnectionsTable(const RangeSet &mapping_descr) const {
        size_t len =  mapping_descr.size();
        TRACE("getting colors, table size "<< len);
        std::vector<std::vector<bool>> cons_table(len);
        for (size_t i = 0; i < len; i++) {
            cons_table[i].resize(len);
            cons_table[i][i] = 0;
        }
        size_t i = 0;
        for (auto i_iter = mapping_descr.begin(); i_iter != mapping_descr.end();
             ++i_iter, ++i) {
            size_t j = i;
            for (auto j_iter = i_iter;
                 j_iter != mapping_descr.end(); ++j_iter, ++j) {
                if (i_iter == j_iter)
                    continue;
                cons_table[i][j] = IsConsistent(*i_iter, *j_iter);
            }
        }
        return cons_table;
    }


//Currently unused but useful for debug purposes
    std::string DebugEmptyBestScoredPath(VertexId start_v, VertexId end_v, EdgeId prev_edge, EdgeId cur_edge,
                                         size_t prev_last_edge_position, size_t cur_first_edge_position, int seq_len) const {
        size_t result = GetDistance(start_v, end_v, /*update cache*/false);
        std::ostringstream ss;
        ss << "Tangled region between edges " << g_.int_id(prev_edge) << " " << g_.int_id(cur_edge) <<  " is not closed, additions from edges: "
           << int(g_.length(prev_edge)) - int(prev_last_edge_position) <<" " << int(cur_first_edge_position)
           << " and seq "<< seq_len << " and shortest path " << result;
        return ss.str();
    }

public:
    PacBioMappingIndex(const Graph &g,
                        debruijn_graph::config::pacbio_processor pb_config, 
                        alignment::BWAIndex::AlignmentMode mode)
            : g_(g),
              pb_config_(pb_config),
              bwa_mapper_(g, mode, pb_config.bwa_length_cutoff) {
        DEBUG("PB Mapping Index construction started");
        DEBUG("Index constructed");
        read_count_ = 0;
    }

//former GetWeightedColors
    std::vector<ColoredRange> GetRangedColors(const io::SingleRead &read) const {
        Sequence s = read.sequence();
        RangeSet mapping_descr = GetBWAClusters(read);

        size_t len = mapping_descr.size();
        std::vector<int> colors(len, UNDEF_COLOR);
        std::vector<double> cluster_size(len) ;
        size_t ii = 0;
        for (const auto &cl : mapping_descr) {
            cluster_size[ii++] = cl.size * cl.quality;
        }
        const auto cons_table = FillConnectionsTable(mapping_descr);

        std::vector<double> max_size(len);
        std::vector<size_t> prev(len);

        int cur_color = 0;
        int num_colors = 0;
        while (true) {
            for (size_t i = 0; i < len; ++i) {
                max_size[i] = 0;
                prev[i] = size_t(-1);
            }
            for (size_t i = 0; i < len; ++i) {
                if (colors[i] != UNDEF_COLOR) continue;
                max_size[i] = cluster_size[i];
                for (size_t j = 0; j < i; j++) {
                    if (colors[j] != -1) continue;
                    if (cons_table[j][i] && math::ls(max_size[i], cluster_size[i] + max_size[j])) {
                        max_size[i] = max_size[j] + cluster_size[i];
                        prev[i] = j;
                    }
                }
            }
            double maxx = 0;
            int maxi = -1;
            for (size_t j = 0; j < len; j++) {
                if (math::gr(max_size[j], maxx)) {
                    maxx = max_size[j];
                    maxi = int(j);
                }
            }
            if (maxi == -1) {
                break;
            }
            num_colors ++;
            cur_color = maxi;
            colors[maxi] = cur_color;
            int real_maxi = maxi, min_i = maxi;

            while (prev[maxi] != -1ul) {
                min_i = maxi;
                maxi = int(prev[maxi]);
                colors[maxi] = cur_color;
            }
            while (real_maxi >= min_i) {
                if (colors[real_maxi] == UNDEF_COLOR) {
                    colors[real_maxi] = DELETED_COLOR;
                }
                real_maxi--;
            }
        }
        DEBUG("Num hits clusters=" << num_colors);
        std::vector<ColoredRange> res;
        size_t ind = 0;
        for ( auto i_iter = mapping_descr.begin(); i_iter != mapping_descr.end(); ++i_iter, ++ind) {
            res.push_back(make_pair(*i_iter, colors[ind]));
        }
        return res;
    }

    std::pair<int, int> GetPathLimits(const QualityRange &a,
                                      const QualityRange &b,
                                      int s_add_len, int e_add_len) const {
        int start_pos = a.sorted_positions[a.last_trustable_index].read_position;
        int end_pos = b.sorted_positions[b.first_trustable_index].read_position;
        int seq_len = -start_pos + end_pos;
        //int new_seq_len =
//TODO::something more reasonable
        int path_min_len = std::max(int(floor((seq_len - int(g_.k())) * pb_config_.path_limit_pressing)), 0);
        int path_max_len = (int) ((double) (seq_len + g_.k() * 2) * pb_config_.path_limit_stretching);
        if (seq_len < 0) {
            DEBUG("suspicious negative seq_len " << start_pos << " " << end_pos << " " << path_min_len << " " << path_max_len);
            if (path_max_len < 0)
            return std::make_pair(-1, -1);
        }
        path_min_len = std::max(path_min_len - int(s_add_len + e_add_len), 0);
        path_max_len = std::max(path_max_len - int(s_add_len + e_add_len), 0);
        return std::make_pair(path_min_len, path_max_len);
    }

    size_t GetDistance(VertexId start_v, VertexId end_v,
                       bool update_cache = true) const {
        size_t result = size_t(-1);
        auto vertex_pair = std::make_pair(start_v, end_v);
        bool not_found;
        auto distance_it = distance_cashed_.begin();
        //FIXME should permit multiple readers
#pragma omp critical(pac_index)
        {
            distance_it = distance_cashed_.find(vertex_pair);
            not_found = (distance_it == distance_cashed_.end());
        }
        if (not_found) {
            omnigraph::DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra dijkstra(
                    omnigraph::DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(g_,
                                                                                            pb_config_.max_path_in_dijkstra,
                                                                                            pb_config_.max_vertex_in_dijkstra));
            dijkstra.Run(start_v);
            if (dijkstra.DistanceCounted(end_v)) {
                result = dijkstra.GetDistance(end_v);
            }
            if (update_cache) {
#pragma omp critical(pac_index)
                distance_cashed_.insert({vertex_pair, result});
            }
        } else {
            TRACE("taking from cashed");
            result = distance_it->second;
        }

        return result;
    }

    bool IsConsistent(const QualityRange &a,
                      const QualityRange &b) const {
        EdgeId a_edge = a.edgeId;
        EdgeId b_edge = b.edgeId;
        DEBUG("Checking consistency: " << g_.int_id(a_edge) << " and " << g_.int_id(b_edge));
        //FIXME: Is this check useful?
        if (a.sorted_positions[a.last_trustable_index].read_position +
                    (int) pb_config_.max_path_in_dijkstra <
            b.sorted_positions[b.first_trustable_index].read_position) {
            DEBUG("Clusters are too far in read");
            return false;
        }
        size_t result = GetDistance(g_.EdgeEnd(a_edge), g_.EdgeStart(b_edge));
        DEBUG ("Distance: " << result);
        if (result == size_t(-1)) {
            return false;
        }
        result += g_.length(a_edge);
        if (similar_in_graph(a.sorted_positions[1], b.sorted_positions[0], (int)result)) {
            DEBUG(" similar! ")
            return true;
        } else {
//FIXME: reconsider this condition! i.e only one large range? That may allow to decrease the bwa length cutoff
            if (- a.sorted_positions[1].edge_position +  (int)result + b.sorted_positions[0].edge_position  <=
                        b.sorted_positions[0].read_position - a.sorted_positions[1].read_position + 2 * int(g_.k())
                        && (pb_config_.bwa_length_cutoff > 2*g_.k() || pb_config_.internal_length_cutoff > 2*g_.k()) ) {
                DEBUG("overlapping range magic worked, " << - a.sorted_positions[1].edge_position +
                                                           (int)result + b.sorted_positions[0].edge_position
                      << " and " <<  b.sorted_positions[0].read_position - a.sorted_positions[1].read_position + 2 * g_.k());
                DEBUG("Ranges:" << a.str(g_) << " " << b.str(g_)
                                <<" llength and dijkstra shift :" << g_.length(a_edge) << " " << (result - g_.length(a_edge)));
                return true;
            } else {
                DEBUG("Not similar")
                return false;
            }
        }
    }

};
}
