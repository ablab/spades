//****************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "simplification/graph_simplification.hpp"
#include "omni/ec_threshold_finder.hpp"

#include "chromosome_removal.hpp"
#include "contig_output.hpp"

namespace debruijn_graph {


void ChromosomeRemoval::CompressAll(Graph &g) {
    for (auto it = g.SmartVertexBegin(); ! it.IsEnd(); ++it) {
        if (g.IsDeadStart(*it) && g.IsDeadEnd(*it)) {
            g.DeleteVertex(*it);
        } else {
            g.CompressVertex(*it);
        }
    }
}

void ChromosomeRemoval::DeleteAndCompress(EdgeId e, Graph &g){
    auto start = g.EdgeStart(e);
    auto end = g.EdgeEnd(e);
    g.DeleteEdge(e);
    bool is_cycle = (start == end || start == g.conjugate(end));
    if (g.IsDeadStart(start) && g.IsDeadEnd(start)) {
        g.DeleteVertex(start);
    } else {
        g.CompressVertex(start);
    }
    if (is_cycle) {
        return;
    }
    if (g.IsDeadStart(end) && g.IsDeadEnd(end)) {
        g.DeleteVertex(end);
    } else {
        g.CompressVertex(end);
    }
}


size_t ChromosomeRemoval::CalculateComponentSize(EdgeId e, Graph &g_) {
    std::stack<EdgeId> next;
    size_t deadend_count = 0;
    next.push(e);
    std::unordered_set<EdgeId> used;
    size_t ans = 0;
    while (!next.empty()){
        auto cur = next.top();
        next.pop();
        if (used.find(cur) != used.end()) {
            continue;
        }
        ans += g_.length(cur);

        used.insert(cur);
        vector<EdgeId> neighbours;
        neighbours.push_back(g_.conjugate(cur));
        auto start = g_.EdgeStart(cur);
        auto tmp = g_.IncidentEdges(start);

        neighbours.insert(neighbours.end(), tmp.begin(), tmp.end());
        auto end = g_.EdgeEnd(cur);
        if (g_.IsDeadStart(start))
            deadend_count++;
        if (g_.IsDeadEnd(end))
            deadend_count++;
        tmp = g_.IncidentEdges(end);
        neighbours.insert(neighbours.end(), tmp.begin(), tmp.end());
        for (auto ee:neighbours) {
            if (used.find(ee) == used.end()) {
                next.push(ee);
            }
        }
    }
    for (auto edge: used) {
        long_component_[edge] = ans;
        long_vertex_component_[g_.EdgeStart(edge)] = ans;
        long_vertex_component_[g_.EdgeEnd(edge)] = ans;
        deadends_count_[edge] = deadend_count;
    }
    return ans;
}

double ChromosomeRemoval::RemoveLongGenomicEdges(conj_graph_pack &gp, size_t long_edge_bound, double coverage_limits, double external_chromosome_coverage){
    INFO("Removing long genomic edge started");
    vector <pair<double, size_t> > coverages;
    size_t total_len = 0, short_len = 0, cur_len = 0;
    for (auto iter = gp.g.ConstEdgeBegin(); ! iter.IsEnd(); ++iter){
        if (gp.g.length(*iter) > cfg::get().pd.edge_length_for_median) {
            coverages.push_back(make_pair(gp.g.coverage(*iter), gp.g.length(*iter)));
            total_len += gp.g.length(*iter);
            long_component_[*iter] = 0;
        } else {
            short_len += gp.g.length(*iter);
        }
    }
    if (total_len == 0) {
        WARN("plasmid detection failed, not enough long edges");
        return 0;
    }
    std::sort(coverages.begin(), coverages.end());
    size_t i = 0;
    while (cur_len < total_len/2 && i <coverages.size()) {
        cur_len += coverages[i].second;
        i++;
    }

    double median_long_edge_coverage;
    if (external_chromosome_coverage < 1.0) {
        median_long_edge_coverage = coverages[i-1].first;
        INFO ("genomic coverage is "<< median_long_edge_coverage << " calculated of length " << total_len);
        for (auto iter = gp.g.ConstEdgeBegin(); ! iter.IsEnd(); ++iter) {
            if (long_component_.find(*iter) == long_component_.end()) {
                CalculateComponentSize(*iter, gp.g);
            }
        }
        INFO("Connected components calculated");
    } else {
        median_long_edge_coverage = external_chromosome_coverage;
    }

    for (auto iter = gp.g.SmartEdgeBegin(); ! iter.IsEnd(); ++iter){
        if (gp.g.length(*iter) > long_edge_bound) {
            if (gp.g.coverage(*iter) < median_long_edge_coverage * (1 + coverage_limits) && gp.g.coverage(*iter)  > median_long_edge_coverage * (1 - coverage_limits)) {
                INFO("Considering long edge: id " << gp.g.int_id(*iter) << " length: " << gp.g.length(*iter) <<" coverage: " << gp.g.coverage(*iter));
                if ( long_component_.find(*iter) != long_component_.end() && 300000 > long_component_[*iter] && deadends_count_[*iter] == 0) {
                    INFO("Edge " << gp.g.int_id(*iter) << " skipped - small nondeadend connected component (" << long_component_[*iter] <<  " )" );
                } else {
                    INFO(" Edge " << gp.g.int_id(*iter) << "  deleted");
                    gp.g.DeleteEdge(*iter);
                }
            }
        }
    }
    CompressAll(gp.g);
    return median_long_edge_coverage;
}

void ChromosomeRemoval::PlasmidSimplify(conj_graph_pack &gp, size_t long_edge_bound,
                                        std::function<void (EdgeId)> removal_handler ) {
    INFO("Simplifying graph for plasmid project");
    size_t iteration_count = 10;
    for (size_t i = 0; i < iteration_count; i++) {
        //pred::TypedPredicate<typename Graph::EdgeId> condition = make_shared<LengthUpperBound<Graph>>(gp.g, long_edge_bound) ;
        omnigraph::EdgeRemovingAlgorithm<Graph> tc(gp.g, pred::And(DeadEndCondition<Graph>(gp.g), LengthUpperBound<Graph>(gp.g, long_edge_bound)),
                                                   removal_handler, true);
        tc.Run();
    }
    gp.EnsureIndex();
}

void ChromosomeRemoval::run(conj_graph_pack &gp, const char*) {
    OutputContigs(gp.g, cfg::get().output_dir + "before_chromosome_removal", false, 0, false);
    INFO("Before iteration " << 0 << " " << gp.g.size() << " vertices in graph");
    double chromosome_coverage = RemoveLongGenomicEdges(gp, cfg::get().pd.long_edge_length, cfg::get().pd.relative_coverage );
    PlasmidSimplify(gp, cfg::get().pd.long_edge_length);
//TODO:: reconsider and move somewhere(not to config)
    size_t max_iteration_count = 30;

    for (size_t i = 0; i < max_iteration_count; i++) {
        size_t graph_size = gp.g.size();
        INFO("Before iteration " << i + 1 << " " << graph_size << " vertices in graph");
        RemoveLongGenomicEdges(gp, cfg::get().pd.long_edge_length, cfg::get().pd.relative_coverage, chromosome_coverage );
        INFO("Before dead_end simplification " << i << " " << gp.g.size() << " vertices in graph");

        PlasmidSimplify(gp, cfg::get().pd.long_edge_length);
        size_t new_graph_size = gp.g.size();
        if (new_graph_size == graph_size) {
            INFO("Iteration " << i << " graph was not changed");
            INFO(new_graph_size << " vertices left");
            break;
        }
    }
//Small repetitive components after filtering
    std::unordered_map<VertexId, size_t> old_vertex_weights (long_vertex_component_.begin(), long_vertex_component_.end());
    for (size_t i = 0; i < max_iteration_count; i++) {
        size_t graph_size = gp.g.size();
        long_vertex_component_.clear();
        long_component_.clear();
        deadends_count_.clear();    
        for (auto iter = gp.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            CalculateComponentSize(*iter, gp.g);
        }

        for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (gp.g.IsDeadEnd(gp.g.EdgeEnd(*iter)) && gp.g.IsDeadStart(gp.g.EdgeStart(*iter))
                && old_vertex_weights.find(gp.g.EdgeStart(*iter)) !=  old_vertex_weights.end()
//* 2 - because all coverages are taken with rc
                && old_vertex_weights[gp.g.EdgeStart(*iter)] > long_component_[*iter] + cfg::get().pd.long_edge_length * 2)  {
                INFO("deleting isolated edge of length" << gp.g.length(*iter));

                gp.g.DeleteEdge(*iter);
            }
        }
        for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (long_component_[*iter] < 2 * cfg::get().pd.small_component_size) {
                if (old_vertex_weights.find(gp.g.EdgeStart(*iter)) != old_vertex_weights.end() &&
                    old_vertex_weights[gp.g.EdgeStart(*iter)] >
                    long_component_[*iter] + cfg::get().pd.long_edge_length * 2 &&
                        gp.g.coverage(*iter) < chromosome_coverage * (1 + cfg::get().pd.small_component_relative_coverage)
                       && gp.g.coverage(*iter) > chromosome_coverage * (1 - cfg::get().pd.small_component_relative_coverage)) {
                    INFO("deleting edge from fake small component, length " << gp.g.length(*iter) << " component_size " << old_vertex_weights[gp.g.EdgeStart(*iter)]) ;
                    gp.g.DeleteEdge(*iter);
                }
            }
        }
        for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (long_component_[*iter] < 2 * cfg::get().pd.min_component_length &&
                                  !(deadends_count_[*iter] == 0 &&
                                    gp.g.length(*iter) > cfg::get().pd.min_isolated_length)) {
                gp.g.DeleteEdge(*iter);
            }
        }

        CompressAll(gp.g);
        PlasmidSimplify(gp, cfg::get().pd.long_edge_length);
        size_t new_graph_size = gp.g.size();
        if (new_graph_size == graph_size) {
            INFO("Iteration " << i << " of small components additional filtering graph was not changed");
            INFO(new_graph_size << "vertices left");
            break;
        }
    }
    INFO("Counting average coverage after genomic edge removal");
    AvgCovereageCounter<Graph> cov_counter(gp.g);
    cfg::get_writable().ds.set_avg_coverage(cov_counter.Count());
    INFO("Average coverage = " << cfg::get().ds.avg_coverage());
}


} //debruijn_graph
