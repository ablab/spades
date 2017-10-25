//****************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "assembly_graph/graph_support/contig_output.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "modules/simplification/ec_threshold_finder.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "chromosome_removal.hpp"

#include "math/xmath.h"


namespace debruijn_graph {

//TODO replace with standard methods
void ChromosomeRemoval::CompressAll(Graph &g) {
    for (auto it = g.SmartVertexBegin(); ! it.IsEnd(); ++it) {
        if (g.IsDeadStart(*it) && g.IsDeadEnd(*it)) {
            g.DeleteVertex(*it);
        } else {
            g.CompressVertex(*it);
        }
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
    INFO("Removing of long chromosomal edges started");
    CoverageUniformityAnalyzer coverage_analyzer(gp.g, long_edge_bound);
    size_t total_len = coverage_analyzer.TotalLongEdgeLength();
    if (total_len == 0) {
        if (external_chromosome_coverage < 1.0) {
            WARN("plasmid detection failed, not enough long edges");
        }
        else {
            INFO("All long edges deleted, stopping detection");
        }
        return 0;
    }

    double median_long_edge_coverage;
    if (external_chromosome_coverage < 1.0) {
        median_long_edge_coverage = coverage_analyzer.CountMedianCoverage();
        double fraction = coverage_analyzer.UniformityFraction(coverage_limits, median_long_edge_coverage);
        if (math::gr(0.8, fraction)) {
            WARN ("More than 20% of long edges have coverage significantly different from median (total " << size_t ((1-fraction) * 0.5 * double(total_len)) <<" of "<< size_t (double(total_len) * 0.5) << " bases).");
            WARN ("In most cases it means that either read coverage is uneven or significant contamination is present - both of these two cases make plasmidSPAdes' results unreliable");
            WARN ("However, that situation may still be OK if you expect to see large plasmids in your dataset, so plasmidSPAdes will continue to work");
        } else {
            INFO(size_t((1 - fraction) * 100) << "% of bases from long edges have coverage significantly different from median");
        }
        for (auto iter = gp.g.ConstEdgeBegin(); ! iter.IsEnd(); ++iter) {
            if (long_component_.find(*iter) == long_component_.end()) {
                CalculateComponentSize(*iter, gp.g);
            }
        }
        INFO("Connected components calculated");
    } else {
        median_long_edge_coverage = external_chromosome_coverage;
    }
    size_t deleted = 0;
    for (auto iter = gp.g.SmartEdgeBegin(); ! iter.IsEnd(); ++iter){
        if (gp.g.length(*iter) > long_edge_bound) {
            if (gp.g.coverage(*iter) < median_long_edge_coverage * (1 + coverage_limits) && gp.g.coverage(*iter)  > median_long_edge_coverage * (1 - coverage_limits)) {
                DEBUG("Considering long edge: id " << gp.g.int_id(*iter) << " length: " << gp.g.length(*iter) <<" coverage: " << gp.g.coverage(*iter));
                if ( long_component_.find(*iter) != long_component_.end() && 300000 > long_component_[*iter] && deadends_count_[*iter] == 0) {
                    DEBUG("Edge " << gp.g.int_id(*iter) << " skipped - because of small nondeadend connected component of size " << long_component_[*iter]);
                } else {
                    DEBUG(" Edge " << gp.g.int_id(*iter) << "  deleted");
                    deleted++;
                    gp.g.DeleteEdge(*iter);
                }
            }
        }
    }
    INFO("Deleted " << deleted <<" long edges");
    CompressAll(gp.g);
    return median_long_edge_coverage;
}

void ChromosomeRemoval::PlasmidSimplify(conj_graph_pack &gp, size_t long_edge_bound,
                                        std::function<void (EdgeId)> removal_handler ) {
    DEBUG("Simplifying graph for plasmid project");
    size_t iteration_count = 10;
    for (size_t i = 0; i < iteration_count; i++) {
        omnigraph::EdgeRemovingAlgorithm<Graph> tc(gp.g, func::And(DeadEndCondition<Graph>(gp.g), LengthUpperBound<Graph>(gp.g, long_edge_bound)),
                                                   removal_handler, true);
        tc.Run();
    }
    gp.EnsureIndex();
}

void ChromosomeRemoval::run(conj_graph_pack &gp, const char*) {
    //FIXME Seriously?! cfg::get().ds like hundred times...
    OutputEdgeSequences(gp.g, cfg::get().output_dir + "before_chromosome_removal");
    INFO("Before iteration " << 0 << ", " << gp.g.size() << " vertices in graph");
    double chromosome_coverage = RemoveLongGenomicEdges(gp, cfg::get().pd->long_edge_length, cfg::get().pd->relative_coverage );
    PlasmidSimplify(gp, cfg::get().pd->long_edge_length);
//TODO:: reconsider and move somewhere(not to config)
    size_t max_iteration_count = 30;

    for (size_t i = 0; i < max_iteration_count; i++) {
        size_t graph_size = gp.g.size();
        INFO("Before iteration " << i + 1 << " " << graph_size << " vertices in graph");
        RemoveLongGenomicEdges(gp, cfg::get().pd->long_edge_length, cfg::get().pd->relative_coverage, chromosome_coverage );
        INFO("Before dead_end simplification " << i << " " << gp.g.size() << " vertices in graph");

        PlasmidSimplify(gp, cfg::get().pd->long_edge_length);
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
                && old_vertex_weights[gp.g.EdgeStart(*iter)] > long_component_[*iter] + cfg::get().pd->long_edge_length * 2)  {
                DEBUG("deleting isolated edge of length" << gp.g.length(*iter));
                gp.g.DeleteEdge(*iter);
            }
        }
        for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (long_component_[*iter] < 2 * cfg::get().pd->small_component_size) {
                if (old_vertex_weights.find(gp.g.EdgeStart(*iter)) != old_vertex_weights.end() &&
                    old_vertex_weights[gp.g.EdgeStart(*iter)] >
                    long_component_[*iter] + cfg::get().pd->long_edge_length * 2 &&
                        gp.g.coverage(*iter) < chromosome_coverage * (1 + cfg::get().pd->small_component_relative_coverage)
                       && gp.g.coverage(*iter) > chromosome_coverage * (1 - cfg::get().pd->small_component_relative_coverage)) {
                    DEBUG("Deleting edge from fake small component, length " << gp.g.length(*iter) << " component_size " << old_vertex_weights[gp.g.EdgeStart(*iter)]) ;
                    gp.g.DeleteEdge(*iter);
                }
            }
        }
        for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (long_component_[*iter] < 2 * cfg::get().pd->min_component_length &&
                                  !(deadends_count_[*iter] == 0 &&
                                    gp.g.length(*iter) > cfg::get().pd->min_isolated_length)) {
                gp.g.DeleteEdge(*iter);
            }
        }

        CompressAll(gp.g);
        PlasmidSimplify(gp, cfg::get().pd->long_edge_length);
        size_t new_graph_size = gp.g.size();
        if (new_graph_size == graph_size) {
            INFO("Iteration " << i << " of small components additional filtering graph was not changed");
            if (new_graph_size == 0) {
                WARN("No putative plasmid contigs found!");
            } else {
                INFO("After plasmidSPAdes subroutine " << new_graph_size << " vertices left");
            }
            break;
        }
    }
    INFO("Counting average coverage after genomic edge removal");
    AvgCovereageCounter<Graph> cov_counter(gp.g);
    cfg::get_writable().ds.average_coverage = cov_counter.Count();
    INFO("Average coverage = " << cfg::get().ds.average_coverage);
}


} //debruijn_graph
