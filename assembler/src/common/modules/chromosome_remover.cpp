//***************************************************************************
//* Copyright (c) 2015-2019 Saint-Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "chromosome_remover.hpp"

#include "assembly_graph/core/graph_iterators.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "assembly_graph/graph_support/contig_output.hpp"
#include "assembly_graph/graph_support/coverage_uniformity_analyzer.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "paired_info/paired_info.hpp"
#include "pipeline/config_struct.hpp"
#include "sequence/genome_storage.hpp"
#include "visualization/position_filler.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "math/xmath.h"
#include "assembly_graph/paths/bidirectional_path_container.hpp"

namespace debruijn_graph {

using VertexSet = std::unordered_set<VertexId>;
using SmartVertexSet = omnigraph::SmartContainer<VertexSet, Graph>;
using std::vector;

//paths should be removed from namespace path_extend
using path_extend::BidirectionalPath;
using path_extend::PathContainer;

//TODO replace with standard methods
void ChromosomeRemover::CompressAll(Graph &g) {
    for (auto it = g.SmartVertexBegin(); ! it.IsEnd(); ++it) {
        if (g.IsDeadStart(*it) && g.IsDeadEnd(*it)) {
            g.DeleteVertex(*it);
        } else {
            g.CompressVertex(*it);
        }
    }
}

void ChromosomeRemover::FillForbiddenSet(Graph &g, VertexSet &forbidden) {
    for (VertexId v : g.vertices()) {
        if (g.IsDeadStart(v) || g.IsDeadEnd(v)) {
            forbidden.insert(v);
        }
    }
}

size_t ChromosomeRemover::CalculateComponentSize(EdgeId e, const Graph &g) {
    std::unordered_set<EdgeId> next;
    size_t deadend_count = 0;
    next.insert(e);
    std::unordered_set<EdgeId> used;
    size_t ans = 0;
    while (!next.empty()) {
        auto cur = *next.begin();
        next.erase(next.begin());
        if (used.count(cur)) {
            continue;
        }
        ans += g.length(cur);

        used.insert(cur);
        std::vector<EdgeId> neighbours;
        neighbours.push_back(g.conjugate(cur));
        auto start = g.EdgeStart(cur);
        auto tmp = g.IncidentEdges(start);

        neighbours.insert(neighbours.end(), tmp.begin(), tmp.end());
        auto end = g.EdgeEnd(cur);
        if (g.IsDeadStart(start))
            deadend_count++;
        if (g.IsDeadEnd(end))
            deadend_count++;
        tmp = g.IncidentEdges(end);
        neighbours.insert(neighbours.end(), tmp.begin(), tmp.end());
        for (auto ee:neighbours) {
            if (!used.count(ee)) {
                next.insert(ee);
            }
        }
    }
    for (EdgeId edge: used) {
        long_component_[edge] = ans;
        long_vertex_component_[g.EdgeStart(edge)] = ans;
        long_vertex_component_[g.EdgeEnd(edge)] = ans;
        deadends_count_[edge] = deadend_count;
    }
    component_list_.push_back(std::vector<EdgeId>(used.begin(), used.end()));
    return ans;
}

double ChromosomeRemover::RemoveLongGenomicEdges(size_t long_edge_bound, double coverage_limits, double external_chromosome_coverage) {
    INFO("Removing of long chromosomal edges started");
    auto& graph = gp_.get_mutable<Graph>();
    CoverageUniformityAnalyzer coverage_analyzer(graph, long_edge_bound);
    size_t total_len = coverage_analyzer.TotalLongEdgeLength();
    if (total_len == 0) {
        if (external_chromosome_coverage < 1.0) {
            WARN("Plasmid detection failed, not enough long edges");
        } else {
            INFO("All long edges deleted, stopping detection");
        }
        return 0;
    }

    double median_long_edge_coverage;
    if (external_chromosome_coverage < 1.0) {
        median_long_edge_coverage = coverage_analyzer.CountMedianCoverage();
        double fraction = coverage_analyzer.UniformityFraction(coverage_limits, median_long_edge_coverage);
        if (math::gr(0.8, fraction)) {
            WARN("More than 20% of long edges have coverage significantly different from median (total " << size_t ((1-fraction) * 0.5 * double(total_len)) <<" of "<< size_t (double(total_len) * 0.5) << " bases).");
            WARN("In most cases it means that either read coverage is uneven or significant contamination is present - both of these two cases make plasmidSPAdes' results unreliable");
            WARN("However, that situation may still be OK if you expect to see large plasmids in your dataset, so plasmidSPAdes will continue to work");
        } else {
            INFO(size_t((1 - fraction) * 100) << "% of bases from long edges have coverage significantly different from median");
        }
        for (EdgeId e : graph.edges()) {
            if (!long_component_.count(e)) {
                CalculateComponentSize(e, graph);
            }
        }
        INFO("Connected components calculated");
    } else {
        median_long_edge_coverage = external_chromosome_coverage;
    }

    size_t deleted = 0;
    for (auto iter = graph.SmartEdgeBegin(); ! iter.IsEnd(); ++iter) {
        EdgeId e = *iter;
        if (graph.length(e) <= long_edge_bound)
            continue;

        if (graph.coverage(e) >= median_long_edge_coverage * (1 + coverage_limits) ||
            graph.coverage(e) <= median_long_edge_coverage * (1 - coverage_limits))
            continue;

        DEBUG("Considering long edge: id " << graph.int_id(e) << " length: " << graph.length(e) << " coverage: " << graph.coverage(e));
        if (long_component_.count(e) && 300000 > long_component_[e] && deadends_count_[e] == 0) {
            DEBUG("Edge " << graph.int_id(e) << " skipped - because of small nondeadend connected component of size " << long_component_[e]);
        } else {
            DEBUG("Edge " << graph.int_id(e) << " deleted");
            deleted += 1;
            graph.DeleteEdge(e);
        }
    }
    INFO("Deleted " << deleted <<" long edges");

    CompressAll(graph);
    return median_long_edge_coverage;
}

void ChromosomeRemover::CoverageFilter(double coverage_cutoff) {
    auto& graph = gp_.get_mutable<Graph>();
    size_t deleted = 0;
    for (auto it = graph.SmartEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        // Never remove cycles
        if (graph.EdgeEnd(e) == graph.EdgeStart(e))
            continue;

        if (math::ge(graph.coverage(e), coverage_cutoff))
            continue;

        DEBUG("Deleting edge " << graph.int_id(e) << " because of low coverage");
        graph.DeleteEdge(e);
        deleted += 1;
    }
    INFO("With coverage cutoff " << coverage_cutoff << " removed " << deleted << " edges");

    CompressAll(graph);
}

void ChromosomeRemover::PlasmidSimplify(size_t long_edge_bound,
                                        std::function<void (EdgeId)> removal_handler ) {
    DEBUG("Simplifying graph for extrachromosomal removal project");
    auto& graph = gp_.get_mutable<Graph>();
    size_t iteration_count = 10;
    const auto &forbidden = gp_.get<SmartVertexSet>("forbidden_vertices");
    size_t forbidden_size = forbidden.size();
    INFO("Blocked vertices: " <<  forbidden_size);
    for (size_t i = 0; i < iteration_count; i++) {
        omnigraph::EdgeRemovingAlgorithm<Graph>
                tc(graph,
                   func::And(func::And(DeadEndCondition<Graph>(graph),
                                       LengthUpperBound<Graph>(graph, long_edge_bound)),
                             IsAllowedCondition<Graph>(graph, forbidden)),
                   removal_handler, true);
        tc.Run();
    }
    gp_.EnsureIndex();
}

//Debug only
void ChromosomeRemover::ReferenceBasedRemoveChromosomal() {
    auto& graph = gp_.get_mutable<Graph>();

    EdgesPositionHandler<Graph> edge_pos(graph, 0, 0);
    visualization::position_filler::PosFiller pos_filler(graph, MapperInstance(gp_), edge_pos);

    for (const auto &chr: gp_.get<GenomeStorage>().GetChromosomes()) {
        pos_filler.Process(chr.sequence, "0_" + chr.name);
        pos_filler.Process(ReverseComplement(chr.sequence), "1_" + chr.name);
    }

    size_t deleted = 0;
    for (auto iter = graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
        if (edge_pos.GetEdgePositions(*iter).size() == 0) {
            deleted++;
            graph.DeleteEdge(*iter);
        }
        DEBUG(*iter << " " << edge_pos.GetEdgePositions(*iter).size());
    }

    CompressAll(graph);
    INFO("Deleted with reference: " << deleted <<" edges");
}

void ChromosomeRemover::RemoveNearlyEverythingByCoverage(double cur_limit) {
    CoverageFilter(cur_limit);
    PlasmidSimplify(plasmid_config_.long_edge_length);
}

vector<vector<EdgeId>> ChromosomeRemover::GetNineShapeComponents () {
    const auto& graph = gp_.get<Graph>();
    long_vertex_component_.clear();
    long_component_.clear();
    deadends_count_.clear();
    component_list_.clear();
    vector<vector<EdgeId>> res;

    for (EdgeId e : graph.edges()) {
        if (long_component_.count(e) == 0) {
            CalculateComponentSize(e, graph);
        }
    }
    size_t count = 0;
    for (auto &comp: component_list_) {
        if (comp.size() == 4) {
            EdgeId first_edge = comp[0];
//conjugate, so /2
            size_t comp_size = (long_component_[first_edge])/2;
            size_t deadends_count = deadends_count_[first_edge] ;
            if (deadends_count != 2)
                break;
            int incoming = -1;
            for (size_t i = 0; i < comp.size(); i++) {
                if (graph.IsDeadStart(graph.EdgeStart(comp[i])) && (double) graph.length(comp[i]) < 0.3 * (double) comp_size) {
                    incoming = (int) i;
                    break;
                }
            }
            if (incoming == -1)
                break;
            int next_circular = -1;
            for (size_t i = 0; i < comp.size(); i++) {
                if (graph.EdgeStart(comp[i]) == graph.EdgeEnd(comp[i]) && graph.EdgeStart(comp[i]) == graph.EdgeEnd(comp[incoming])) {
                    next_circular = (int) i;
                    break;
                }
            }
            if (next_circular == -1)
                break;
            res.push_back(vector<EdgeId>{comp[incoming], comp[next_circular]});
            count ++;
        }
    }
    return res;
}

void ChromosomeRemover::OutputSuspiciousComponents () {
    auto& graph = gp_.get_mutable<Graph>();
    long_vertex_component_.clear();
    long_component_.clear();
    deadends_count_.clear();
    component_list_.clear();
    size_t component_size_max = 200000;
    size_t component_size_min = 1000;
    std::string tmp = std::to_string(ext_limit_);
    while (tmp.length() < 4) tmp = "_" + tmp;
    std::string out_file = "components" + tmp + ".fasta";
    double var = 0.3;
    DEBUG("calculating component sizes");
    for (EdgeId e: graph.canonical_edges()) {
        if (!long_component_.count(e)) {
            CalculateComponentSize(e, graph);
        }
    }
    CoverageUniformityAnalyzer coverage_analyzer(graph, 0);
    std::ofstream is(cfg::get().output_dir + out_file);
    size_t component_count = 1;
    const auto& used_edges = gp_.get<SmartContainer<std::unordered_set<EdgeId>, Graph>>("used_edges");
    for (auto &comp: component_list_) {
        VERIFY(comp.size() > 0);
        EdgeId first_edge = comp[0];
//conjugate, so /2
        size_t comp_size = (long_component_[first_edge])/2;
        size_t deadends_count = deadends_count_[first_edge] ;
        if (comp_size > component_size_min && comp_size < component_size_max &&
            (deadends_count <= 4)) {
            DEBUG("Checking component size " << comp_size);
            std::vector<std::pair<double, size_t>> coverages;
            size_t total_len = 0;
            size_t used_len = 0;
            for (auto edge : comp) {
                coverages.emplace_back(graph.coverage(edge), graph.length(edge));
                total_len += graph.length(edge);

                if (used_edges.count(edge) > 0)
                    used_len += graph.length(edge);
            }
            if (2 * used_len > total_len) {
                DEBUG("Already found circular path ");
                continue;
            }

            double average_cov = coverage_analyzer.CountMedianCoverage(coverages, total_len);
            size_t good_len = 0;
            for (auto edge : comp) {
                if (graph.coverage(edge) > (1 - var) * average_cov &&
                    graph.coverage(edge) < (1 + var) * average_cov)
                    good_len += graph.length(edge);
            }

            if (math::ls (average_cov, (double) ext_limit_ * 1.3)) {
                DEBUG("component coverage close to current limit");
            } else if (math::ls((double) good_len, 0.8 * (double) total_len)) {
                DEBUG("component coverage too variable: fraction close to average " << (double) good_len / (double) total_len);
            } else {
                DEBUG("Component is good!");
                size_t count = 1;
                for (auto edge : comp) {
                    if (!(edge <= graph.conjugate(edge)))
                        continue;

                    is << ">CUTOFF_" << ext_limit_ << "_COMPONENT_" << component_count << "_EDGE_" << count <<
                            "_length_" << graph.length(edge) << "_cov_" << graph.coverage(edge) << "_id_" << edge.int_id() << std::endl;
                    is << graph.EdgeNucls(edge) << std::endl;
                    count++;
//                        gp_.g.DeleteEdge(edge);
                }
                component_count += 1;
            }

        }
    }
}

void ChromosomeRemover::RunMetaPipeline() {
    using namespace omnigraph;
    using namespace omnigraph::de;
    const auto& graph = gp_.get<Graph>();
    using UsedEdges = omnigraph::SmartContainer<std::unordered_set<EdgeId>, Graph>;
    if (!gp_.count<UsedEdges>("used_edges"))
        gp_.add("used_edges", UsedEdges(graph));
    if (plasmid_config_.reference_removal != "") {
        VERIFY_MSG(false, "Reference-based chromosome removal is switched off");
        INFO("Removing all edges with no genomic sequence");
        ReferenceBasedRemoveChromosomal();
        return;
    }

    if (gp_.count<PairedInfoIndicesHandlerT<Graph>>("paired_handlers") == 0) {
        PairedInfoIndicesHandlerT<Graph> paired_handlers;
        PairedInfoIndicesHandlerT<Graph> scaffolding_handlers;
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            paired_handlers.emplace_back(PairedInfoIndexHandlerT<Graph>(gp_.get_mutable<PairedInfoIndicesT<Graph>>("clustered_indices")[i]));
            scaffolding_handlers.emplace_back(PairedInfoIndexHandlerT<Graph>(gp_.get_mutable<PairedInfoIndicesT<Graph>>("scaffolding_indices")[i]));
        }
        gp_.add("paired_handlers", paired_handlers);
        gp_.add("scaffolding_handlers", scaffolding_handlers);
        OutputEdgeSequences(graph, cfg::get().output_dir + "before_chromosome_removal");
    }
//first iteration of coverage-based chromosome removal
    if (gp_.count<SmartVertexSet>("forbidden_vertices") == 0) {
        INFO("Forbidding tip ends.. ");
        VertexSet forb;
        FillForbiddenSet(gp_.get_mutable<Graph>(), forb);
        gp_.add("forbidden_vertices", make_smart_container<VertexSet, Graph>(graph, forb));
    }
    size_t forbidden_size = gp_.get<SmartVertexSet>("forbidden_vertices").size();
    INFO("Forbidden (initial tip ends) vertex size: " << forbidden_size);
    OutputSuspiciousComponents ();
    std::string suffix = std::to_string(ext_limit_);
    while (suffix.length() < 4) suffix = "_" + suffix;

    OutputEdgesByID(graph, cfg::get().output_dir + "edges_before" + suffix);
    RemoveNearlyEverythingByCoverage((double) ext_limit_);
    FilterSmallComponents(); 
//Graph is not changed after this line and before next chromosome remover iteration

    if (plasmid_config_.output_linear) {
        auto nine_components = GetNineShapeComponents();
        if (gp_.count<PathContainer>("Plasmid paths")){
            gp_.get_mutable<PathContainer>("Plasmid paths").clear();
        } else {
            gp_.emplace_with_key<PathContainer>("Plasmid paths");
        }
        auto &pathContainer = gp_.get_mutable<PathContainer>("Plasmid paths");
        for (auto path: nine_components){
            pathContainer.CreatePair(graph, path);
        }
    }

}

void ChromosomeRemover::RunIsolatedPipeline() {
    //SmartContainer<std::unordered_set<VertexId>, Graph> forb(graph);
//currently not forbidding anything...
    VertexSet forb;
    const auto& graph = gp_.get<Graph>();
    gp_.add("forbidden_vertices", make_smart_container<VertexSet, Graph>(graph, forb));
    chromosome_coverage_ = RemoveLongGenomicEdges(plasmid_config_.long_edge_length,
                                                  plasmid_config_.relative_coverage);
    PlasmidSimplify(plasmid_config_.long_edge_length);
    for (size_t i = 0; i < max_iteration_count; i++) {
        size_t graph_size = graph.size();
        RemoveLongGenomicEdges(plasmid_config_.long_edge_length, plasmid_config_.relative_coverage,
                               chromosome_coverage_);
        INFO("Before dead_end simplification " << i << " " << graph.size() << " vertices in graph");
        PlasmidSimplify(plasmid_config_.long_edge_length);
        size_t new_graph_size = graph.size();
        if (new_graph_size == graph_size) {
            INFO("At iteration " << i << " graph was not changed");
            INFO(new_graph_size << " vertices left");
            break;
        }
    }
    FilterSmallComponents();
}

void ChromosomeRemover::FilterSmallComponents() {
    auto& graph = gp_.get_mutable<Graph>();
    //Small repetitive components after filtering
    std::unordered_map<VertexId, size_t> old_vertex_weights(long_vertex_component_.begin(), long_vertex_component_.end());
    for (size_t i = 0; i < max_iteration_count; i++) {
        DEBUG("Iteration " << i);
        size_t graph_size = graph.size();
        long_vertex_component_.clear();
        long_component_.clear();
        deadends_count_.clear();
        DEBUG("Calculating component sizes");
        for (EdgeId e: graph.canonical_edges()) {
            if (!long_component_.count(e)) {
                CalculateComponentSize(e, graph);
            }
        }
        DEBUG("Component sizes calculated");
//removing edges of coverage ~chromosome coverage that before this iteration were in relatively large components and now are in relatively small ones - both isolated and small components.
        for (auto iter = graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            EdgeId e = *iter;
            if (long_component_[e] >= 2 * plasmid_config_.small_component_size)
                continue;

            if (graph.IsDeadEnd(graph.EdgeEnd(e)) && graph.IsDeadStart(graph.EdgeStart(e)) &&
                old_vertex_weights.count(graph.EdgeStart(e)) &&
                // * 2 - because all coverages are taken with rc
                old_vertex_weights[graph.EdgeStart(e)] > long_component_[e] + plasmid_config_.long_edge_length * 2)  {
                DEBUG("Deleting isolated edge of length" << graph.length(e));
                graph.DeleteEdge(e);
            }
        }
        DEBUG("isolated deleted");
        for (auto iter = graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            EdgeId e = *iter;
            if (long_component_[e] >= 2 * plasmid_config_.small_component_size)
                continue;

            if (old_vertex_weights.count(graph.EdgeStart(e)) &&
                old_vertex_weights[graph.EdgeStart(e)] > plasmid_config_.small_component_size * 4 &&
                graph.coverage(e) < chromosome_coverage_ * (1 + plasmid_config_.small_component_relative_coverage) &&
                graph.coverage(e) > chromosome_coverage_ * (1 - plasmid_config_.small_component_relative_coverage)) {
                DEBUG("Deleting edge from fake small component, length " << graph.length(e) << " id " << graph.int_id(e) << " coverage " << graph.coverage(e) << " component_size " << old_vertex_weights[graph.EdgeStart(e)]);
                graph.DeleteEdge(e);
            }
        }
        DEBUG("components deleted");

// Small components with dead-ends and relatively short edges.
// TODO:: think, whether it may be bad in viral setting.
        for (auto iter = graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            EdgeId e = *iter;
            bool should_leave = deadends_count_[e] == 0;
            should_leave &= graph.length(e) > plasmid_config_.min_isolated_length;
            if (long_component_[e] < 2 * plasmid_config_.min_component_length &&
                !should_leave) {
                graph.DeleteEdge(e);
            }
        }

        CompressAll(graph);
        PlasmidSimplify(plasmid_config_.long_edge_length);

        size_t new_graph_size = graph.size();
        if (new_graph_size == graph_size) {
            INFO("At iteration " << i << " of small components additional filtering graph was not changed");
            if (new_graph_size == 0) {
                INFO("No putative extrachromosomal sequence remained!");
            } else {
                INFO("After chromosome removal subroutine " << new_graph_size << " vertices left");
            }
            break;
        }
    }
}



} //debruijn_graph
