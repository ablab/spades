//***************************************************************************
//* Copyright (c) 2015-2019 Saint-Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "assembly_graph/graph_support/contig_output.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "chromosome_remover.hpp"
#include "common/paired_info/paired_info.hpp"
#include "pipeline/config_struct.hpp"
#include "assembly_graph/graph_support/coverage_uniformity_analyzer.hpp"

#include "utils/filesystem/path_helper.hpp"
#include "math/xmath.h"

namespace debruijn_graph {
using namespace std;

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

void ChromosomeRemover::FillForbiddenSet(Graph &g, unordered_set<VertexId> &forbidden) {
    for (VertexId v: g.vertices()) {
        if (g.IsDeadStart(v) || g.IsDeadEnd(v)) {
            forbidden.insert(v);
        }
    }
}

size_t ChromosomeRemover::CalculateComponentSize(EdgeId e, Graph &g) {
    std::unordered_set<EdgeId> next;
    size_t deadend_count = 0;
    next.insert(e);
    std::unordered_set<EdgeId> used;
    size_t ans = 0;
    while (!next.empty()){
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
    for (auto edge: used) {
        long_component_[edge] = ans;
        long_vertex_component_[g.EdgeStart(edge)] = ans;
        long_vertex_component_[g.EdgeEnd(edge)] = ans;
        deadends_count_[edge] = deadend_count;
    }
    component_list_.push_back(std::vector<EdgeId>(used.begin(), used.end()));
    return ans;
}

double ChromosomeRemover::RemoveLongGenomicEdges(size_t long_edge_bound, double coverage_limits, double external_chromosome_coverage){
    INFO("Removing of long chromosomal edges started");
    CoverageUniformityAnalyzer coverage_analyzer(gp_.g, long_edge_bound);
    size_t total_len = coverage_analyzer.TotalLongEdgeLength();
    if (total_len == 0) {
        if (external_chromosome_coverage < 1.0) {
            WARN("Plasmid detection failed, not enough long edges");
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
        for (auto iter = gp_.g.ConstEdgeBegin(); ! iter.IsEnd(); ++iter) {
            if (long_component_.find(*iter) == long_component_.end()) {
                CalculateComponentSize(*iter, gp_.g);
            }
        }
        INFO("Connected components calculated");
    } else {
        median_long_edge_coverage = external_chromosome_coverage;
    }
    size_t deleted = 0;
    for (auto iter = gp_.g.SmartEdgeBegin(); ! iter.IsEnd(); ++iter){
        if (gp_.g.length(*iter) > long_edge_bound) {
            if (gp_.g.coverage(*iter) < median_long_edge_coverage * (1 + coverage_limits) && gp_.g.coverage(*iter)  > median_long_edge_coverage * (1 - coverage_limits)) {
                DEBUG("Considering long edge: id " << gp_.g.int_id(*iter) << " length: " << gp_.g.length(*iter) <<" coverage: " << gp_.g.coverage(*iter));
                if ( long_component_.find(*iter) != long_component_.end() && 300000 > long_component_[*iter] && deadends_count_[*iter] == 0) {
                    DEBUG("Edge " << gp_.g.int_id(*iter) << " skipped - because of small nondeadend connected component of size " << long_component_[*iter]);
                } else {
                    DEBUG("Edge " << gp_.g.int_id(*iter) << "  deleted");
                    deleted++;
                    gp_.g.DeleteEdge(*iter);
                }
            }
        }
    }
    INFO("Deleted " << deleted <<" long edges");
    CompressAll(gp_.g);
    return median_long_edge_coverage;
}

void ChromosomeRemover::CoverageFilter(double coverage_cutoff) {
    size_t deleted = 0;
    for (auto it = gp_.g.SmartEdgeBegin(true); !it.IsEnd(); ++it) {
        if (math::ls(gp_.g.coverage(*it), coverage_cutoff) && gp_.g.EdgeEnd(*it) != gp_.g.EdgeStart(*it)) {
            DEBUG("Deleting edge " << gp_.g.int_id(*it) << " because of low coverage");
            gp_.g.DeleteEdge(*it);
            deleted ++;
        }
    }
    INFO("With limit " <<coverage_cutoff << "deleted  " << deleted <<" edges");
    CompressAll(gp_.g);
}

void ChromosomeRemover::PlasmidSimplify(size_t long_edge_bound,
                                        std::function<void (EdgeId)> removal_handler ) {
    DEBUG("Simplifying graph for plasmid project");
    size_t iteration_count = 10;
    INFO("Blocked vertices: " <<  gp_.get_const<std::unordered_set<VertexId>>("forbidden_vertices").size());
    for (size_t i = 0; i < iteration_count; i++) {
        omnigraph::EdgeRemovingAlgorithm<Graph> tc(gp_.g, func::And(func::And(DeadEndCondition<Graph>(gp_.g),
                LengthUpperBound<Graph>(gp_.g, long_edge_bound)), IsAllowedCondition<Graph>(gp_.g,  gp_.get_const<std::unordered_set<VertexId>>("forbidden_vertices"))),
                                                   removal_handler, true);
        tc.Run();
    }
    gp_.EnsureIndex();
}

//Debug only
void ChromosomeRemover::ReferenceBasedRemoveChromosomal() {
    EdgesPositionHandler<Graph> edge_pos(gp_.g, 0, 0);
    visualization::position_filler::PosFiller<Graph> pos_filler(gp_.g, MapperInstance(gp_), edge_pos);

    for (const auto &chr: gp_.genome.GetChromosomes()) {
        pos_filler.Process(chr.sequence, "0_" + chr.name);
        pos_filler.Process(ReverseComplement(chr.sequence), "1_" + chr.name);
    }
    size_t deleted = 0;
    for (auto iter = gp_.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {

        if (edge_pos.GetEdgePositions(*iter).size() == 0) {
            deleted++;
            gp_.g.DeleteEdge(*iter);
        }
        DEBUG(*iter << " " << edge_pos.GetEdgePositions(*iter).size());
    }
    CompressAll(gp_.g);
    INFO("Deleted with reference: " << deleted <<" edges");
}

void ChromosomeRemover::RemoveNearlyEverythingByCoverage(double cur_limit) {
    CoverageFilter(cur_limit);
    PlasmidSimplify(plasmid_config_.long_edge_length);
}

void ChromosomeRemover::OutputNineComponents (conj_graph_pack &gp, size_t ext_limit_) {
    long_vertex_component_.clear();
    long_component_.clear();
    deadends_count_.clear();
    component_list_.clear();
    string tmp = std::to_string(ext_limit_);
    while (tmp.length() < 4) tmp = "_" + tmp;
    std::string out_file = "final_contigs" + tmp + ".linear_repeat.fasta";
    std::ofstream is(cfg::get().output_dir + out_file);

    for (auto iter = gp.g.ConstEdgeBegin(true); !iter.IsEnd(); ++iter) {
        if (long_component_.find(*iter) == long_component_.end()) {
            CalculateComponentSize(*iter, gp.g);
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
                if (gp.g.IsDeadStart(gp.g.EdgeStart(comp[i])) && gp.g.length(comp[i]) < 0.3 * comp_size) {
                    incoming = i;
                    break;
                }
            }
            if (incoming == -1)
                break;
            int next_circular = -1;
            for (size_t i = 0; i < comp.size(); i++) {
                if (gp.g.EdgeStart(comp[i]) == gp.g.EdgeEnd(comp[i]) && gp.g.EdgeStart(comp[i]) == gp.g.EdgeEnd(comp[incoming])) {
                    next_circular = i;
                    break;
                }
            }
            if (next_circular == -1)
                break;
            stringstream ss;
            ss << gp.g.EdgeNucls(comp[incoming]);
            ss << gp.g.EdgeNucls(comp[next_circular]).Subseq(gp.g.k());
            string seq = ss.str();
            double cov = (gp.g.coverage(comp[incoming]) * gp.g.length(comp[incoming]) + gp.g.coverage(comp[next_circular]) * gp.g.length(comp[next_circular]))/(gp.g.length(comp[incoming]) + gp.g.length(comp[next_circular]));
            is << ">CUTOFF_" << ext_limit_ <<"_NINE_" << count <<
               "_length_"<< seq.length() <<"_cov_" << cov << "_id_" <<comp[incoming].int_id() << "_" << comp[next_circular].int_id() << endl;
            is <<seq << endl;
        }
    }
}

void ChromosomeRemover::OutputSuspiciousComponents () {
    long_vertex_component_.clear();
    long_component_.clear();
    deadends_count_.clear();
    component_list_.clear();
    size_t component_size_max = 200000;
    size_t component_size_min = 1000;
    string tmp = std::to_string(ext_limit_);
    while (tmp.length() < 4) tmp = "_" + tmp;
    std::string out_file = "components" + tmp + ".fasta";
    double var = 0.3;
    DEBUG("calculating component sizes");
    for (EdgeId e: gp_.g.canonical_edges()) {
        if (!long_component_.count(e)) {
            CalculateComponentSize(e, gp_.g);
        }
    }
    CoverageUniformityAnalyzer coverage_analyzer(gp_.g, 0);
    std::ofstream is(cfg::get().output_dir + out_file);
    size_t component_count = 1;
    VERIFY(gp_.count<unordered_set<EdgeId>>("used_edges") > 0);
    unordered_set<EdgeId> used_edges = gp_.get_const<unordered_set<EdgeId>>("used_edges");
    for (auto &comp: component_list_) {
        VERIFY(comp.size() > 0);
        EdgeId first_edge = comp[0];
//conjugate, so /2
        size_t comp_size = (long_component_[first_edge])/2;
        size_t deadends_count = deadends_count_[first_edge] ;
        if (comp_size > component_size_min && comp_size < component_size_max &&
            ( deadends_count <= 4)) {
            DEBUG("Checking component size " << comp_size);
            vector<pair<double, size_t>> coverages;
            size_t total_len = 0;
            size_t used_len = 0;
            for (auto edge:comp) {
                coverages.push_back (make_pair(gp_.g.coverage(edge), gp_.g.length(edge)));
                total_len += gp_.g.length(edge);

                if (used_edges.count(edge) > 0)
                    used_len += gp_.g.length(edge);
            }
            if (2 * used_len > total_len) {
                DEBUG("Already found circular path ");
                continue;
            }
            double average_cov = coverage_analyzer.CountMedianCoverage(coverages, total_len);
            size_t good_len = 0;
            for (auto edge:comp) {
                if (gp_.g.coverage(edge) > (1 - var) * average_cov && gp_.g.coverage(edge) < (1 + var) * average_cov)
                    good_len += gp_.g.length(edge);
            }
            if (math::ls (average_cov, (double) ext_limit_ * 1.3)) {
                DEBUG ("component coverage close to current limit");
            } else if (math::ls((double) good_len, 0.8 * (double) total_len)) {
                DEBUG ("component coverage too variable: fraction close to average" << (double) good_len / (double) total_len);
            } else {
                DEBUG("Component is good!");
                size_t count = 1;
                for (auto edge: comp) {
                    if (edge <= gp_.g.conjugate(edge)) {
                        is << ">CUTOFF_" << ext_limit_ <<"_COMPONENT_" << component_count << "_EDGE_" << count <<
                           "_length_"<< gp_.g.length(edge) <<"_cov_" << gp_.g.coverage(edge) << "_id_" <<edge.int_id() << endl;
                        is << gp_.g.EdgeNucls(edge) << endl;
                        count++;
//                        gp_.g.DeleteEdge(edge);
                    }
                }

                component_count ++;
            }

        }
    }
}

void ChromosomeRemover::RunMetaPipeline() {
    if (plasmid_config_.reference_removal != "") {
        VERIFY_MSG(false, "Reference-based chromosome removal is switched off");
        INFO("Removing all edges with no genomic sequence");
        ReferenceBasedRemoveChromosomal();
        return ;
    }

    if (gp_.count<omnigraph::de::PairedInfoIndicesHandlerT<Graph>>("paired_handlers") == 0) {
        omnigraph::de::PairedInfoIndicesHandlerT<Graph> paired_handlers;
        omnigraph::de::PairedInfoIndicesHandlerT<Graph> scaffolding_handlers;
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            paired_handlers.emplace_back(omnigraph::de::PairedInfoIndexHandlerT<Graph>(gp_.clustered_indices[i]));
            scaffolding_handlers.emplace_back(omnigraph::de::PairedInfoIndexHandlerT<Graph>(gp_.scaffolding_indices[i]));
        }
        gp_.add("paired_handlers", paired_handlers);
        gp_.add("scaffolding_handlers", scaffolding_handlers);
        OutputEdgeSequences(gp_.g, cfg::get().output_dir + "before_chromosome_removal");
    }
//first iteration of coverage-based chromosome removal
    if (gp_.count<std::unordered_set<VertexId>>("forbidden_vertices") == 0) {
        INFO("Forbidding tip ends.. ");
        unordered_set<VertexId> forb;
        FillForbiddenSet(gp_.g, forb);
        gp_.add("forbidden_vertices", forb);
    }
    INFO("Forbidden (initial tip ends) vertex size: " << gp_.get_const<std::unordered_set<VertexId>>("forbidden_vertices").size() );
    OutputSuspiciousComponents ();
    OutputNineComponents (gp_, ext_limit_);
    string tmp = std::to_string(ext_limit_);
    while (tmp.length() < 4) tmp = "_" + tmp;

    OutputEdgesByID(gp_.g, cfg::get().output_dir + "edges_before" + tmp);
    RemoveNearlyEverythingByCoverage((double) ext_limit_);
}

void ChromosomeRemover::RunIsolatedPipeline() {
    unordered_set<VertexId> forb;
//currently not forbidding anything...
    gp_.add("forbidden_vertices", forb);
    chromosome_coverage_ = RemoveLongGenomicEdges(plasmid_config_.long_edge_length,
                                                 plasmid_config_.relative_coverage);
    PlasmidSimplify(plasmid_config_.long_edge_length);
    for (size_t i = 0; i < max_iteration_count; i++) {
        size_t graph_size = gp_.g.size();
        RemoveLongGenomicEdges(plasmid_config_.long_edge_length, plasmid_config_.relative_coverage,
                               chromosome_coverage_);
        INFO("Before dead_end simplification " << i << " " << gp_.g.size() << " vertices in graph");
        PlasmidSimplify(plasmid_config_.long_edge_length);
        size_t new_graph_size = gp_.g.size();
        if (new_graph_size == graph_size) {
            INFO("Iteration " << i << " graph was not changed");
            INFO(new_graph_size << " vertices left");
            break;
        }
    }
}

void ChromosomeRemover::FilterSmallComponents() {
    //Small repetitive components after filtering
    std::unordered_map<VertexId, size_t> old_vertex_weights (long_vertex_component_.begin(), long_vertex_component_.end());
    for (size_t i = 0; i < max_iteration_count; i++) {
        DEBUG("Iteration " << i);
        size_t graph_size = gp_.g.size();
        long_vertex_component_.clear();
        long_component_.clear();
        deadends_count_.clear();
        DEBUG("Calculating component sizes");
        for (EdgeId e: gp_.g.canonical_edges()) {
            if (!long_component_.count(e)) {
                CalculateComponentSize(e, gp_.g);
            }
        }
        DEBUG("Component sizes calculated");
//removing edges of coverage ~chromosome coverage that before this iteration were in relatively large components and now are in relatively small ones - both isolated and small components.
        for (auto iter = gp_.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (gp_.g.IsDeadEnd(gp_.g.EdgeEnd(*iter)) && gp_.g.IsDeadStart(gp_.g.EdgeStart(*iter))
                && old_vertex_weights.find(gp_.g.EdgeStart(*iter)) !=  old_vertex_weights.end()
                // * 2 - because all coverages are taken with rc
                && old_vertex_weights[gp_.g.EdgeStart(*iter)] > long_component_[*iter] + plasmid_config_.long_edge_length * 2)  {
                DEBUG("Deleting isolated edge of length" << gp_.g.length(*iter));
                gp_.g.DeleteEdge(*iter);
            }
        }
        for (auto iter = gp_.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (long_component_[*iter] < 2 * plasmid_config_.small_component_size) {
                if (old_vertex_weights.find(gp_.g.EdgeStart(*iter)) != old_vertex_weights.end() &&
                    old_vertex_weights[gp_.g.EdgeStart(*iter)] >
                    plasmid_config_.small_component_size * 4 &&
                    gp_.g.coverage(*iter) < chromosome_coverage_ * (1 + plasmid_config_.small_component_relative_coverage)
                    && gp_.g.coverage(*iter) > chromosome_coverage_ * (1 - plasmid_config_.small_component_relative_coverage)) {
                    DEBUG("Deleting edge from fake small component, length " << gp_.g.length(*iter) << " id " << gp_.g.int_id(*iter) << " coverage " << gp_.g.coverage(*iter) << " component_size " << old_vertex_weights[gp_.g.EdgeStart(*iter)]) ;
                    gp_.g.DeleteEdge(*iter);
                }
            }
        }


// Small components with dead-ends and relatively short edges.
// TODO:: think, whether it may be bad in viral setting.
        for (auto iter = gp_.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            bool should_leave = deadends_count_[*iter] == 0;
            should_leave &= gp_.g.length(*iter) > plasmid_config_.min_isolated_length;
            if (long_component_[*iter] < 2 * plasmid_config_.min_component_length &&
                !should_leave) {
                gp_.g.DeleteEdge(*iter);
            }
        }
        CompressAll(gp_.g);
        PlasmidSimplify(plasmid_config_.long_edge_length);
        size_t new_graph_size = gp_.g.size();
        if (new_graph_size == graph_size) {
            INFO("Iteration " << i << " of small components additional filtering graph was not changed");
            if (new_graph_size == 0) {
                WARN("No putative plasmid contigs found!");
            } else {
                INFO("After chromosome removal subroutine " << new_graph_size << " vertices left");
            }
            break;
        }
    }
}



} //debruijn_graph
