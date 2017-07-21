//****************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "assembly_graph/graph_support/contig_output.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "modules/simplification/ec_threshold_finder.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "assembly_graph/handlers/edge_fate_position_handler.hpp"
#include "chromosome_removal.hpp"
#include "common/paired_info/paired_info.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "pipeline/config_struct.cpp"

#include "utils/filesystem/path_helper.hpp"
#include <libgen.h>
#include "math/xmath.h"

namespace debruijn_graph {
using namespace std;
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

void ChromosomeRemoval::FillForbiddenSet(Graph &g, std::set<VertexId> &forbidden_) {
    for (auto it = g.SmartVertexBegin(); ! it.IsEnd(); ++it) {
        if (g.IsDeadStart(*it) || g.IsDeadEnd(*it)) {
            forbidden_.insert(*it);
        }
    }
    INFO(forbidden_.size());
}

size_t ChromosomeRemoval::CalculateComponentSize(EdgeId e, Graph &g_) {
    std::unordered_set<EdgeId> next;
    size_t deadend_count = 0;
    next.insert(e);
    std::unordered_set<EdgeId> used;
    size_t ans = 0;
    while (!next.empty()){
        auto cur = *next.begin();
        next.erase(next.begin());
        if (used.find(cur) != used.end()) {
            continue;
        }
        ans += g_.length(cur);

        used.insert(cur);
        std::vector<EdgeId> neighbours;
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
                next.insert(ee);
            }
        }
    }
    for (auto edge: used) {
        long_component_[edge] = ans;
        long_vertex_component_[g_.EdgeStart(edge)] = ans;
        long_vertex_component_[g_.EdgeEnd(edge)] = ans;
        deadends_count_[edge] = deadend_count;
    }
    component_list_.push_back(std::vector<EdgeId>(used.begin(), used.end()));
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

void ChromosomeRemoval::CoverageFilter(conj_graph_pack &gp, double coverage_cutoff) {
    size_t deleted = 0;
    for (auto it = gp.g.SmartEdgeBegin(true); !it.IsEnd(); ++it) {
        if (gp.g.coverage(*it) < coverage_cutoff && gp.g.EdgeEnd(*it) != gp.g.EdgeStart(*it)) {
            DEBUG("Deleting edge " << gp.g.int_id(*it) << " because of low coverage");
            gp.g.DeleteEdge(*it);
            deleted ++;
        }
    }
    INFO("With limit " <<coverage_cutoff << "deleted  " << deleted <<" edges");
    CompressAll(gp.g);
}

void ChromosomeRemoval::PlasmidSimplify(conj_graph_pack &gp, size_t long_edge_bound,
                                        std::function<void (EdgeId)> removal_handler ) {
    if (!cfg::get().pd->circular_removal)
        return;
    DEBUG("Simplifying graph for plasmid project");
    size_t iteration_count = 10;
    INFO("blocked vertices: " << gp.forbidden_vertices.size());
    for (size_t i = 0; i < iteration_count; i++) {
        omnigraph::EdgeRemovingAlgorithm<Graph> tc(gp.g, func::And(func::And(DeadEndCondition<Graph>(gp.g), LengthUpperBound<Graph>(gp.g, long_edge_bound)), IsAllowedCondition<Graph>(gp.g, gp.forbidden_vertices)),
                                                   removal_handler, true);
        tc.Run();
    }
    gp.EnsureIndex();
}

void ChromosomeRemoval::ReferenceBasedRemoveChromosomal(conj_graph_pack &gp) {
    EdgesPositionHandler<Graph> edge_pos(gp.g, 0, 0);
    visualization::position_filler::PosFiller<Graph> pos_filler(gp.g, MapperInstance(gp), edge_pos);

    for (const auto &chr: gp.genome.GetChromosomes()) {
        pos_filler.Process(chr.sequence, "0_" + chr.name);
        pos_filler.Process(ReverseComplement(chr.sequence), "1_" + chr.name);
    }
    size_t deleted = 0;
    for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {

        if (edge_pos.GetEdgePositions(*iter).size() == 0) {
            deleted++;
            gp.g.DeleteEdge(*iter);
        }
        DEBUG(*iter << " " << edge_pos.GetEdgePositions(*iter).size());
    }
    CompressAll(gp.g);
    INFO("deleted with reference: " << deleted <<" edges");
}

void ChromosomeRemoval::RemoveNearlyEverythingByCoverage(conj_graph_pack &gp, size_t cur_limit) {
    CoverageFilter(gp, cur_limit);
    PlasmidSimplify(gp, cfg::get().pd->long_edge_length);
}

void ChromosomeRemoval::OutputSuspiciousComponents (conj_graph_pack &gp, size_t ext_limit_) {
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
    for (auto iter = gp.g.ConstEdgeBegin(true); ! iter.IsEnd(); ++iter) {
        if (long_component_.find(*iter) == long_component_.end()) {
            CalculateComponentSize(*iter, gp.g);
        }
    }
    CoverageUniformityAnalyzer coverage_analyzer(gp.g, 0);
    std::ofstream is(cfg::get().output_dir + out_file);
    size_t component_count = 1;
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
                coverages.push_back (make_pair(gp.g.coverage(edge), gp.g.length(edge)));
                total_len += gp.g.length(edge);
                if (gp.used_edges.find(edge) != gp.used_edges.end())
                    used_len += gp.g.length(edge);
            }
            if (used_len > 0.5 * total_len) {
                DEBUG("Already found circular path ");
                continue;
            }
            double average_cov = coverage_analyzer.CountMedianCoverage(coverages, total_len);
            size_t good_len = 0;
            for (auto edge:comp) {
                if (gp.g.coverage(edge) > (1 - var) * average_cov && gp.g.coverage(edge) < (1 + var) * average_cov)
                    good_len += gp.g.length(edge);
            }
            if (average_cov < ext_limit_ * 1.3) {
                DEBUG ("component coverage close to current limit");
            } else if (good_len < 0.8 * total_len) {
                DEBUG ("component coverage too variable: fraction close to average" << good_len *1.0/total_len);
            } else {
                DEBUG("Component is good!");
                size_t count = 1;
                for (auto edge: comp) {
                    if (edge <= gp.g.conjugate(edge)) {
                        is << ">CUTOFF_" << ext_limit_ <<"_COMPONENT_" << component_count << "_EDGE_" << count <<
                           "_length_"<< gp.g.length(edge) <<"_cov_" << gp.g.coverage(edge) << "_id_" <<edge.int_id() << endl;
                        is << gp.g.EdgeNucls(edge) << endl;
                        count++;
//                        gp.g.DeleteEdge(edge);
                    }
                }

                component_count ++;
            }

        }
    }
}


void ChromosomeRemoval::run(conj_graph_pack &gp, const char*) {
    if (gp.count<omnigraph::de::PairedInfoIndicesHandlerT<Graph>>("paired_handlers") == 0) {
        omnigraph::de::PairedInfoIndicesHandlerT<Graph> paired_handlers;
        omnigraph::de::PairedInfoIndicesHandlerT<Graph> scaffolding_handlers;
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            paired_handlers.emplace_back(omnigraph::de::PairedInfoIndexHandlerT<Graph>(gp.clustered_indices[i]));
            scaffolding_handlers.emplace_back(omnigraph::de::PairedInfoIndexHandlerT<Graph>(gp.scaffolding_indices[i]));
        }
        gp.add("paired_handlers", paired_handlers);
        gp.add("scaffolding_handlers", scaffolding_handlers);
        OutputEdgeSequences(gp.g, cfg::get().output_dir + "before_chromosome_removal");
    }
    INFO("Before iteration " << 0 << ", " << gp.g.size() << " vertices in graph");
    std::string additional_list = cfg::get().pd->remove_list;
    bool use_chromosomal_list = (additional_list != "" && cfg::get().pd->HMM_filtration != "none");


    if (cfg::get().pd->reference_removal != "") {
        VERIFY_MSG(false, "Reference-based chromosome removal is switched off");
        INFO("Removing all edges with no genomic sequence");
        ReferenceBasedRemoveChromosomal(gp);
        return ;
    }

    double chromosome_coverage;
    if (cfg::get().pd->meta_mode) {
        if (cfg::get().pd->iterative_coverage_elimination) {
//TODO::very WTF
            if (ext_limit_ == 5) {
                INFO("forbidding tip ends");
                FillForbiddenSet(gp.g, gp.forbidden_vertices);
                INFO("FORBIDDING SIZE: " << gp.forbidden_vertices.size());
            }
            OutputSuspiciousComponents (gp, ext_limit_);
            string tmp = std::to_string(ext_limit_);
            while (tmp.length() < 4) tmp = "_" + tmp;
            
            OutputEdgesByID(gp.g, cfg::get().output_dir + "edges_before" + tmp + ".fasta");
            RemoveNearlyEverythingByCoverage(gp, ext_limit_);
        }  else {
            VERIFY(false);
        }
    }
    else
        chromosome_coverage = RemoveLongGenomicEdges(gp, cfg::get().pd->long_edge_length,
                                                     cfg::get().pd->relative_coverage);

    PlasmidSimplify(gp, cfg::get().pd->long_edge_length);
//TODO:: reconsider and move somewhere(not to config)
    size_t max_iteration_count = 30;

    for (size_t i = 0; i < max_iteration_count; i++) {
        size_t graph_size = gp.g.size();
        if (!cfg::get().pd->meta_mode)
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
        DEBUG("Iteration " << i);
        size_t graph_size = gp.g.size();
        long_vertex_component_.clear();
        long_component_.clear();
        deadends_count_.clear();    
        DEBUG("calculating component sizes");
        for (auto iter = gp.g.ConstEdgeBegin(); ! iter.IsEnd(); ++iter) {
            if (long_component_.find(*iter) == long_component_.end()) {
                CalculateComponentSize(*iter, gp.g);
            }
        }
        DEBUG("component sizes calculated");

        for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (gp.g.IsDeadEnd(gp.g.EdgeEnd(*iter)) && gp.g.IsDeadStart(gp.g.EdgeStart(*iter))
                && old_vertex_weights.find(gp.g.EdgeStart(*iter)) !=  old_vertex_weights.end()
//* 2 - because all coverages are taken with rc
                && old_vertex_weights[gp.g.EdgeStart(*iter)] > long_component_[*iter] + cfg::get().pd->long_edge_length * 2)  {
                DEBUG("deleting isolated edge of length" << gp.g.length(*iter));
                gp.g.DeleteEdge(*iter);
            }
        }
        if (! use_chromosomal_list) {
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
        }
//      TODO: some similar ideas for list-based removal

        for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            bool should_leave = deadends_count_[*iter] == 0;
            if (!use_chromosomal_list)
                should_leave &= gp.g.length(*iter) > cfg::get().pd->min_isolated_length;
            if (long_component_[*iter] < 2 * cfg::get().pd->min_component_length &&
                                  !should_leave) {
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
