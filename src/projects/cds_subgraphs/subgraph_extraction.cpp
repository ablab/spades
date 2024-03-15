//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "subgraph_extraction.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/paths/path_utils.hpp"

namespace cds_subgraphs {

omnigraph::GraphComponent<Graph>
MinDistRelevantComponentFinder::RelevantComponent(size_t cds_len_est,
                                        size_t base_len,
                                        const DistInfo &starts,
                                        const DistInfo &ends) const {
    INFO("Extracting relevant component");
    if (starts.size() * ends.size() == 0) {
        ERROR("Set of stop codons to the left/right was empty");
        VERIFY(false);
    }
    std::unordered_map<VertexId, DistInfo> s_dists;
    std::unordered_map<VertexId, DistInfo> e_dists;
    const size_t max_cds_len = RoundedProduct(cds_len_est, max_len_coeff_);

    //TODO check Dijkstra return codes
    size_t global_bound = RoundedProduct(std::max(MaxDist(starts), MaxDist(ends)) + base_len, max_len_coeff_);
    DEBUG("Max starts dist " << MaxDist(starts));
    DEBUG("Max ends dist " << MaxDist(ends));
    DEBUG("Base len " << base_len);
    DEBUG("Global length bound " << global_bound);

    //TODO can optimize! First process direction with fewer vertices, then block unreached vertices
    //TODO we can ignore light edges
    const size_t MAX_VERTEX_BOUND = 10000;
    auto fwd_dijkstra = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, global_bound,
                                                                     MAX_VERTEX_BOUND);
    DEBUG("Launching forward Dijkstras from ends of " << starts.size() << " edges");
    //TODO can be slightly optimized for several starts on same edge
    for (const auto &p : utils::key_set(starts)) {
        fwd_dijkstra.Run(g_.EdgeEnd(p.first));
        //if (fwd_dijkstra.VertexLimitExceeded()) {
        //    WARN("Vertex limit was exceeded while launching Dijkstra "
        //            "from end of edge " << g_.str(p.first));
        //    return GraphComponent<Graph>(g_);
        //}
        for (VertexId v : fwd_dijkstra.ReachedVertices()) {
            s_dists[v][p] = fwd_dijkstra.GetDistance(v) + g_.length(p.first) - p.second;
        }
    }

    if (fwd_dijkstra.VertexLimitExceeded()) {
        DEBUG("Vertex limit was exceeded while launching forward Dijkstras");
    }

    DEBUG("Launching backward Dijkstras from starts of " << ends.size() << " edges");
    //TODO deduplicate?!!
    auto bwd_dijkstra = omnigraph::DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(g_, global_bound,
                                                                             MAX_VERTEX_BOUND);
    for (const auto &p : utils::key_set(ends)) {
        bwd_dijkstra.Run(g_.EdgeStart(p.first));
        //if (bwd_dijkstra.VertexLimitExceeded()) {
        //    WARN("Vertex limit was exceeded while launching Dijkstra "
        //            "from start of edge " << g_.str(p.first));
        //    return GraphComponent<Graph>(g_);
        //}
        for (VertexId v : bwd_dijkstra.ReachedVertices()) {
            e_dists[v][p] = bwd_dijkstra.GetDistance(v) + p.second;
        }
    }

    if (bwd_dijkstra.VertexLimitExceeded()) {
        DEBUG("Vertex limit was exceeded while launching backward Dijkstras");
    }

    DEBUG("Intersecting results");
    //TODO correct usage of max_len_coeff to nucleote length
    BaseDistF base_dist_f =
            [this, &starts, &ends, base_len](GraphPos s, GraphPos e) {
                return RoundedProduct(utils::get(starts, s) + utils::get(ends, e) + base_len, max_len_coeff_);
            };

    //vertices with minimal distance to any end among valid paths
    std::vector<std::pair<size_t, VertexId>> dist_vertices;

    for (const auto &v_di: s_dists) {
        VertexId v = v_di.first;
        const auto &s_dist = v_di.second;

        auto it = e_dists.find(v);
        if (it == e_dists.end())
            continue;

        const auto &e_dist = it->second;
        size_t min_e_dist = MinExitDist(s_dist, e_dist, base_dist_f);
        if (min_e_dist < std::numeric_limits<size_t>::max()) {
            dist_vertices.emplace_back(min_e_dist, v);
        } else {
            TRACE("All paths that vertex " << g_.str(v) << " belongs to are much longer than the shortest path between the same start/end pairs");
        }
    }
    //Adding starts of the edges with stop codons backward into consideration
    for (const auto &p : utils::key_set(starts)) {
        dist_vertices.emplace_back(std::numeric_limits<size_t>::max(), g_.EdgeStart(p.first));
    }

    DEBUG("Creating component");
    std::sort(dist_vertices.begin(), dist_vertices.end());
    std::set<VertexId> within_cds_limit;

    //FIXME what if corresponding stop codon is too far?
    //Probably need to check stop codon distance first!

    //TODO add finding potential start codons?
    for (auto d_v : dist_vertices) {
        size_t min_end_dist = d_v.first;
        VertexId v = d_v.second;
        if (min_end_dist < max_cds_len) {
            within_cds_limit.insert(v);
        }
    }

    auto gc = omnigraph::GraphComponent<Graph>::FromVertices(g_, within_cds_limit.begin(), within_cds_limit.end());

    //Adding edges upstream
    std::set<EdgeId> revised_edges(gc.edges());
    for (auto d_v : dist_vertices) {
        size_t min_end_dist = d_v.first;
        VertexId v = d_v.second;
        if (min_end_dist >= max_cds_len) {
            for (EdgeId e : g_.OutgoingEdges(v)) {
                if (within_cds_limit.count(g_.EdgeEnd(e))) {
                    revised_edges.insert(e);
                }
             }
         }
     }

    //Adding edges containing stop codons forward (if length check passed)
    for (const auto &p_d : ends) {
        EdgeId e = p_d.first.first;
        if (base_len + p_d.second < max_cds_len) {
            if (!CheckConnectedTo(g_.EdgeStart(e), within_cds_limit)) {
                //FIXME add entire path?
                WARN("Edge " << edge_naming_f_(g_, e) << " being added as containing stop codon "
                                             "is disconnected from inner component vertices");
            }

            revised_edges.insert(e);
        } else {
            WARN("Edge " << edge_naming_f_(g_, e) << " containing stop codon "
                                         "on position " << p_d.first.second << " was not added to component");
            INFO("Max cds len " << max_cds_len << " ; base_len " << base_len << " ; dist to stop " << p_d.second);
        }
    }

    DEBUG("Relevant component extracted");
    return omnigraph::GraphComponent<Graph>::FromEdges(g_, revised_edges.begin(), revised_edges.end());
}

omnigraph::GraphComponent<Graph>
PartialGenePathProcessor::ProcessPartialGenePath(const EdgePath &gene_path,
                                                 size_t cds_len_est,
                                                 std::set<GraphPos> *stop_codon_poss) const {
    if (gene_path.size() == 0) {
        WARN("Partial gene path was empty");
        return omnigraph::GraphComponent<Graph>(g_);
    }
    INFO("Processing partial gene path " << PrintEdgePath(gene_path));
    //DEBUG("Nucls: " << PathSeq(gene_path));

    INFO("Searching for stop codon forward");
    EdgeId last_e = gene_path.sequence().back();
    CodonFinder fwd_finder(g_, GraphPos(last_e, gene_path.end_pos() - 1), STOP_CODONS);
    auto paths_fwd = fwd_finder.Go();
    INFO("Codons forward found: " << paths_fwd.size());
    auto ends = fwd_finder.Terminates(paths_fwd);

    INFO("Added to set of potential stop codons for the gene");
    for (const auto &gp_d : ends) {
        //NB: Inserted number gives the (k+1-mer) coordinate beyond the
        //    k+1-mer whose sequence ends with stop codon
        stop_codon_poss->insert(gp_d.first);
    }

    INFO("Paths forward:");
    for (const auto &p : paths_fwd) {
        INFO("Path: " << PrintEdgePath(p));
        Sequence s = PathSeq(p).Subseq(g_.k());
        TRACE("Sequence: len=" << s.size() << "; " << s);
    }

    INFO("Searching for stop codon backward");
    EdgeId first_e = gene_path.sequence().front();
    CodonFinder bwd_finder(g_, GraphPos(g_.conjugate(first_e), g_.length(first_e) - gene_path.start_pos() - 1),
                           RC_STOP_CODONS);
    auto paths_bwd = bwd_finder.Go();
    INFO("Codons backward found: " << paths_bwd.size());
    auto starts = RCTerminates(bwd_finder.Terminates(paths_bwd));

    INFO("Paths backward (reversed): ");
    for (const auto &p : RCEdgePaths(paths_bwd)) {
        INFO("Path: " << PrintEdgePath(p));
        Sequence s = PathSeq(p);
        TRACE("Sequence: len=" << s.size() - g_.k() << "; " << s.Subseq(0, s.size() - g_.k()));
    }

    INFO("Graph sinks to consider: " << ends.size());
    //TODO optimize debug logging
    for (const auto &gpos_d: ends) {
        EdgeId e = gpos_d.first.first;
        DEBUG("Edge " << edge_naming_f_(g_, e) << " pos " << gpos_d.first.second);
    }
    //TODO optimize debug logging
    INFO("Graph sources to consider: " << starts.size());
    for (const auto &gpos_d: starts) {
        EdgeId e = gpos_d.first.first;
        DEBUG("Edge " << edge_naming_f_(g_, e) << " pos " << gpos_d.first.second);
    }

    return (starts.size() * ends.size() == 0) ? omnigraph::GraphComponent<Graph>(g_) :
           rel_comp_finder_.RelevantComponent(cds_len_est, PathLength(g_, gene_path), starts, ends);
}

EdgePath
CDSSubgraphExtractor::ExtractPath(const std::string &partial_cds, const SequenceMapper<Graph> &mapper) const {
    io::SingleRead r("empty_name", partial_cds);
    //TODO output more info about particular prediction mapping!
    EdgePath path = path_finder_.FindDetailedReadPath(mapper.MapRead(r));
    if (!CheckContiguous(g_, path.sequence())) {
//            TRACE("Path for one of predictions for " << gene_id << " not contiguous");
        TRACE("Path not contiguous");
    }
    if (PathLength(g_, path) + g_.k() < RoundedProduct(r.size(), 0.7)) {
//            WARN("Path for one of predictions for " << gene_id << " was much shorter than query");
        WARN("Path much shorter than query");
        return EdgePath();
    } else {
        return path;
    }
}

omnigraph::GraphComponent<Graph>
CDSSubgraphExtractor::ProcessPartialCDS(const std::string &partial_cds,
                                        size_t cds_length_estimate,
                                        std::set<GraphPos> *stop_codon_poss,
                                        size_t min_len_to_explore,
                                        double frac_to_explore) const {
    CHECK_FATAL_ERROR(partial_cds.size() % 3 == 0, "Size of partial CDS prediction not divisible by 3");
    VERIFY(math::le(frac_to_explore, 1.));
    INFO("Original query length " << partial_cds.size());
    INFO("CDS length estimate " << cds_length_estimate);
    if (partial_cds.size() < min_len_to_explore) {
        WARN("Size of partial CDS prediction shorter than " << min_len_to_explore);
        min_len_to_explore = partial_cds.size();
    }
    size_t crop_len = std::max(min_len_to_explore, RoundedProduct(partial_cds.size(), frac_to_explore));
    //aligning on frame
    crop_len = crop_len / 3 * 3;
    INFO("Will explore left/right parts of length " << crop_len);

    std::set<EdgeId> edges;

    VERIFY(crop_len <= partial_cds.size());

    INFO("Processing rightmost part of the query of length " << crop_len);
    utils::insert_all(edges, cds_path_processor_.ProcessPartialGenePath(
                ExtractPath(partial_cds.substr(partial_cds.size() - crop_len), mapper_),
                cds_length_estimate, stop_codon_poss).edges());

    INFO("Processing leftmost part of the query of length " << crop_len);
    utils::insert_all(edges, cds_path_processor_.ProcessPartialGenePath(
                ExtractPath(partial_cds.substr(0, crop_len), mapper_),
                cds_length_estimate, stop_codon_poss).edges());

    if (edges.size() == 0) {
        INFO("Couldn't gather any reasonable component");
        return omnigraph::GraphComponent<Graph>(g_);
    }

    return omnigraph::GraphComponent<Graph>::FromEdges(g_, edges.begin(), edges.end());
}

}
