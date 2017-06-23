#pragma once
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include "assembly_graph/graph_support/scaff_supplementary.hpp"

class RawContractedGraphBuilder {
    const Graph& g_;
    const size_t min_length_;
    path_extend::ScaffoldingUniqueEdgeStorage unique_storage_;

public:

    RawContractedGraphBuilder(const Graph& g, const size_t min_length) : g_(g), min_length_(min_length) {}

    void BuildRawContractedGraph() {
        typedef std::map<VertexId,std::size_t> rank_t; // => order on Element
        typedef std::map<VertexId, VertexId> parent_t;
        rank_t rank_map;
        parent_t parent_map;
        boost::associative_property_map<rank_t>   rank_pmap(rank_map);
        boost::associative_property_map<parent_t> parent_pmap(parent_map);
        omnigraph::IterationHelper <Graph, VertexId> vertex_iteration_helper(g_);

        BuildDSU(g_, rank_pmap, parent_pmap, min_length_);
    }

private:

    template <typename Rank, typename Parent>
            void BuildDSU(const Graph& g, Rank& r, Parent& p, const size_t min_length) {
        omnigraph::IterationHelper <Graph, EdgeId> edge_iteration_helper(g);

        boost::disjoint_sets<Rank,Parent> dsets(r, p);
        for (auto it = g.vertices().begin(); it != g.vertices().end(); ++it) {
            dsets.make_set(*it);
        }
        INFO("Overall vertices: " << dsets.count_sets(g.vertices().begin(), g.vertices().end()));
        size_t short_edges = 0;
        size_t self_linkages = 0;
        size_t contracted_long = 0;
        unordered_set<VertexId> long_roots;
        for (auto it = edge_iteration_helper.begin(); it != edge_iteration_helper.end(); ++it) {
            auto edge = *it;
            if (g.length(edge) < min_length) {
                ++short_edges;
                VertexId start = g.EdgeStart(edge);
                VertexId end = g.EdgeEnd(edge);
                if (dsets.find_set(start) == dsets.find_set(end)) {
                    ++self_linkages;
                }
                dsets.union_set(g.EdgeStart(edge), g.EdgeEnd(edge));
            } else {
                if (dsets.find_set(g.EdgeStart(edge)) == dsets.find_set(g.EdgeEnd(edge))) {
                    contracted_long++;
                }
                long_roots.insert(dsets.find_set(g.EdgeEnd(edge)));
                long_roots.insert(dsets.find_set(g.EdgeStart(edge)));
            }
        }
        INFO("Short edges: " << short_edges);
        INFO("Self linkages: " << self_linkages);
        INFO("Vertices in contracted: " << dsets.count_sets(g.vertices().begin(), g.vertices().end()));
        INFO("Contracted long: " << contracted_long);
        INFO("Long vertices: " << long_roots.size());
    };
};

//struct ReachableInfo {
//    EdgeId edge_;
//    size_t distance_;
//};
//
//class ContractedGraph {
//    std::unordered_map<EdgeId, vector<ReachableInfo>> edge_to_info_;
//
//
//    ContractedGraph()
//};