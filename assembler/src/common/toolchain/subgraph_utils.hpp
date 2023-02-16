//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "io/graph/gfa_writer.hpp"

#include <fstream>

#pragma once

namespace toolchain {

class ComponentExpander {
    const Graph &g_;

    bool IsInnerVertex(VertexId v, const std::unordered_set<EdgeId> &edges) const {
        auto in_f = [&edges](EdgeId e) {
            return edges.count(e);
        };
        return std::any_of(g_.in_begin(v), g_.in_end(v), in_f) &&
               std::any_of(g_.out_begin(v), g_.out_end(v), in_f);
    }

public:
    explicit ComponentExpander(const Graph &g) : g_(g) {}

    omnigraph::GraphComponent<Graph> Expand(const omnigraph::GraphComponent<Graph> &component) const {
        INFO("Expanding component to include incident edges of all 'inner' vertices");
        std::unordered_set<EdgeId> expanded_edges(component.edges());
        for (VertexId v : component.vertices()) {
            if (IsInnerVertex(v, component.edges())) {
                utils::insert_all(expanded_edges, g_.IncidentEdges(v));
            }
        }
        return omnigraph::GraphComponent<Graph>::FromEdges(g_,
                                                           expanded_edges.begin(),
                                                           expanded_edges.end());
    }

private:
    DECL_LOGGER("ComponentExpander");
};

static bool IsDeadEnd(const Graph &g, VertexId v) {
    return g.IncomingEdgeCount(v) * g.OutgoingEdgeCount(v) == 0;
}

static std::vector<EdgeId> CollectDeadEnds(const omnigraph::GraphComponent<Graph> &gc, bool canonical_only = true) {
    const auto &g = gc.g();
    std::vector<EdgeId> answer;
    for (EdgeId e : gc.edges()) {
        VERIFY(gc.edges().count(g.conjugate(e)));
        if (canonical_only && g.conjugate(e) < e)
            continue;
        if (IsDeadEnd(g, g.EdgeStart(e)) || IsDeadEnd(g, g.EdgeEnd(e))) {
            answer.push_back(e);
        }
    }
    return answer;
}

static void WriteComponentWithDeadends(const omnigraph::GraphComponent<Graph> &component,
                                       const std::string &prefix,
                                       const io::EdgeNamingF<Graph> &naming_f) {
    INFO("Writing GFA to " << prefix << ".gfa")
    std::ofstream os(prefix + ".gfa");
    const auto &g = component.g();
    gfa::GFAWriter writer(g, os, naming_f);
    writer.WriteSegmentsAndLinks(component);

    INFO("Writing dead-ends to " << prefix << ".deadends")
    std::ofstream deadends_os(prefix + ".deadends");
    for (EdgeId e : CollectDeadEnds(component)) {
        deadends_os << naming_f(g, e) << "\n";
    }
}

}