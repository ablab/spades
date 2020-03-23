//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "io/graph/gfa_writer.hpp"

#include <fstream>

namespace toolchain {

class ComponentExpander {
    const debruijn_graph::Graph &g_;

    bool IsInnerVertex(debruijn_graph::VertexId v, const std::set<debruijn_graph::EdgeId> &edges) const {
        auto in_f = [&edges](debruijn_graph::EdgeId e) {
            return edges.count(e);
        };
        return std::any_of(g_.in_begin(v), g_.in_end(v), in_f) &&
               std::any_of(g_.out_begin(v), g_.out_end(v), in_f);
    }

public:
    explicit ComponentExpander(const debruijn_graph::Graph &g) : g_(g) {}

    omnigraph::GraphComponent<debruijn_graph::Graph> Expand(const omnigraph::GraphComponent<debruijn_graph::Graph> &component) const {
        INFO("Expanding component to include incident edges of all 'inner' vertices");
        std::set<debruijn_graph::EdgeId> expanded_edges(component.edges());
        for (debruijn_graph::VertexId v : component.vertices()) {
            if (IsInnerVertex(v, component.edges())) {
                utils::insert_all(expanded_edges, g_.IncidentEdges(v));
            }
        }
        return omnigraph::GraphComponent<debruijn_graph::Graph>::FromEdges(g_,
                                                           expanded_edges.begin(),
                                                           expanded_edges.end());
    }

private:
    DECL_LOGGER("ComponentExpander");
};

static bool IsDeadEnd(const debruijn_graph::Graph &g, debruijn_graph::VertexId v) {
    return g.IncomingEdgeCount(v) * g.OutgoingEdgeCount(v) == 0;
}

static std::vector<debruijn_graph::EdgeId> CollectDeadEnds(const omnigraph::GraphComponent<debruijn_graph::Graph> &gc, bool canonical_only = true) {
    const auto &g = gc.g();
    std::vector<debruijn_graph::EdgeId> answer;
    for (debruijn_graph::EdgeId e : gc.edges()) {
        VERIFY(gc.edges().count(g.conjugate(e)));
        if (canonical_only && g.conjugate(e) < e)
            continue;
        if (IsDeadEnd(g, g.EdgeStart(e)) || IsDeadEnd(g, g.EdgeEnd(e))) {
            answer.push_back(e);
        }
    }
    return answer;
}

static void WriteComponentWithDeadends(const omnigraph::GraphComponent<debruijn_graph::Graph> &component,
                                       const std::string &prefix,
                                       const io::EdgeNamingF<debruijn_graph::Graph> &naming_f) {
    INFO("Writing GFA to " << prefix << ".gfa")
    std::ofstream os(prefix + ".gfa");
    const auto &g = component.g();
    gfa::GFAWriter writer(g, os, naming_f);
    writer.WriteSegmentsAndLinks(component);

    INFO("Writing dead-ends to " << prefix << ".deadends")
    std::ofstream deadends_os(prefix + ".deadends");
    for (debruijn_graph::EdgeId e : CollectDeadEnds(component)) {
        deadends_os << naming_f(g, e) << "\n";
    }
}

}
