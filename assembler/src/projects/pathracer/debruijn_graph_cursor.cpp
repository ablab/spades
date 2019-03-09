#include "debruijn_graph_cursor.hpp"

#include "utils.hpp"

#include "debug_assert/debug_assert.hpp"

#include "assembly_graph/core/graph.hpp"

using namespace debruijn_graph;

struct omnigraph_assert : debug_assert::default_handler,
                          debug_assert::set_level<1> {};

std::vector<DebruijnGraphCursor> DebruijnGraphCursor::prev(DebruijnGraphCursor::Context context) const {
    const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);

    // Case 1: edge is a tip and we're inside the terminal vertex
    if (position() == 0) {
        // assert(pg_->ingoing_[edge_id_].size() == 0);
        return {};
    }

    // Case 2: move backwards possibly going inside the terminal vertex of a
    // tip
    if (position() != g.k() ||
        g.IncomingEdgeCount(g.EdgeStart(edge())) == 0) {
        return { DebruijnGraphCursor(edge(), position() - 1) };
    }

    // Case 3: go into incoming edges
    DEBUG_ASSERT(position() == g.k(), omnigraph_assert{});
    DEBUG_ASSERT(g.IncomingEdgeCount(g.EdgeStart(edge())), omnigraph_assert{});
    std::vector<DebruijnGraphCursor> result;
    for (EdgeId in : g.IncomingEdges(g.EdgeStart(edge())))
        result.emplace_back(in, g.length(in) + g.k() - 1);

    return result;
}


std::vector<DebruijnGraphCursor> DebruijnGraphCursor::next(DebruijnGraphCursor::Context context) const {
    const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);

    // Common case: we have not reached the end of the edge (in nucls)
    if (position() + 1 < g.length(edge()) + g.k())
        return { DebruijnGraphCursor(edge(), position() + 1) };

    // Otherwise we're inside the vertex and need to go out of it
    DEBUG_ASSERT(position() + 1 == g.length(edge()) + g.k(), omnigraph_assert{});
    std::vector<DebruijnGraphCursor> result;
    result.reserve(4);
    for (EdgeId out : g.OutgoingEdges(g.EdgeEnd(edge())))
        result.emplace_back(out, g.k());  // Vertices are k-mers

    return result;
}

std::vector<DebruijnGraphCursor> DebruijnGraphCursor::all(const ConjugateDeBruijnGraph &g) {
    std::vector<DebruijnGraphCursor> result;
    for (auto it = g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        for (size_t pos = 0; pos < g.length(e) + g.k(); ++pos) {
            DebruijnGraphCursor cursor{e, pos};
            cursor.normalize_prefix_to_suffix(&g);
            result.push_back(std::move(cursor));
        }
    }
    // TODO Filter duplicates
    return result;
}
