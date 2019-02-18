#include "omnigraph_wrapper.hpp"

#include "cursor.hpp"
#include "hmmpath.hpp"
#include "utils.hpp"

#include "debug_assert/debug_assert.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/components/graph_component.hpp"

using namespace debruijn_graph;

struct omnigraph_assert : debug_assert::default_handler,
                          debug_assert::set_level<1> {};

std::vector<DebruijnGraphCursor> DebruijnGraphCursor::prev(const void *context) const {
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


inline std::vector<DebruijnGraphCursor> DebruijnGraphCursor::next(const void *context) const {
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

std::vector<DebruijnComponentCursor> DebruijnComponentCursor::prev(const void *context) const {
    const ConjugateDeBruijnGraph &g = this->g(context);

    // Case 1: edge is a tip and we're inside the terminal vertex
    if (position() == 0) {
        // assert(pc_->g().ingoing_[edge_id_].size() == 0);
        return {};
    }

    // Case 2: move backwards possibly going inside the terminal vertex of a
    // tip
    if (position() != g.k() ||
        (g.IncomingEdgeCount(g.EdgeStart(edge())) == 0)) {
        return { DebruijnComponentCursor(edge(), position() - 1) };
    }

    // Case 3: go into incoming edges
    DEBUG_ASSERT(position() == g.k(), omnigraph_assert{});
    DEBUG_ASSERT(g.IncomingEdgeCount(g.EdgeStart(edge())), omnigraph_assert{});
    std::vector<DebruijnComponentCursor> result;
    for (EdgeId in : g.IncomingEdges(g.EdgeStart(edge()))) {
        if (!this->c(context).contains(in))
            continue; // We need to stay in component
        result.emplace_back(in, g.length(in) + g.k() - 1);
    }

    return result;
}


inline std::vector<DebruijnComponentCursor> DebruijnComponentCursor::next(const void *context) const {
    const ConjugateDeBruijnGraph &g = this->g(context);

    // Common case: we have not reached the end of the edge (in nucls)
    if (position() + 1 < g.length(edge()) + g.k())
        return { DebruijnComponentCursor(edge(), position() + 1) };

    // Otherwise we're inside the vertex and need to go out of it
    DEBUG_ASSERT(position() + 1 == g.length(edge()) + g.k(), omnigraph_assert{});
    std::vector<DebruijnComponentCursor> result;
    result.reserve(4);
    for (EdgeId out : g.OutgoingEdges(g.EdgeEnd(edge()))) {
        if (!this->c(context).contains(out))
            continue; // We need to stay in component
        result.emplace_back(out, g.k());  // Vertices are k-mers
    }

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

std::vector<DebruijnComponentCursor> DebruijnComponentCursor::all(const omnigraph::GraphComponent<ConjugateDeBruijnGraph> &c) {
    std::vector<DebruijnComponentCursor> result;
    for (auto it = c.e_begin(); it != c.e_end(); ++it) {
        EdgeId e = *it;
        for (size_t pos = 0; pos < c.g().length(e) + c.g().k(); ++pos) {
            DebruijnComponentCursor cursor{e, pos};
            cursor.normalize_prefix_to_suffix(&c);
            result.push_back(std::move(cursor));
        }
    }
    // TODO Filter duplicates
    return result;
}

PathSet<DebruijnGraphCursor>
find_best_path(const hmm::Fees &fees,
               const std::vector<DebruijnGraphCursor> &initial,
               const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<RestrictedGraphCursor<DebruijnGraphCursor>>
find_best_path(const hmm::Fees &fees,
               const std::vector<RestrictedGraphCursor<DebruijnGraphCursor>> &initial,
               const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<DebruijnComponentCursor>
find_best_path(const hmm::Fees &fees,
               const std::vector<DebruijnComponentCursor> &initial,
               const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<ReversalGraphCursor<DebruijnGraphCursor>>
find_best_path_rev(const hmm::Fees &fees,
                   const std::vector<ReversalGraphCursor<DebruijnGraphCursor>> &initial,
                   const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<ReversalGraphCursor<DebruijnComponentCursor>>
find_best_path_rev(const hmm::Fees &fees,
                   const std::vector<ReversalGraphCursor<DebruijnComponentCursor>> &initial,
                   const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<AAGraphCursor<DebruijnComponentCursor>>
find_best_path(const hmm::Fees &fees,
               const std::vector<AAGraphCursor<DebruijnComponentCursor>> &initial,
               const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<AAGraphCursor<DebruijnGraphCursor>>
find_best_path(const hmm::Fees &fees,
               const std::vector<AAGraphCursor<DebruijnGraphCursor>> &initial,
               const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<AAGraphCursor<RestrictedGraphCursor<DebruijnGraphCursor>>>
find_best_path(const hmm::Fees &fees,
               const std::vector<AAGraphCursor<RestrictedGraphCursor<DebruijnGraphCursor>>> &initial,
               const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>
find_best_path(const hmm::Fees &fees,
               const std::vector<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>> &initial,
               const void *context) {
    return impl::find_best_path(fees, initial, context);
}

PathSet<AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>>
find_best_path(const hmm::Fees &fees,
               const std::vector<AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>> &initial,
               const void *context) {
    return impl::find_best_path(fees, initial, context);
}

PathSet<StringCursor>
find_best_path(const hmm::Fees &fees,
               const std::vector<StringCursor> &initial,
               const void *context) {
    return impl::find_best_path(fees, initial, context);
}
