//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "pe_utils.hpp"
#include "assembly_graph/core/graph.hpp"

using namespace debruijn_graph;

namespace path_extend {

bool GetLoopAndExit(const Graph& g, EdgeId forward_cycle_edge, EdgeId& back_cycle_edge, EdgeId& loop_outgoing, EdgeId& loop_incoming) {
    VertexId loop_end = g.EdgeEnd(forward_cycle_edge);
    VertexId loop_start = g.EdgeStart(forward_cycle_edge);
    if (g.OutgoingEdgeCount(loop_end) != 2 || g.IncomingEdgeCount(loop_end) != 1 ||
        g.OutgoingEdgeCount(loop_start) != 1 || g.IncomingEdgeCount(loop_start) != 2) {
        return false;
    }

    auto edges = g.OutgoingEdges(loop_end);
    EdgeId edge1 = *edges.begin();
    EdgeId edge2 = *(std::next(edges.begin()));
    if (g.EdgeEnd(edge1) == g.EdgeEnd(edge2)) {
//Patologic situation, two glued loops
        return false;
    }
    if (g.EdgeEnd(edge1) == g.EdgeStart(forward_cycle_edge)) {
        back_cycle_edge = edge1;
        loop_outgoing = edge2;
    } else if (g.EdgeEnd(edge2) == g.EdgeStart(forward_cycle_edge)) {
        back_cycle_edge = edge2;
        loop_outgoing = edge1;
    } else {
        return false;
    }

    for (EdgeId edge : g.IncomingEdges(loop_start)) {
        if (edge != back_cycle_edge) {
            loop_incoming = edge;
        }
    }

    return true;
}

void RemoveLoop(BidirectionalPath &path, const GraphCoverageMap& cov_map,
                size_t skip_identical_edges, bool fullRemoval) {
    size_t toRemove = LoopDetector(path, cov_map).EdgesToRemove(skip_identical_edges, fullRemoval);
    for (size_t i = 0; i < toRemove; ++i) {
        path.PopBack();
    }
}


LoopDetector::LoopDetector(const BidirectionalPath &p, const GraphCoverageMap& cov_map)
        : path_(p),
          cov_map_(cov_map) {
}

size_t LoopDetector::LoopEdges(size_t skip_identical_edges, size_t min_cycle_appearences) const {
    if (path_.Size() == 0) {
        return 0;
    }
    EdgeId e = path_.Back();
    size_t count = cov_map_.Count(e, path_);
    if (count <= 1 || count < min_cycle_appearences * (skip_identical_edges + 1)) {
        return 0;
    }
    auto edge_positions = path_.FindAll(e);
    VERIFY(edge_positions.size() == count);
    VERIFY(edge_positions.size() >= skip_identical_edges);
    size_t loopSize = edge_positions.back() - edge_positions[edge_positions.size() - 1 - (skip_identical_edges + 1)];
    return loopSize;
}

bool LoopDetector::PathIsLoop(size_t edges) const {
    if (edges == 0 || path_.Size() <= 1)
        return false;

    for (size_t i = 0; i < edges; ++i) {
        EdgeId e = path_.At(i);
        for (int j = (int) path_.Size() - ((int) edges - (int) i); j >= 0; j -= (int) edges) {
            if (path_.operator [](j) != e) {
                return false;
            }
        }
    }
    return true;
}

size_t LoopDetector::LastLoopCount(size_t skip_identical_edges, size_t min_cycle_appearences) const {
    size_t edges = LoopEdges(skip_identical_edges, min_cycle_appearences);
    return LastLoopCount(edges);
}

size_t LoopDetector::LastLoopCount(size_t edges) const {
    if (edges == 0)
        return 0;

    BidirectionalPath loop = path_.SubPath(path_.Size() - edges);
    size_t count = 0;
    int i = (int) path_.Size() - (int) edges;
    int delta = -(int) edges;

    while (i >= 0) {
        if (!path_.CompareFrom(i, loop)) {
            break;
        }
        ++count;
        i += delta;
    }

    return count;
}

bool LoopDetector::IsCycled(size_t loopLimit, size_t& skip_identical_edges) const {
    if (path_.Size() == 0 || cov_map_.Count(path_.Back(), path_) < loopLimit)
        return false;

    skip_identical_edges = 0;
    size_t loop_count = LastLoopCount(skip_identical_edges, loopLimit);
    while (loop_count > 0) {
        if (loop_count >= loopLimit)
            return true;

        loop_count = LastLoopCount(++skip_identical_edges, loopLimit);
    }
    return false;
}

size_t LoopDetector::EdgesToRemove(size_t skip_identical_edges, bool fullRemoval) const {
    size_t edges = LoopEdges(skip_identical_edges, 1);
    size_t count = LastLoopCount(edges);
    bool onlyCycle = PathIsLoop(edges);
    int result;

    if (onlyCycle || path_.Size() <= count * edges) {
        result = (int) path_.Size() - (int) edges;
    } else if (fullRemoval) {
        result = (int) count * (int) edges;
    } else {
        result = (int) (count - 1) * (int) edges;
    }

    return result < 0 ? 0 : result;
}

bool LoopDetector::EdgeInShortLoop(EdgeId e) const {
    EdgeId back_cycle_edge;
    EdgeId loop_exit;
    EdgeId loop_in;
    return GetLoopAndExit(path_.graph(), e, back_cycle_edge, loop_exit, loop_in);
}

bool LoopDetector::PrevEdgeInShortLoop() const {
    if (path_.Size() <= 2)
        return false;

    const Graph& g = path_.graph();
    EdgeId e2 = path_.At(path_.Size() - 1);
    EdgeId e1 = path_.At(path_.Size() - 2);
    VertexId v2 = g.EdgeEnd(e1);
    if (g.OutgoingEdgeCount(v2) == 2 &&
        g.EdgeEnd(e2) == g.EdgeStart(e1) &&
        g.EdgeEnd(e1) == g.EdgeStart(e2)) {
        return EdgeInShortLoop(e1);
    }
    
    return false;
}

}
