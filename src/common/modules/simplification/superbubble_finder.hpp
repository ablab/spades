//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "sequence/range.hpp"
#include "assembly_graph/components/graph_component.hpp"

#include "utils/stl_utils.hpp"
#include "utils/logger/logger.hpp"
#include "math/xmath.h"
#include <queue>
#include <unordered_set>

namespace omnigraph {

template<class Graph>
class SuperbubbleFinder {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    VertexId start_vertex_;
    size_t max_length_;
    size_t max_count_;

    size_t cnt_;
    //TODO think of alternative definitions of weight (currently: total k-mer multiplicity)
    //vertex to heaviest path weight / path length range
    std::unordered_map<VertexId, std::pair<size_t, Range>> superbubble_vertices_;
    std::unordered_map<VertexId, EdgeId> heaviest_backtrace_;
    VertexId end_vertex_;

    bool CheckCanBeProcessed(VertexId v) const {
        DEBUG("Check if vertex " << g_.str(v) << " is dominated close neighbour");
        for (EdgeId e : g_.IncomingEdges(v)) {
            if (superbubble_vertices_.count(g_.EdgeStart(e)) == 0) {
                DEBUG("Blocked by external vertex " << g_.int_id(g_.EdgeStart(e)) << " that starts edge " << g_.int_id(e));
                DEBUG("Check fail");
                return false;
            }
        }
        DEBUG("Check ok");
        return true;
    }

    void UpdateCanBeProcessed(VertexId v,
                              std::unordered_set<VertexId>& can_be_processed,
                              std::unordered_set<VertexId>& border) const {
        DEBUG("Updating can be processed");
        for (EdgeId e : g_.OutgoingEdges(v)) {
            DEBUG("Considering edge " << g_.str(e));
            VertexId neighbour_v = g_.EdgeEnd(e);
            if (neighbour_v == start_vertex_) {
                VERIFY(v == start_vertex_);
                continue;
            }
            VERIFY(superbubble_vertices_.count(neighbour_v) == 0);
            DEBUG("Adding vertex " << g_.str(neighbour_v) << " to border");
            border.insert(neighbour_v);
            if (CheckCanBeProcessed(neighbour_v)) {
                DEBUG("Adding vertex " << g_.str(neighbour_v) << " to 'can be processed' set");
                can_be_processed.insert(neighbour_v);
            }
        }
    }

    bool CheckNoEdgeToStart(VertexId v) {
        for (EdgeId e : g_.OutgoingEdges(v)) {
            if (g_.EdgeEnd(e) == start_vertex_) {
                return false;
            }
        }
        return true;
    }

public:
    SuperbubbleFinder(const Graph& g, VertexId v, size_t max_length = -1ul,
                      size_t max_count = -1ul)
            : g_(g),
              start_vertex_(v),
              max_length_(max_length),
              max_count_(max_count),
              cnt_(0) {

    }

    //todo handle case when first/last vertex have other outgoing/incoming edges
    //true if no thresholds exceeded
    bool FindSuperbubble() {
        if (g_.OutgoingEdgeCount(start_vertex_) < 2) {
            return false;
        }
        DEBUG("Adding starting vertex " << g_.str(start_vertex_) << " to dominated set");
        superbubble_vertices_[start_vertex_] = std::make_pair(0, Range(0, 0));
        heaviest_backtrace_[start_vertex_] = EdgeId();
        cnt_++;
        std::unordered_set<VertexId> can_be_processed;
        std::unordered_set<VertexId> border;
        UpdateCanBeProcessed(start_vertex_, can_be_processed, border);
        while (!can_be_processed.empty()) {
            //finish after checks and adding the vertex
            bool final = (border.size() == 1 && can_be_processed.size() == 1);
            DEBUG("Final: " << final);
            if (++cnt_ > max_count_) {
                return false;
            }
            VertexId v = *can_be_processed.begin();
            can_be_processed.erase(can_be_processed.begin());

            DEBUG("Counting distance range for vertex " << g_.str(v));
            size_t min_d = std::numeric_limits<size_t>::max();
            size_t max_d = 0;
            size_t max_w = 0;
            EdgeId entry;

            VERIFY(g_.IncomingEdgeCount(v) > 0);
            VERIFY(CheckCanBeProcessed(v));
            for (EdgeId e : g_.IncomingEdges(v)) {
                //in case of dominated_only == false
                if (superbubble_vertices_.count(g_.EdgeStart(e)) == 0)
                    continue;
                size_t weight;
                Range range;
                std::tie(weight, range) = utils::get(superbubble_vertices_, g_.EdgeStart(e));
                range.shift((int) g_.length(e));
                DEBUG("Edge " << g_.str(e) << " provide distance range " << range);
                if (range.start_pos < min_d)
                    min_d = range.start_pos;
                if (range.end_pos > max_d)
                    max_d = range.end_pos;

                weight += size_t(math::round(double(g_.length(e)) * g_.coverage(e)));
                if (weight > max_w) {
                    max_w = weight;
                    entry = e;
                }
            }
            VERIFY((max_d > 0) && (min_d < std::numeric_limits<size_t>::max()) && (min_d <= max_d));
            DEBUG("Range " << Range(min_d, max_d));
            Range r(min_d, max_d);
//            return std::make_pair(std::make_pair(max_w, entry), Range(min_d, max_d));
            if (r.start_pos > max_length_) {
                return false;
            }
            //Inner vertices cannot have edge to start vertex
            //TODO Also all added edges have to have an outgoing edge
            if (!final && (!CheckNoEdgeToStart(v) || g_.OutgoingEdgeCount(v) == 0)) {
                return false;
            }

            DEBUG("Adding vertex " << g_.str(v) << " to dominated set");
            superbubble_vertices_[v] = std::make_pair(max_w, r);
            heaviest_backtrace_[v] = entry;
            border.erase(v);
            if (final) {
                end_vertex_ = v;
                return true;
            } else {
                UpdateCanBeProcessed(v, can_be_processed, border);
            }
        }
        DEBUG("Finished search for starting vertex " << g_.str(start_vertex_));
        return false;
    }

    const std::unordered_map<VertexId, Range>& bubble_vertices() const {
        return superbubble_vertices_;
    }

    GraphComponent<Graph> AsGraphComponent() const {
        return GraphComponent<Graph>::FromVertices(g_, utils::key_set(superbubble_vertices_));
    }

    Range PathLengthRange() const {
        return end_vertex_ == VertexId() ? Range() :
               utils::get(superbubble_vertices_, end_vertex_);
    }

    VertexId end_vertex() const {
        return end_vertex_;
    }

    const std::vector<EdgeId> HeaviestPath() const {
        VERIFY(end_vertex_ != VertexId());
        std::vector<EdgeId> rev_edges;
        VertexId v = end_vertex_;
        while (v != VertexId()) {
            EdgeId e = heaviest_backtrace_[v];
            rev_edges.push_back(e);
            v = g_.EdgeStart(e);
        }
        return std::vector<EdgeId>(rev_edges.rbegin(), rev_edges.rend());
    }

private:
    DECL_LOGGER("SuperbubbleFinder");
};
}
