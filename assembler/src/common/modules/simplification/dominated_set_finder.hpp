#pragma once

namespace omnigraph {
template<class Graph>
class DominatedSetFinder {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    VertexId start_vertex_;
    size_t max_length_;
    size_t max_count_;

    size_t cnt_;
    std::map<VertexId, Range> dominated_;

    bool CheckCanBeProcessed(VertexId v) const {
        DEBUG("Check if vertex " << g_.str(v) << " is dominated close neighbour");
        for (EdgeId e : g_.IncomingEdges(v)) {
            if (dominated_.count(g_.EdgeStart(e)) == 0) {
                DEBUG("Blocked by external vertex " << g_.int_id(g_.EdgeStart(e)) << " that starts edge " << g_.int_id(e));
                DEBUG("Check fail");
                return false;
            }
        }
        DEBUG("Check ok");
        return true;
    }

    void UpdateCanBeProcessed(VertexId v,
                              std::queue<VertexId>& can_be_processed) const {
        DEBUG("Updating can be processed");
        for (EdgeId e : g_.OutgoingEdges(v)) {
            DEBUG("Considering edge " << g_.str(e));
            VertexId neighbour_v = g_.EdgeEnd(e);
            if (CheckCanBeProcessed(neighbour_v)) {
                can_be_processed.push(neighbour_v);
            }
        }
    }

    Range NeighbourDistanceRange(VertexId v, bool dominated_only = true) const {
        DEBUG("Counting distance range for vertex " << g_.str(v));
        size_t min = numeric_limits<size_t>::max();
        size_t max = 0;
        VERIFY(g_.IncomingEdgeCount(v) > 0);
        VERIFY(!dominated_only || CheckCanBeProcessed(v));
        for (EdgeId e : g_.IncomingEdges(v)) {
            //in case of dominated_only == false
            if (dominated_.count(g_.EdgeStart(e)) == 0)
                continue;
            Range range = dominated_.find(g_.EdgeStart(e))->second;
            range.shift((int) g_.length(e));
            DEBUG("Edge " << g_.str(e) << " provide distance range " << range);
            if (range.start_pos < min)
                min = range.start_pos;
            if (range.end_pos > max)
                max = range.end_pos;
        }
        VERIFY((max > 0) && (min < numeric_limits<size_t>::max()) && (min <= max));
        Range answer(min, max);
        DEBUG("Range " << answer);
        return answer;
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
    DominatedSetFinder(const Graph& g, VertexId v, size_t max_length = -1ul,
                       size_t max_count = -1ul)
            : g_(g),
              start_vertex_(v),
              max_length_(max_length),
              max_count_(max_count),
              cnt_(0) {

    }

    //true if no thresholds exceeded
    bool FillDominated() {
        DEBUG("Adding starting vertex " << g_.str(start_vertex_) << " to dominated set");
        dominated_.insert(make_pair(start_vertex_, Range(0, 0)));
        cnt_++;
        std::queue<VertexId> can_be_processed;
        UpdateCanBeProcessed(start_vertex_, can_be_processed);
        while (!can_be_processed.empty()) {
            if (++cnt_ > max_count_) {
                return false;
            }
            VertexId v = can_be_processed.front();
            can_be_processed.pop();
            Range r = NeighbourDistanceRange(v);
            if (r.start_pos > max_length_) {
                return false;
            }
            //Currently dominated vertices cannot have edge to start vertex
            if (CheckNoEdgeToStart(v)) {
                DEBUG("Adding vertex " << g_.str(v) << " to dominated set");
                dominated_.insert(make_pair(v, r));
                UpdateCanBeProcessed(v, can_be_processed);
            }
        }
        return true;
    }

    const map<VertexId, Range>& dominated() const {
        return dominated_;
    }

    GraphComponent<Graph> AsGraphComponent() const {
        return GraphComponent<Graph>::FromVertices(g_, utils::key_set(dominated_));
    }

    //little meaning if FillDominated returned false
    const map<VertexId, Range> CountBorder() const {
        map<VertexId, Range> border;
        for (VertexId v : utils::key_set(border)) {
            for (EdgeId e : g_.OutgoingEdges(v)) {
                VertexId e_end = g_.EdgeEnd(e);
                if (dominated_.count(e_end) == 0) {
                    border[e_end] = NeighbourDistanceRange(e_end, false);
                }
            }
        }
        return border;
    }

};
}
