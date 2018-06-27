#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/mapping_path.hpp"

#include "sequence/sequence_tools.hpp"
#include "utils/perf/perfcounter.hpp"

namespace graph_aligner {
using debruijn_graph::EdgeId;

enum DijkstraReturnCode {
    OK = 0,
    NOT_CONNECTED = 1,
    TOO_LONG_GAP = 2,
    TOO_MANY_VERTICES = 4,
    QUEUE_LIMIT = 8,
    ITERATION_LIMIT = 16,
    NO_PATH = 32
};

struct GapClosingConfig {
    bool run_dijkstra = false;
    bool restore_ends = false;
    // dijkstra run parameters
    size_t max_vertex_in_gap = 0;
    size_t queue_limit = 0;
    size_t iteration_limit = 0;
    bool find_shortest_path = false;
    bool restore_mapping = false;
    int penalty_interval = 200;

    int max_ed_proportion = 3;
    int ed_lower_bound = 200;
    int ed_upper_bound = 1000;
    int max_restorable_end_length = 3000;
};


struct GraphState {
    debruijn_graph::EdgeId e;
    int start_pos;
    int end_pos;

    GraphState(debruijn_graph::EdgeId e_, int start_pos_, int end_pos_)
        : e(e_), start_pos(start_pos_), end_pos(end_pos_)
    {}

    GraphState(const GraphState &state):
        e(state.e), start_pos(state.start_pos), end_pos(state.end_pos)
    {}

    bool operator == (const GraphState &state) const {
        return (this->e == state.e && this->start_pos == state.start_pos && this->end_pos == state.end_pos);
    }

    bool operator != (const GraphState &state) const {
        return (this->e != state.e || this->start_pos != state.start_pos || this->end_pos != state.end_pos);
    }

    bool operator < (const GraphState &state) const {
        return (this->e < state.e
                || (this->e == state.e && this->start_pos < state.start_pos)
                || (this->e == state.e && this->start_pos == state.start_pos && this->end_pos < state.end_pos));
    }

    string str() const {
        return "edge_id=" + std::to_string(e.int_id())
               + " " + std::to_string(start_pos)
               + "-" + std::to_string(end_pos);
    }
};


struct QueueState {
    GraphState gs;
    int i;
    bool is_empty;

    QueueState()
        : gs(debruijn_graph::EdgeId(), 0, 0), i(0), is_empty(true)
    {}

    QueueState(debruijn_graph::EdgeId &e_, int start_pos_, int end_pos_, int i_)
        : gs(e_, start_pos_, end_pos_), i(i_), is_empty(false)
    {}

    QueueState(GraphState gs_, int i_)
        : gs(gs_), i(i_), is_empty(false)
    {}

    QueueState(const QueueState &state):
        gs(state.gs), i(state.i), is_empty(state.is_empty)
    {}

    bool operator == (const QueueState &state) const {
        return (this->gs == state.gs && this->i == state.i);
    }

    bool operator != (const QueueState &state) const {
        return (this->gs != state.gs || this->i != state.i);
    }

    bool operator < (const QueueState &state) const {
        return (this->gs < state.gs || (this->gs == state.gs && this->i < state.i) );
    }

    bool empty() const {
        return is_empty;
    }

    string str() const {
        return gs.str() + " seq_ind=" + std::to_string(i) + " empty=" + std::to_string(is_empty);
    }
};

struct StateHasher {
    std::size_t operator()(const QueueState& k) const
    {
        using std::hash;

        return ((hash<int>()(k.gs.e.int_id())
                 ^ (hash<int>()(k.gs.start_pos) << 1)) >> 1)
               ^ (hash<int>()(k.gs.end_pos) << 1)
               ^ (hash<int>()(k.i) << 17);
    }
};


class DijkstraGraphSequenceBase {

protected:
    bool IsBetter(int seq_ind, int ed);

    void Update(const QueueState &state, const QueueState &prev_state, int score);

    void AddNewEdge(const GraphState &gs, const QueueState &prev_state, int ed);

    bool QueueLimitsExceeded(size_t iter);

    bool RunDijkstra();

    virtual bool AddState(const QueueState &cur_state, debruijn_graph::EdgeId e, int ed) = 0;

    virtual bool IsEndPosition(const QueueState &cur_state) = 0;

public:
    DijkstraGraphSequenceBase(const debruijn_graph::Graph &g,
                              const GapClosingConfig &gap_cfg,
                              string ss,
                              debruijn_graph::EdgeId start_e, int start_p, int path_max_length)
        : g_(g)
        , gap_cfg_(gap_cfg)
        , ss_(ss)
        , start_e_(start_e)
        , start_p_(start_p)
        , path_max_length_(path_max_length)
        , min_score_(std::numeric_limits<int>::max())
        , return_code_(0)
        , queue_limit_(gap_cfg_.queue_limit)
        , iter_limit_(gap_cfg_.iteration_limit)
        , updates_(0) {
        best_ed_.resize(ss_.size(), path_max_length_);
        AddNewEdge(GraphState(start_e_, start_p_, (int) g_.length(start_e_)), QueueState(), 0);
    }

    void CloseGap();

    vector<debruijn_graph::EdgeId> path() const {
        return mapping_path_.simple_path();
    }

    omnigraph::MappingPath<debruijn_graph::EdgeId> mapping_path() const {
        return mapping_path_;
    }

    string path_str() const {
        string result = "";
        auto gap_path = mapping_path_.simple_path();
        for (EdgeId e : gap_path) {
            result += std::to_string(e.int_id()) + ",";
        }
        return result;
    }

    int edit_distance() const {
        return min_score_;
    }

    int return_code() const {
        return return_code_;
    }

    int path_end_position() const {
        DEBUG("End position edge=" << end_qstate_.gs.e.int_id() << " end_pos=" << end_qstate_.gs.end_pos
              << " seq_pos=" << end_qstate_.i << " s_len=" << ss_.size())
        return end_qstate_.gs.end_pos;
    }

    int seq_end_position() const {
        DEBUG("End position edge=" << end_qstate_.gs.e.int_id() << " end_pos=" << end_qstate_.gs.end_pos
              << " seq_pos=" << end_qstate_.i << " s_len=" << ss_.size())
        return end_qstate_.i;
    }

    ~DijkstraGraphSequenceBase() {
        DEBUG("updates=" << updates_)
    }

protected:
    set<pair<int, QueueState> > q_;
    unordered_map<QueueState, int, StateHasher> visited_;
    unordered_map<QueueState, QueueState, StateHasher> prev_states_;

    omnigraph::MappingPath<debruijn_graph::EdgeId> mapping_path_;

    vector<int> best_ed_;

    const debruijn_graph::Graph &g_;
    const GapClosingConfig gap_cfg_;
    const string ss_;
    const debruijn_graph::EdgeId start_e_;
    const int start_p_;
    int path_max_length_;
    int min_score_;
    QueueState end_qstate_;
    int return_code_;

    const size_t queue_limit_;
    const size_t iter_limit_;
    int updates_;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DijkstraGapFiller: public DijkstraGraphSequenceBase {

private:

    virtual bool AddState(const QueueState &cur_state, debruijn_graph::EdgeId e, int ed);

    virtual bool IsEndPosition(const QueueState &cur_state);

public:
    DijkstraGapFiller(const debruijn_graph::Graph &g, const GapClosingConfig &gap_cfg, string ss,
                      debruijn_graph::EdgeId start_e, debruijn_graph::EdgeId end_e,
                      int start_p, int end_p, int path_max_length,
                      const map<debruijn_graph::VertexId, size_t> &reachable_vertex)
        : DijkstraGraphSequenceBase(g, gap_cfg, ss, start_e, start_p, path_max_length)
        , end_e_(end_e) , end_p_(end_p)
        , reachable_vertex_(reachable_vertex) {
        GraphState end_gstate(end_e_, 0, end_p_);
        end_qstate_ = QueueState(end_gstate, (int) ss_.size());
        if (start_e_ == end_e_ && end_p_ - start_p_ > 0) {
            string seq_str = ss_;
            string tmp = g_.EdgeNucls(start_e_).str();
            string edge_str = tmp.substr(start_p_, end_p_ - start_p_);
            int score = StringDistance(seq_str, edge_str, path_max_length_);
            if (score != std::numeric_limits<int>::max()) {
                path_max_length_ = min(path_max_length_, score);
                QueueState state(GraphState(start_e_, start_p_, end_p_), (int) ss_.size());
                Update(state, QueueState(), score);
                if (score == path_max_length_) {
                    end_qstate_ = state;
                }
            }
        }
    }

protected:
    debruijn_graph::EdgeId end_e_;
    const int end_p_;
    const map<debruijn_graph::VertexId, size_t> &reachable_vertex_;
};



class DijkstraEndsReconstructor: public DijkstraGraphSequenceBase {

private:
    virtual bool AddState(const QueueState &cur_state, debruijn_graph::EdgeId e, int ed);

    virtual bool IsEndPosition(const QueueState &cur_state);
public:
    DijkstraEndsReconstructor(const debruijn_graph::Graph &g, const GapClosingConfig &gap_cfg, string ss,
                              debruijn_graph::EdgeId start_e, int start_p, int path_max_length)
        : DijkstraGraphSequenceBase(g, gap_cfg, ss, start_e, start_p, path_max_length) {
        end_qstate_ = QueueState();
        if (g_.length(start_e_) + g_.k() - start_p_ + path_max_length_ > ss_.size()) {
            string seq_str = ss_;
            string tmp = g_.EdgeNucls(start_e_).str();
            string edge_str = tmp.substr(start_p_);
            int position = -1;
            int score = SHWDistance(seq_str, edge_str, path_max_length_, position);
            if (score != std::numeric_limits<int>::max()) {
                path_max_length_ = min(path_max_length_, score);
                QueueState state(GraphState(start_e_, start_p_, start_p_ + position + 1), (int) ss_.size() );
                Update(state, QueueState(), score);
                if (score == path_max_length_) {
                    end_qstate_ = state;
                }
            }
        }

    }
};

} // namespace graph_aligner