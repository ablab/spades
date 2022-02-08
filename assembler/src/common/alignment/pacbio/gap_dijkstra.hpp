//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/mapping_path.hpp"

#include "sequence/sequence_tools.hpp"
#include "utils/perf/perfcounter.hpp"

namespace sensitive_aligner {

using debruijn_graph::EdgeId;

union DijkstraReturnCode {
    struct {
        bool queue_limit: 1;
        bool iter_limit: 1;
        bool not_connected: 1;
        bool no_path: 1;
        bool long_gap: 1;
        bool wide_gap: 1;
    };
    unsigned status = 0;
};

struct DijkstraParams {
    size_t queue_limit = 1000000;
    size_t iteration_limit = 1000000;
    size_t updates_limit = 1000000;
    bool find_shortest_path = true;
    bool restore_mapping = false;
    float penalty_ratio = 200;

    int max_ed_proportion = 3;
    int ed_lower_bound = 200;
    int ed_upper_bound = 1000;
};

struct GapClosingConfig: public DijkstraParams {
    bool run_dijkstra = false;
    size_t max_gs_states = 120000000;
};

struct EndsClosingConfig: public DijkstraParams {
    int max_restorable_length = 0;
};


struct GraphState {
    EdgeId e;
    int start_pos;
    int end_pos;

    GraphState(EdgeId e_, int start_pos_, int end_pos_)
        : e(e_), start_pos(start_pos_), end_pos(end_pos_)
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

    std::string str() const {
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
        : gs(EdgeId(), 0, 0), i(0), is_empty(true)
    {}

    QueueState(EdgeId e_, int start_pos_, int end_pos_, int i_)
        : gs(e_, start_pos_, end_pos_), i(i_), is_empty(false)
    {}

    QueueState(GraphState gs_, int i_)
        : gs(gs_), i(i_), is_empty(false)
    {}

    bool operator == (const QueueState &state) const {
        return (this->gs == state.gs && this->i == state.i);
    }

    bool operator != (const QueueState &state) const {
        return !(*this == state);
    }

    bool operator < (const QueueState &state) const {
        return (this->gs < state.gs || (this->gs == state.gs && this->i < state.i) );
    }

    bool empty() const {
        return is_empty;
    }

    std::string str() const {
        return gs.str() + " seq_ind=" + std::to_string(i) + " empty=" + std::to_string(is_empty);
    }
};

} // namespace sensitive_aligner

namespace std {

template <>
struct hash<sensitive_aligner::QueueState> {
    inline size_t hash_size_t_pair(size_t s0, size_t s1) const {
         s1 ^= s1 << 23;  // a
         return (s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26)) + s0;
    }

    size_t operator()(const sensitive_aligner::QueueState& k) const {
        return std::hash<size_t>()(hash_size_t_pair(
                                    hash_size_t_pair(k.gs.e.int_id(), k.gs.start_pos), 
                                    hash_size_t_pair(k.gs.end_pos, k.i)
                                ));
    }
};

}  // namespace std

namespace sensitive_aligner {

class DijkstraGraphSequenceBase {
  public:
    DijkstraGraphSequenceBase(const debruijn_graph::Graph &g,
                              const DijkstraParams &gap_cfg,
                              const std::string &ss,
                              EdgeId start_e, int start_p, int path_max_length)
        : g_(g)
        , gap_cfg_(gap_cfg)
        , ss_(ss)
        , start_e_(start_e)
        , start_p_(start_p)
        , path_max_length_(path_max_length)
        , min_score_(std::numeric_limits<int>::max())
        , queue_limit_(gap_cfg_.queue_limit)
        , iter_limit_(gap_cfg_.iteration_limit)
        , updates_(0) {
        best_ed_.resize(ss_.size(), path_max_length_);
        AddNewEdge(GraphState(start_e_, start_p_, (int) g_.length(start_e_)), QueueState(), 0);
    }

    void CloseGap();

    std::vector<EdgeId> path() const {
        return mapping_path_.simple_path();
    }

    omnigraph::MappingPath<EdgeId> mapping_path() const {
        return mapping_path_;
    }

    std::string path_str() const {
        std::string result = "";
        auto gap_path = mapping_path_.simple_path();
        for (EdgeId e : gap_path) {
            result += std::to_string(e.int_id()) + ",";
        }
        return result;
    }

    int edit_distance() const {
        return min_score_;
    }

    DijkstraReturnCode return_code() const {
        return return_code_;
    }

    int path_end_position() const {
        return end_qstate_.gs.end_pos;
    }

    int seq_end_position() const {
        return end_qstate_.i;
    }

    ~DijkstraGraphSequenceBase() {}

  protected:
    bool IsBetter(int seq_ind, int ed);

    void Update(const QueueState &state, const QueueState &prev_state, int score);

    void AddNewEdge(const GraphState &gs, const QueueState &prev_state, int ed);

    bool QueueLimitsExceeded(size_t iter);

    bool RunDijkstra();

    virtual bool AddState(const QueueState &cur_state, EdgeId e, int ed) = 0;

    virtual bool IsEndPosition(const QueueState &cur_state) = 0;

    omnigraph::MappingPath<EdgeId> mapping_path_;


    const debruijn_graph::Graph &g_;
    const DijkstraParams gap_cfg_;
    const std::string ss_;
    const EdgeId start_e_;
    const int start_p_;

    QueueState end_qstate_;

    int path_max_length_;
    int min_score_;
    DijkstraReturnCode return_code_;

  private:
    static const int SHORT_SEQ_LENGTH = 100;
    static const int ED_DEVIATION = 20;

    std::set<std::pair<int, QueueState>> q_;
    std::unordered_map<QueueState, int> visited_;
    std::unordered_map<QueueState, QueueState> prev_states_;
    std::vector<int> best_ed_;

    const size_t queue_limit_;
    const size_t iter_limit_;
    size_t updates_;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DijkstraGapFiller: public DijkstraGraphSequenceBase {
  public:
    DijkstraGapFiller(const debruijn_graph::Graph &g,
                      const GapClosingConfig &gap_cfg,
                      const std::string &ss,
                      EdgeId start_e, EdgeId end_e,
                      int start_p, int end_p, int path_max_length,
                      const std::unordered_map<debruijn_graph::VertexId, size_t> &reachable_vertex)
        : DijkstraGraphSequenceBase(g, gap_cfg, ss, start_e, start_p, path_max_length)
        , end_e_(end_e) , end_p_(end_p)
        , reachable_vertex_(reachable_vertex) {
        GraphState end_gstate(end_e_, 0, end_p_);
        end_qstate_ = QueueState(end_gstate, (int) ss_.size());
        if (start_e_ == end_e_ && end_p_ - start_p_ > 0) {
            std::string edge_str = g_.EdgeNucls(start_e_).Subseq(start_p_, end_p_).str();
            int score = StringDistance(ss_, edge_str, path_max_length_);
            if (score != std::numeric_limits<int>::max()) {
                path_max_length_ = std::min(path_max_length_, score);
                QueueState state(GraphState(start_e_, start_p_, end_p_), (int) ss_.size());
                Update(state, QueueState(), score);
                if (score == path_max_length_) {
                    end_qstate_ = state;
                }
            }
        }
    }

  private:

    bool AddState(const QueueState &cur_state, EdgeId e, int ed) override;

    bool IsEndPosition(const QueueState &cur_state) override;

    EdgeId end_e_;
    const int end_p_;
    const std::unordered_map<debruijn_graph::VertexId, size_t> &reachable_vertex_;
};



class DijkstraEndsReconstructor: public DijkstraGraphSequenceBase {
  public:
    DijkstraEndsReconstructor(const debruijn_graph::Graph &g,
                              const EndsClosingConfig &gap_cfg,
                              const std::string &ss,
                              EdgeId start_e, int start_p, int path_max_length)
        : DijkstraGraphSequenceBase(g, gap_cfg, ss, start_e, start_p, path_max_length) {
        end_qstate_ = QueueState();
        if (g_.length(start_e_) + g_.k() - start_p_ + path_max_length_ > ss_.size()) {
            std::string edge_str = g_.EdgeNucls(start_e_).Subseq(start_p_).str();
            int position = -1;
            int score = SHWDistance(ss_, edge_str, path_max_length, position);
            if (score != std::numeric_limits<int>::max()) {
                path_max_length_ = score;
                QueueState state(GraphState(start_e_, start_p_, start_p_ + position + 1), (int) ss_.size() );
                Update(state, QueueState(), score);
                end_qstate_ = state;
            }
        }

    }

  private:
    bool AddState(const QueueState &cur_state, EdgeId e, int ed) override;

    bool IsEndPosition(const QueueState &cur_state) override;
};

} // namespace sensitive_aligner
