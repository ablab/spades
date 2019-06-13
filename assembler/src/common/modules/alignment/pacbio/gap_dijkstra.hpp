//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <deque>

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/mapping_path.hpp"

#include "sequence/sequence_tools.hpp"
#include "utils/perf/perfcounter.hpp"


namespace sensitive_aligner {
using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;

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
    size_t queue_limit = 0;
    size_t iteration_limit = 0;
    size_t updates_limit = 0;
    bool find_shortest_path = false;
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
    bool restore_ends = false;
    int max_restorable_length = 0;
};

struct ProteinAlignmentConfig: public DijkstraParams {
    double min_alignment_len = 0.2;
    std::string penalty_matrix = "blosum62";
    bool stop_codon = true;
    int max_restorable_length = 1000;
    int indel_score = 5;
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

    std::string str() const {
        return gs.str() + " seq_ind=" + std::to_string(i) + " empty=" + std::to_string(is_empty);
    }
};

struct ProteinQueueState: public QueueState {
    std::string offset;
    
    ProteinQueueState()
        : offset("")
    {}  

    ProteinQueueState(debruijn_graph::EdgeId e_, int start_pos_, int end_pos_, int i_, std::string offset_)
        : QueueState(e_, start_pos_, end_pos_, i_), offset(offset_)
    {}

    ProteinQueueState(GraphState gs_, int i_, std::string offset_)
        : QueueState(gs_, i_), offset(offset_)
    {} 

    bool operator == (const ProteinQueueState &state) const {
        return (this->gs == state.gs && this->i == state.i && this->offset == state.offset);
    }

    bool operator != (const ProteinQueueState &state) const {
        return (this->gs != state.gs || this->i != state.i || this->offset != state.offset);
    }

    bool operator < (const ProteinQueueState &state) const {
        return (this->gs < state.gs || (this->gs == state.gs && this->i < state.i)
               || (this->gs == state.gs && this->i < state.i && this->offset < state.offset) );
    }

    bool empty() const {
        return is_empty;
    }

    std::string str() const {
        return gs.str() + " seq_ind=" + std::to_string(i) + " empty=" + std::to_string(is_empty) + " offset=" + offset;
    }

};

}

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

template <>
struct hash<sensitive_aligner::ProteinQueueState> {
    inline size_t hash_size_t_pair(size_t s0, size_t s1) const {
         s1 ^= s1 << 23;  // a
         return (s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26)) + s0;
    }

    size_t operator()(const sensitive_aligner::ProteinQueueState& k) const {
        return std::hash<size_t>()(hash_size_t_pair(
                                    hash_size_t_pair(
                                        hash_size_t_pair(k.gs.e.int_id(), k.gs.start_pos), 
                                        hash_size_t_pair(k.gs.end_pos, k.i)),
                                    k.offset.size()
                                 ));
    }
};

}  // namespace std

namespace sensitive_aligner {

template<class T, typename U>
class DijkstraGraphSequenceBase {
public:
    DijkstraGraphSequenceBase(const debruijn_graph::Graph &g,
                              const DijkstraParams &gap_cfg,
                              std::string ss,
                              debruijn_graph::EdgeId start_e, int start_p, int path_max_length)
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
        best_ed_.resize(ss_.size() + 1, std::pair<T, int> (T(), path_max_length_));
        state_time =0;
        newstates_time = 0;
        m_time = 0;
        d_time = 0;
        i_time = 0;
    }   

    void CloseGap();

    std::vector<debruijn_graph::EdgeId> path() const {
        return mapping_path_.simple_path();
    }

    omnigraph::MappingPath<debruijn_graph::EdgeId> mapping_path() const {
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

protected:

    virtual bool DetermineBestPrefix(bool found_path) = 0;

    virtual void AddNewStatesFrom(const T &state, int score) = 0;

    bool QueueLimitsExceeded(size_t iter);

    bool RunDijkstra();

    virtual bool IsEndPosition(const T &cur_state) = 0;

    virtual bool IsGoodEdge(const EdgeId e) = 0;

    virtual void PopFront() = 0;

    U q_;
    std::unordered_map<T, int> visited_;
    std::unordered_map<T, T> prev_states_;

    omnigraph::MappingPath<debruijn_graph::EdgeId> mapping_path_;

    std::vector<std::pair<T, int>> best_ed_;

    const debruijn_graph::Graph &g_;
    const DijkstraParams gap_cfg_;
    
    const std::string ss_;
    const debruijn_graph::EdgeId start_e_;
    const int start_p_;
    int path_max_length_;
    int min_score_;
    T end_qstate_;
    DijkstraReturnCode return_code_;

    const size_t queue_limit_;
    const size_t iter_limit_;
    size_t updates_;
    double state_time;
    double newstates_time;
    double m_time;
    double i_time;
    double d_time;
};
////
class DijkstraReadGraph: public DijkstraGraphSequenceBase<QueueState, std::deque<std::pair<int, QueueState> > > {

public:
    DijkstraReadGraph(const debruijn_graph::Graph &g,
                              const DijkstraParams &gap_cfg,
                              std::string ss,
                              debruijn_graph::EdgeId start_e, int start_p, int path_max_length)
    : DijkstraGraphSequenceBase(g, gap_cfg, ss, start_e, start_p, path_max_length)
    {
        AddState(QueueState(GraphState(start_e_, start_p_, start_p_), 0), 0, QueueState(), true);
    }

    ~DijkstraReadGraph() {}

protected:

    void AddState(const QueueState &state, int score, const QueueState &prev_state, bool front);

    virtual void AddNewStatesFrom(const QueueState &state, int score);

    virtual void PopFront() {
        q_.pop_front();
    };
};



class DijkstraProteinGraph: public DijkstraGraphSequenceBase<ProteinQueueState, std::set<std::pair<int, ProteinQueueState> > > {

public:
    DijkstraProteinGraph(const debruijn_graph::Graph &g,
                              const ProteinAlignmentConfig &gap_cfg,
                              std::string ss,
                              debruijn_graph::EdgeId start_e, int start_p, int path_max_length, bool is_reverse)
    :DijkstraGraphSequenceBase(g, gap_cfg, ss, start_e, start_p, path_max_length),
    is_reverse_(is_reverse),
    stop_stop_codon_(gap_cfg.stop_codon),
    indel_score_(gap_cfg.indel_score),
    matrix_(NULL),
    min_value_(0)
    {
        const parasail_matrix_t *matrix = NULL;
        matrix = parasail_matrix_lookup(gap_cfg.penalty_matrix.c_str());
        matrix_ = parasail_matrix_copy(matrix);
        const int proteins_num = 24;
        for (int i = 0; i < proteins_num; ++ i){
            for (int j = 0; j < proteins_num; ++ j){
                if (-matrix_->matrix[i*proteins_num + j] < min_value_) {
                    min_value_ = -matrix_->matrix[i*proteins_num + j];
                }
            }
        }
        min_value_ = abs(min_value_);
        DEBUG("min_value=" << min_value_)
        path_max_length_ = min_value_* ((int) ss_.size());
        DEBUG("max_path_len=" << path_max_length_)
        for (size_t i = 0; i < best_ed_.size(); ++ i) {
            best_ed_[i].second = path_max_length_;
        }
        AddState(ProteinQueueState(GraphState(start_e_, start_p_, start_p_), 0, ""), 0, ProteinQueueState());
    }

    ~DijkstraProteinGraph() {
        parasail_matrix_free(matrix_);
    }

protected:

    void AddState(const ProteinQueueState &state, int score, const ProteinQueueState &prev_state);

    void AddNewStatesFrom(const ProteinQueueState &state, int score);

    bool StopCodon(const std::string &a) const {
        if (!stop_stop_codon_) return false;
        if (!is_reverse_) return stop_codons.count(a) > 0;
        return cstop_codons.count(a) > 0;
    }

    int deletion_score() { return min_value_ + indel_score_;}

    int insertion_score() { return indel_score_;}

    int mm_score(const std::string &a, const std::string &b) const {
        if (!is_reverse_) {
            return min_value_ + PenaltyMatrixTriplets(a, b, matrix_);
        }
        std::string new_a = (!Sequence(a)).str();
        std::string new_b = (!Sequence(b)).str();
        return min_value_ + PenaltyMatrixTriplets(new_a, new_b, matrix_);
    }

    virtual void PopFront() {
        q_.erase(q_.begin());
    };

    const bool is_reverse_;
    const bool stop_stop_codon_;
    const int indel_score_;
    parasail_matrix_t *matrix_;
    int min_value_;
    static const std::set<std::string> stop_codons;
    static const std::set<std::string> cstop_codons;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DijkstraGapFiller: public DijkstraReadGraph {

public:
    DijkstraGapFiller(const debruijn_graph::Graph &g, const GapClosingConfig &gap_cfg, std::string ss,
                      debruijn_graph::EdgeId start_e, debruijn_graph::EdgeId end_e,
                      int start_p, int end_p, int path_max_length,
                      const std::map<debruijn_graph::VertexId, size_t> &reachable_vertex)
        : DijkstraReadGraph(g, gap_cfg, ss, start_e, start_p, path_max_length)
        , end_e_(end_e) , end_p_(end_p)
        , reachable_vertex_(reachable_vertex) {
        GraphState end_gstate(end_e_, 0, end_p_);
        end_qstate_ = QueueState(end_gstate, (int) ss_.size());
    }

private:

    bool DetermineBestPrefix(bool found_path) override;

    virtual bool IsEndPosition(const QueueState &cur_state);

    virtual bool IsGoodEdge(const EdgeId e) {
        return reachable_vertex_.count(g_.EdgeEnd(e)) > 0 || end_e_ == e;
    }

    debruijn_graph::EdgeId end_e_;
    const int end_p_;
    const std::map<debruijn_graph::VertexId, size_t> &reachable_vertex_;
};



class DijkstraEndsReconstructor: public DijkstraReadGraph {
public:
    DijkstraEndsReconstructor(const debruijn_graph::Graph &g, const EndsClosingConfig &gap_cfg, std::string ss,
                              debruijn_graph::EdgeId start_e, int start_p, int path_max_length)
        : DijkstraReadGraph(g, gap_cfg, ss, start_e, start_p, path_max_length) {
        end_qstate_ = QueueState();
    }

 private:

    bool DetermineBestPrefix(bool found_path) override;

    virtual bool IsEndPosition(const QueueState &cur_state);

    virtual bool IsGoodEdge(const EdgeId e) {
        (void) e;
        return true;
    }
};


class DijkstraProteinEndsReconstructor: public DijkstraProteinGraph {
public:
    DijkstraProteinEndsReconstructor(const debruijn_graph::Graph &g, const ProteinAlignmentConfig &gap_cfg, std::string ss,
                              debruijn_graph::EdgeId start_e, int start_p, int path_max_length, bool is_reverse = false)
        : DijkstraProteinGraph(g, gap_cfg, ss, start_e, start_p, path_max_length, is_reverse) {
        end_qstate_ = ProteinQueueState();
    }

 private:

    bool DetermineBestPrefix(bool found_path) override;

    virtual bool IsEndPosition(const ProteinQueueState &cur_state);

    virtual bool IsGoodEdge(const EdgeId e) {
        (void) e;
        return true;
    }
};

} // namespace sensitive_aligner