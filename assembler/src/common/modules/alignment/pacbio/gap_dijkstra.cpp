//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "modules/alignment/pacbio/gap_dijkstra.hpp"

namespace sensitive_aligner {

using namespace std;

const int DijkstraGraphSequenceBase::SHORT_SEQ_LENGTH;
const int DijkstraGraphSequenceBase::ED_DEVIATION;

bool DijkstraGraphSequenceBase::IsBetter(int seq_ind, int ed) {
    if (seq_ind == (int) ss_.size() ) {
        if (ed <= path_max_length_) {
            return true;
        }
        return false;
    }
    VERIFY(seq_ind < (int) ss_.size())
    VERIFY(seq_ind >= 0)
    if (seq_ind < SHORT_SEQ_LENGTH ||
            max(best_ed_[seq_ind] + (int) ((double)seq_ind * gap_cfg_.penalty_ratio), ED_DEVIATION) >= ed) {
        best_ed_[seq_ind] = min(best_ed_[seq_ind], ed);
        return true;
    }
    return false;
}

void DijkstraGraphSequenceBase::Update(const QueueState &state, const QueueState &prev_state, int score) {
    if (visited_.count(state) > 0) {
        if (visited_[state] >= score) {
            ++ updates_;
            q_.erase(make_pair(visited_[state], state));
            if (IsBetter(state.i, score)) {
                q_.insert(make_pair(score, state));
                visited_[state] = score;
                prev_states_[state] = prev_state;
            }
        }
    } else {
        if (IsBetter(state.i, score)) {
            ++ updates_;
            visited_.insert(make_pair(state, score));
            prev_states_.insert(make_pair(state, prev_state));
            q_.insert(make_pair(score, state));
        }
    }
}

void DijkstraGraphSequenceBase::AddNewEdge(const GraphState &gs, const QueueState &prev_state, int ed) {
    string edge_str = g_.EdgeNucls(gs.e).Subseq(gs.start_pos, gs.end_pos).str();
    if (0 == edge_str.size()) {
        QueueState state(gs, prev_state.i);
        Update(state, prev_state,  ed);
        return;
    }
    if (path_max_length_ - ed >= 0) {
        if (path_max_length_ - ed >= (int) edge_str.size()) {
            QueueState state(gs, prev_state.i);
            Update(state, prev_state,  ed + (int) edge_str.size());
        }
    }
    if (ss_.size() - prev_state.i > 0) {
        // len - is a maximum length of substring to align on current edge
        int len = min( (int) g_.length(gs.e) - gs.start_pos + path_max_length_, // length of current edge + maximum insertion size
                       (int) ss_.size() - prev_state.i  ); // length of suffix left
        string seq_str = ss_.substr(prev_state.i, len);
        vector<int> positions;
        vector<int> scores;
        if (path_max_length_ - ed >= 0) {
            SHWDistanceExtended(seq_str, edge_str, path_max_length_ - ed, positions, scores);
            int prev_score = numeric_limits<int>::max();
            for (size_t i = 0; i < positions.size(); ++ i) {
                if (positions[i] >= 0 && scores[i] >= 0) {
                    int next_score = i + 1 >= positions.size() || positions[i + 1] < 0 || scores[i + 1] < 0 ?
                                     numeric_limits<int>::max() : scores[i + 1];
                    if (scores[i] <= prev_score && scores[i] <= next_score) {
                        QueueState state(gs, prev_state.i + positions[i] + 1);
                        Update(state, prev_state, ed + scores[i]);
                    }
                    prev_score = scores[i];
                } else {
                    prev_score = numeric_limits<int>::max();
                }
            }
        }
    }
}

bool DijkstraGraphSequenceBase::QueueLimitsExceeded(size_t iter) {
    return_code_.queue_limit = q_.size() > queue_limit_;
    return_code_.iter_limit = iter > iter_limit_;
    return return_code_.status;
}

bool DijkstraGraphSequenceBase::RunDijkstra() {
    bool found_path = false;
    size_t iter = 0;
    QueueState cur_state;
    int ed = 0;
    while (q_.size() > 0 &&
            !QueueLimitsExceeded(iter) &&
            ed <= path_max_length_ &&
            updates_ < gap_cfg_.updates_limit) {
        cur_state = q_.begin()->second;
        ed = visited_[cur_state];
        ++ iter;
        q_.erase(q_.begin());
        if (visited_.count(end_qstate_) > 0) {
            found_path = true;
        }
        if (IsEndPosition(cur_state)) {
            end_qstate_ = cur_state;
            return true;
        }
        for (const EdgeId &e : g_.OutgoingEdges(g_.EdgeEnd(cur_state.gs.e))) {
            found_path = AddState(cur_state, e, ed);
            if (!gap_cfg_.find_shortest_path && found_path) return true;
        }
        if (!gap_cfg_.find_shortest_path && found_path) return true;
    }
    return found_path;
}

void DijkstraGraphSequenceBase::CloseGap() {
    bool found_path = RunDijkstra();
    DEBUG("updates=" << updates_)

    if (!found_path) {
        return_code_.no_path = true;
    }
    if (found_path) {
        QueueState state(end_qstate_);
        while (!state.empty()) {
            min_score_ = visited_[end_qstate_];
            int start_edge = prev_states_[state].i;
            int end_edge =  state.i;
            mapping_path_.push_back(state.gs.e,
                                    omnigraph::MappingRange(Range(start_edge, end_edge),
                                            Range(state.gs.start_pos, state.gs.end_pos) ));
            state = prev_states_[state];
        }
        mapping_path_.reverse();
    }
    return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool DijkstraGapFiller::AddState(const QueueState &cur_state, EdgeId e, int ed) {
    if (reachable_vertex_.size() == 0 || reachable_vertex_.count(g_.EdgeEnd(cur_state.gs.e)) > 0) {
        if (reachable_vertex_.size() == 0 || reachable_vertex_.count(g_.EdgeEnd(e)) > 0) {
            GraphState next_state(e, 0, (int) g_.length(e) );
            AddNewEdge(next_state, cur_state, ed);
        }
        if (e == end_e_ && path_max_length_ - ed >= 0) {
            string seq_str = ss_.substr(cur_state.i);
            string edge_str = g_.EdgeNucls(e).Subseq(0, end_p_).str();
            int score = StringDistance(seq_str, edge_str, path_max_length_ - ed);
            if (score != numeric_limits<int>::max()) {
                path_max_length_ = min(path_max_length_, ed + score);
                QueueState state(GraphState(e, 0, end_p_), (int)ss_.size());
                Update(state, cur_state, ed + score);
                if (ed + score == path_max_length_) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool DijkstraGapFiller::IsEndPosition(const QueueState &cur_state) {
    if (cur_state.i == (int) ss_.size() &&
            cur_state.gs.e == end_qstate_.gs.e &&
            cur_state.gs.end_pos == end_qstate_.gs.end_pos) {
        return true;
    }
    return false;
}



bool DijkstraEndsReconstructor::AddState(const QueueState &cur_state, EdgeId e, int ed) {
    if (IsEndPosition(cur_state)) {
        return false;
    }
    GraphState next_state(e, 0, (int) g_.length(e));
    AddNewEdge(next_state, cur_state, ed);
    VERIFY(ss_.size() >= (size_t) cur_state.i)
    size_t remaining = ss_.size() - cur_state.i;
    if (g_.length(e) + g_.k() + path_max_length_ - ed > remaining && path_max_length_ - ed >= 0) {
        string seq_str = ss_.substr(cur_state.i);
        string edge_str = g_.EdgeNucls(e).str();
        int position = -1;
        int score = SHWDistance(seq_str, edge_str, path_max_length_ - ed, position);
        if (score != numeric_limits<int>::max()) {
            path_max_length_ = min(path_max_length_, ed + score);
            QueueState state(GraphState(e, 0, position + 1), (int) ss_.size());
            Update(state, cur_state, ed + score);
            if (ed + score == path_max_length_) {
                end_qstate_ = state;
                return true;
            }
        }
    }
    return false;
}

bool DijkstraEndsReconstructor::IsEndPosition(const QueueState &cur_state) {
    return (cur_state.i == (int)ss_.size());
}

} // namespace sensitive_aligner
