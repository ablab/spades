#include "modules/alignment/pacbio/gap_dijkstra.hpp"

namespace gap_dijkstra {


bool DijkstraGraphSequenceBase::ShouldUpdateQueue(int seq_ind, int ed)
{
    if (seq_ind == ss_.size() ) {
        if (ed <= path_max_length_) {
            return true;
        } else {
            return false;
        }
    }
    VERIFY(seq_ind < (int) ss_.size())
    VERIFY(seq_ind >= 0)
    if (best_ed_[seq_ind] + gap_cfg_.penalty_interval >= ed) {
        best_ed_[seq_ind] = min(best_ed_[seq_ind], ed);
        return true;
    } else {
        return false;
    }
}

void DijkstraGraphSequenceBase::Update(const QueueState &state, const QueueState &prev_state, int score)
{
    if (visited_.count(state) > 0) {
        if (visited_[state] >= score) {
            q_.erase(make_pair(visited_[state], state));
            visited_[state] = score;
            prev_states_[state] = prev_state;
            if (ShouldUpdateQueue(state.i, score)) {
                q_.insert(make_pair(score, state));
            }
        }
    } else {
        if (ShouldUpdateQueue(state.i, score)) {
            visited_.insert(make_pair(state, score));
            prev_states_.insert(make_pair(state, prev_state));
            q_.insert(make_pair(score, state));
        }
    }
}

void DijkstraGraphSequenceBase::AddNewEdge(const GraphState &gs, const QueueState &prev_state, int ed)
{
    string tmp = g_.EdgeNucls(gs.e).str();
    string edge_str = tmp.substr(gs.start_pos, gs.end_pos - gs.start_pos);
    if ( 0 == (int) edge_str.size()) {
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
            SHWDistanceFull(seq_str, edge_str, path_max_length_ - ed, positions, scores);
            for (size_t i = 0; i < positions.size(); ++ i) {
                if (positions[i] >= 0 && scores[i] >= 0) {
                    QueueState state(gs, prev_state.i + positions[i] + 1);
                    Update(state, prev_state, ed + scores[i]);
                }
            }
        }
    }
}

bool DijkstraGraphSequenceBase::QueueLimitsExceeded(size_t iter) {
    if (q_.size() > queue_limit_) {
        return_code_ += DijkstraReturnCode::QUEUE_LIMIT;
    }
    if (iter > iter_limit_) {
        return_code_ += DijkstraReturnCode::ITERATION_LIMIT;
    }
    return q_.size() > queue_limit_ || iter > iter_limit_;
}

bool DijkstraGraphSequenceBase::RunDijkstra() {
    bool found_path = false;
    size_t iter = 0;
    QueueState cur_state;
    int ed = 0;
    while (q_.size() > 0 
           && !QueueLimitsExceeded(iter)
           && ed <= path_max_length_) 
    {
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
        for (const debruijn_graph::EdgeId &e : g_.OutgoingEdges(g_.EdgeEnd(cur_state.gs.e))) {
            found_path = AddState(cur_state, e, ed);
            if (!gap_cfg_.find_shortest_path && found_path) return true;
        }
        if (!gap_cfg_.find_shortest_path && found_path) return true;
    }
    return found_path;
}

void DijkstraGraphSequenceBase::CloseGap() {
    bool found_path = RunDijkstra();
    if (!found_path) {
        return_code_ += DijkstraReturnCode::NO_PATH;
    }
    if (found_path) {
        QueueState state(end_qstate_);
        while (!state.empty()) {
            gap_path_.push_back(state.gs.e);
            min_score_ = visited_[end_qstate_];
            int start_edge = prev_states_[state].i;
            int end_edge =  state.i;
            mapping_path_.push_back(state.gs.e, omnigraph::MappingRange(Range(start_edge, end_edge),
                                    Range(state.gs.start_pos, state.gs.end_pos) ));
            state = prev_states_[state];
        }
        reverse(gap_path_.begin(), gap_path_.end());
        mapping_path_.reverse();
    }
    return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool DijkstraGapFiller::AddState(const QueueState &cur_state, debruijn_graph::EdgeId e, int ed) {
    bool found_path = false;
    if (reachable_vertex_.size() == 0 || reachable_vertex_.count(g_.EdgeEnd(cur_state.gs.e)) > 0) {
        if (reachable_vertex_.size() == 0 || reachable_vertex_.count(g_.EdgeEnd(e)) > 0) {
            GraphState next_state(e, 0, (int) g_.length(e) );
            AddNewEdge(next_state, cur_state, ed);
        }
        if (e == end_e_ && path_max_length_ - ed >= 0) {
            string seq_str = ss_.substr(cur_state.i);
            string tmp = g_.EdgeNucls(e).str();
            string edge_str = tmp.substr(0, end_p_);
            int score = StringDistance(seq_str, edge_str, path_max_length_ - ed);
            if (score != std::numeric_limits<int>::max()) {
                path_max_length_ = min(path_max_length_, ed + score);
                QueueState state(GraphState(e, 0, end_p_), (int)ss_.size());
                Update(state, cur_state, ed + score);
                if (ed + score == path_max_length_) {
                    found_path = true;
                }
            }
        }
    }
    return found_path;
}

bool DijkstraGapFiller::IsEndPosition(const QueueState &cur_state) {
    if (cur_state.i == ss_.size() && cur_state.gs.e == end_qstate_.gs.e && cur_state.gs.end_pos == end_qstate_.gs.end_pos) {
        return true;
    }
    return false;
}



bool DijkstraEndsReconstructor::AddState(const QueueState &cur_state, debruijn_graph::EdgeId e, int ed) {
    bool found_path = false;
    if (!IsEndPosition(cur_state)) {
        GraphState next_state(e, 0, (int) g_.length(e));
        AddNewEdge(next_state, cur_state, ed);
        int remaining = (int) ss_.size() - cur_state.i;
        if ((int) g_.length(e) + (int) g_.k() + path_max_length_ - ed > remaining && path_max_length_ - ed >= 0) {
            string seq_str = ss_.substr(cur_state.i);
            string tmp = g_.EdgeNucls(e).str();
            int position = -1;
            int score = SHWDistance(seq_str, tmp, path_max_length_ - ed, position);
            if (score != std::numeric_limits<int>::max()) {
                path_max_length_ = min(path_max_length_, ed + score);
                QueueState state(GraphState(e, 0, position + 1), (int) ss_.size());
                Update(state, cur_state, ed + score);
                if (ed + score == path_max_length_) {
                    found_path = true;
                    end_qstate_ = state;
                }
            }
        }
    }
    return found_path;
}

bool DijkstraEndsReconstructor::IsEndPosition(const QueueState &cur_state) {
    if (cur_state.i == ss_.size()) {
        return true;
    }
    return false;
}

} // namespace gap_dijsktra
