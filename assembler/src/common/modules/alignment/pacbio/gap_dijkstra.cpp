#include "gap_dijkstra.hpp"

namespace pacbio {


bool DijkstraGraphSequenceBase::ShouldUpdateQueue(int seq_ind, int ed)
{
    //return true;
    if (seq_ind == -1){
        return true;
    }
    VERIFY(seq_ind < (int) ss_.size())
    VERIFY(seq_ind >= 0)
    if (best_ed_[seq_ind] + gap_cfg_.penalty_interval >= ed) {
        if (seq_ind != (int) ss_.size() - 1) {
            best_ed_[seq_ind] = min(best_ed_[seq_ind], ed);
        }  
        return true; 
    } else {
        return false;
    }
}

void DijkstraGraphSequenceBase::Update(const QueueState &state, const QueueState &prev_state, int score)
{
    if (visited_.count(state) > 0) {
        if (visited_[state] >= score){
            q_.erase(make_pair(visited_[state], state));
            visited_[state] = score;
            prev_states_[state] = prev_state;
            if (ShouldUpdateQueue(state.i, score)){
                q_.insert(make_pair(score, state));
            }
        }
    } else{
        if (ShouldUpdateQueue(state.i, score)){
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
    int len = min( (int) g_.length(gs.e) - gs.start_pos + path_max_length_, (int) ss_.size() - (prev_state.i + 1) );
    string seq_str = ss_.substr(prev_state.i + 1, len);
    vector<int> positions;
    vector<int> scores;
    if (path_max_length_ - ed >= 0) {
        if (path_max_length_ - ed >= (int) edge_str.size()) {
            QueueState state(gs, prev_state.i);
            Update(state, prev_state,  ed + (int) edge_str.size());
        }
        if (ss_.size() - (prev_state.i + 1) > 0) {
            SHWDistance(seq_str, edge_str, path_max_length_ - ed, positions, scores);
            for (size_t i = 0; i < positions.size(); ++ i) {
                if (positions[i] >= 0 && scores[i] >= 0) {
                    QueueState state(gs, prev_state.i + 1 + positions[i]);
                    Update(state, prev_state, ed + scores[i]);
                }
            }
        }
    }
}

void DijkstraGraphSequenceBase::CloseGap() {
    bool found_path = false;
    size_t i = 0;
    while (q_.size() > 0) {
        QueueState cur_state = q_.begin()->second;
        int ed = visited_[cur_state];
        if (q_.size() > queue_limit_ || i > iter_limit_) {
            if (visited_.count(end_qstate_) > 0){
                found_path = true;
                min_score_ = visited_[end_qstate_];
            }
            break;
        }
        if (IsEndPosition(cur_state)) {
            found_path = true;
            break;   
        }
        if (ed > path_max_length_) {
            break;
        }
        i ++;
        q_.erase(q_.begin());
        for (const debruijn_graph::EdgeId &e: g_.OutgoingEdges(g_.EdgeEnd(cur_state.gs.e))) {
            found_path = AddState(cur_state, e, ed);
            if (!gap_cfg_.find_shortest_path && found_path) break;
        }
        if (!gap_cfg_.find_shortest_path && found_path) break;
    }
    if (found_path) {
        QueueState state = end_qstate_;
        while (!state.empty()) {
            gap_path_.push_back(state.gs.e);
            int start_edge = (int) max(0, (int) prev_states_[state].i);
            int end_edge = max(0, (int) state.i);
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
        if (e == end_e_ && path_max_length_ - ed >= 0 && cur_state.i + 1 < (int) ss_.size()){
            string seq_str = ss_.substr(cur_state.i + 1);
            string tmp = g_.EdgeNucls(e).str();
            string edge_str = tmp.substr(0, end_p_);
            int score = NWDistance(seq_str, edge_str, path_max_length_ - ed);
            if (score != -1) {
                path_max_length_ = min(path_max_length_, ed + score);
                QueueState state(GraphState(e, 0, end_p_), (int)ss_.size() - 1);
                Update(state, cur_state, ed + score);
                if (ed + score == path_max_length_) {
                     min_score_ = ed + score;
                     found_path = true;
                }
            }
        }
    }
    return found_path;
}


bool DijkstraEndsReconstructor::AddState(const QueueState &cur_state, debruijn_graph::EdgeId e, int ed) {
    bool found_path = false;
    GraphState next_state(e, 0, (int) g_.length(e));
    AddNewEdge(next_state, cur_state, ed);
    int remaining = (int) ss_.size() - cur_state.i;
    if ((int) g_.length(e) + (int) g_.k() + path_max_length_ - ed > remaining && path_max_length_ - ed >= 0 && cur_state.i + 1 < (int) ss_.size()){
        string seq_str = ss_.substr(cur_state.i + 1);
        string tmp = g_.EdgeNucls(e).str();
        int position = -1;
        int score = SHWDistance2(seq_str, tmp, path_max_length_ - ed, position);
        if (score >= 0) {
            path_max_length_ = min(path_max_length_, ed + score);
            QueueState state(GraphState(e, 0, position + 1), (int) ss_.size() - 1);
            Update(state, cur_state, ed + score);
            if (ed + score == path_max_length_) {
                 min_score_ = ed + score;
                 DEBUG("++=Final ed1=" << ed + score);
                 DEBUG("EdgeDijkstra1: path was found ed=" << ed + score << " q_.size=" << q_.size() << " s_len=" << ss_.size() )
                 found_path = true;
                 end_qstate_ = state;
            }
        }
    }
    return found_path;
}


}