//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "modules/alignment/pacbio/gap_dijkstra.hpp"

namespace sensitive_aligner {

using namespace std;

const set<string> DijkstraProteinGraph::stop_codons = {"TAA", "TAG", "TGA"};
const set<string> DijkstraProteinGraph::cstop_codons = {"TTA", "CTA", "TCA"};

void DijkstraReadGraph::AddState(const QueueState &state, int ed,const QueueState &prev, bool front) {
    utils::perf_counter pc;
    if (ed < path_max_length_) {
        if (visited_.count(state) == 0 || (visited_.count(state) > 0 && visited_[state] > ed)){
            ++ updates_;
            if (front) {
                q_.push_front(make_pair(ed, state));
            } else {
                q_.push_back(make_pair(ed, state));
            }
            visited_[state] = ed;
            prev_states_[state] = prev;
            if (best_ed_[state.i].second > ed) {
                best_ed_[state.i] = std::pair<QueueState, int> (state, ed);
            }
        }
    }
    state_time += pc.time();
}

void DijkstraReadGraph::AddNewStatesFrom(const QueueState &state, int ed) {
    // deletion
    utils::perf_counter pc;
    if (state.i + 1 <= (int) ss_.size()) {
        AddState(QueueState(state.gs, state.i + 1), ed + 1, state, false);
    }
    d_time += pc.time();
    pc.reset();

    if ((int) g_.length(state.gs.e) <= state.gs.end_pos) {
        for (const debruijn_graph::EdgeId &e : g_.OutgoingEdges(g_.EdgeEnd(state.gs.e))) {
            if (IsGoodEdge(e) > 0) {
                int new_pos = state.gs.end_pos - (int) g_.length(state.gs.e);
                AddState(QueueState(GraphState(e, new_pos, new_pos), state.i), ed, state, false);
            }
        }
    } else {
        // insertion
        AddState(QueueState(GraphState(state.gs.e, state.gs.end_pos, state.gs.end_pos + 1), state.i), ed + 1, state, false);
        i_time += pc.time();
        pc.reset();

        // match/mismatch
        if (state.i + 1 <= (int) ss_.size()) {
            string s = g_.EdgeNucls(state.gs.e).Subseq(state.gs.end_pos, state.gs.end_pos + 1).str();
            bool front = s[0] == ss_[state.i] ? true: false;
            int penalty = s[0] == ss_[state.i] ? 0: 1;
            AddState(QueueState(GraphState(state.gs.e, state.gs.end_pos, state.gs.end_pos + 1), state.i + 1), ed + penalty, state, front);
        }
        m_time += pc.time();
        pc.reset();
    }
}

template<class T, typename U>
bool DijkstraGraphSequenceBase<T,U>::QueueLimitsExceeded(size_t iter) {
    return_code_.queue_limit = q_.size() > queue_limit_;
    return_code_.iter_limit = iter > iter_limit_;
    return return_code_.status;
}

template<class T, typename U>
bool DijkstraGraphSequenceBase<T,U>::RunDijkstra() {
    bool found_path = false;
    size_t iter = 0;
    T cur_state;
    int ed = 0;
    while (q_.size() > 0 &&
            !QueueLimitsExceeded(iter) &&
            ed <= path_max_length_ &&
            updates_ < gap_cfg_.updates_limit) {
        cur_state = q_.begin()->second;
        ed = visited_[cur_state];
        PopFront();
        if (IsEndPosition(cur_state)) {
            end_qstate_ = cur_state;
            return true;
        }
        utils::perf_counter pc;
        if (ed == visited_[cur_state]) {
            AddNewStatesFrom(cur_state, ed);
            ++ iter;
        }
        newstates_time += pc.time();
    }
    return found_path;
}

bool DijkstraProteinEndsReconstructor::DetermineBestPrefix(bool found_path) {
    bool res = false;
    for (size_t i = 1; i < best_ed_.size(); ++ i) {
        if (i % 3 == 0 && best_ed_[i].second < (min_value_ - 4)*((int) i/3)) {
            end_qstate_ = best_ed_[i].first;
            res = true;
        }
    }
    return res || found_path;
}

bool DijkstraEndsReconstructor::DetermineBestPrefix(bool found_path) {
    double t = (double) path_max_length_/(double) ss_.size();
    bool res = false;
    if (!found_path) {
        for (size_t i = 1; i < best_ed_.size(); ++ i) {
            if (best_ed_[i].second < t*(int)i) {
                end_qstate_ = best_ed_[i].first;
                res = true;
            }
        }
    }
    return res || found_path;
}

bool DijkstraGapFiller::DetermineBestPrefix(bool found_path) {
    return found_path;
}

template<class T, typename U>
void DijkstraGraphSequenceBase<T,U>::CloseGap() {
    bool found_path = RunDijkstra();
    found_path = DetermineBestPrefix(found_path);
    if (!found_path) {
        return_code_.no_path = true;
    }
    if (found_path) {
        T state(end_qstate_);
        int end_on_seq =  state.i;
        int end_on_edge = state.gs.end_pos;
        while (!state.empty()) {
            min_score_ = visited_[end_qstate_];
            int start_on_seq = prev_states_[state].i;
            if (prev_states_[state].empty() ||  
                    prev_states_[state].gs.e != state.gs.e  ||  state.gs.end_pos < prev_states_[state].gs.end_pos)  {
                mapping_path_.push_back(state.gs.e,
                                        omnigraph::MappingRange(Range(start_on_seq, end_on_seq),
                                                                Range(state.gs.start_pos, end_on_edge) ));
                end_on_seq = start_on_seq;
                end_on_edge = prev_states_[state].gs.end_pos;
            }    
            state = prev_states_[state];
        }
        mapping_path_.reverse();
    }
    return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DijkstraProteinGraph::AddState(const ProteinQueueState &state, int ed,const ProteinQueueState &prev) {
    utils::perf_counter pc;
    if (ed < path_max_length_) {
        if (visited_.count(state) == 0 || (visited_.count(state) > 0 && visited_[state] > ed)){
            ++ updates_;
            if (visited_.count(state) > 0) {
                q_.erase(std::pair<int, ProteinQueueState>(visited_[state], state));
            }
            q_.insert(std::pair<int, ProteinQueueState>(ed, state));
            visited_[state] = ed;
            prev_states_[state] = prev;
            if (best_ed_[state.i].second > ed && best_ed_[state.i].first.offset.size() == 0) {
                best_ed_[state.i] = std::pair<ProteinQueueState, int> (state, ed);
            }
        }
    }
    state_time += pc.time();
}

void DijkstraProteinGraph::AddNewStatesFrom(const ProteinQueueState &state, int ed) {

    utils::perf_counter pc;
    // deletion
    if (state.offset.size() == 0) {
        if (state.i + 3 <= (int) ss_.size()) {
            AddState(ProteinQueueState(state.gs, state.i + 3, ""), ed + deletion_score(), state);
        }
        d_time += pc.time();
        pc.reset();
    }

    if (state.offset.size() == 3) {
        // insertion
        if (!StopCodon(state.offset)) {
            AddState(ProteinQueueState(state.gs, state.i, ""), ed + insertion_score(), state);
            i_time += pc.time();
            pc.reset();
        }

        //match/mismatch
        if (state.i + 3 <= (int) ss_.size()) {
            string aa = ss_.substr(state.i, 3);
            int score = mm_score(aa, state.offset);
            if (state.i + 3 == (int) ss_.size() || !StopCodon(state.offset)) {
                AddState(ProteinQueueState(state.gs, state.i + 3, ""), ed + score, state);
            }
        }
        m_time += pc.time();
        pc.reset();
    }

    if (state.offset.size() < 3) {
        if ((int) g_.length(state.gs.e) <= state.gs.end_pos) {
            for (const debruijn_graph::EdgeId &e : g_.OutgoingEdges(g_.EdgeEnd(state.gs.e))) {
                if (IsGoodEdge(e) > 0) {
                    int new_pos = state.gs.end_pos - (int) g_.length(state.gs.e);
                    AddState(ProteinQueueState(
                        GraphState(e, new_pos, new_pos), state.i, state.offset), ed, state);
                }
            }
        } else {
            string s = g_.EdgeNucls(state.gs.e).Subseq(state.gs.end_pos, state.gs.end_pos + 1).str();
            AddState(ProteinQueueState(GraphState(state.gs.e, state.gs.end_pos, state.gs.end_pos + 1), state.i, state.offset + s[0]), ed, state);
        }
    } 

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DijkstraGapFiller::IsEndPosition(const QueueState &cur_state) {
    if (cur_state.i == (int) ss_.size() &&
            cur_state.gs.e == end_qstate_.gs.e &&
            cur_state.gs.end_pos == end_qstate_.gs.end_pos) {
        return true;
    }
    return false;
}

bool DijkstraEndsReconstructor::IsEndPosition(const QueueState &cur_state) {
    return (cur_state.i == (int)ss_.size());
}

bool DijkstraProteinEndsReconstructor::IsEndPosition(const ProteinQueueState &cur_state) {
    return (cur_state.i == (int)ss_.size() && cur_state.offset.size() == 0);
}

template class DijkstraGraphSequenceBase<QueueState, std::deque<std::pair<int, QueueState> > >;
template class DijkstraGraphSequenceBase<ProteinQueueState, std::set<std::pair<int, ProteinQueueState> > >;

} // namespace sensitive_aligner