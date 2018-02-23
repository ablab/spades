#include "assembly_graph/index/edge_multi_index.hpp"
#include "modules/alignment/edge_index_refiller.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"
#include "utils/perf/perfcounter.hpp"
#include "string_distance.hpp"

namespace pacbio {


struct GraphState {
    EdgeId e;
    int start_pos;
    int end_pos;

    GraphState(EdgeId e_, int start_pos_, int end_pos_)
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
        return (this->e < state.e || (this->e == state.e && this->start_pos < state.start_pos) 
                                  || (this->e == state.e && this->start_pos == state.start_pos && this->end_pos < state.end_pos));
    }
};


struct QueueState {
    GraphState gs;
    int i;

    QueueState()
            : gs(EdgeId(), -1, -1), i(-1)
    {}

    QueueState(EdgeId &e_, int start_pos_, int end_pos_, int i_)
                : gs(e_, start_pos_, end_pos_), i(i_)
    {}

    QueueState(GraphState gs_, int i_)
                : gs(gs_), i(i_)
    {}

    QueueState(const QueueState &state):
                gs(state.gs), i(state.i)
    {}

    bool operator == (const QueueState &state) const {
        return (this->gs == state.gs && this->i == state.i);
    }

    bool operator != (const QueueState &state) const {
        return (this->gs != state.gs || this->i != state.i);
    }

    bool operator < (const QueueState &state) const {
        return (this->i < state.i || (this->i == state.i && this->gs < state.gs) );
    }

    bool isempty()
    {
        return *this == QueueState();
    }
};


class DijkstraGraphSequenceVirtual {

protected:
        bool ShouldUpdateQueue(int seq_ind, int ed)
        {
            if (seq_ind == -1){
                return true;
            }
            VERIFY(seq_ind < (int) ss_.size())
            VERIFY(seq_ind >= 0)
            if (best_ed_[seq_ind] + 20 >= ed) {
                best_ed_[seq_ind] = min(best_ed_[seq_ind], ed);  
                if (seq_ind == ss_.size() - 1) {
                    path_max_length_ = min(path_max_length_, best_ed_[seq_ind] + 20);
                }
                return true; 
            } else {
                return false;
            }
        }

        void Update(QueueState &state, const QueueState &prev_state, int score, int ed, bool is_end_position=false)
        {
            if (visited_.count(state) > 0) {
                if (dist_[state] >= ed){
                    q_.erase(make_pair(visited_[state], state));
                    visited_[state] = score;
                    dist_[state] = ed;
                    prev_states_[state] = prev_state;
                    if (!is_end_position && ShouldUpdateQueue(state.i, ed)){
                        q_.insert(make_pair(score, state));
                    }
                }
            } else{
                if (ShouldUpdateQueue(state.i, ed)){
                    visited_.insert(make_pair(state, score));
                    dist_.insert(make_pair(state, ed));
                    prev_states_.insert(make_pair(state, prev_state));
                    if (!is_end_position) {
                        q_.insert(make_pair(score, state));
                    }
                }
            }
        }

        void AddNewEdge(GraphState gs, QueueState prev_state, int ed)
        {
            int end = min( (int) g_.length(gs.e) - gs.start_pos + path_max_length_, (int) ss_.size() - (prev_state.i + 1) );
            string seq_str = ss_.substr(prev_state.i + 1, end);
            string tmp = g_.EdgeNucls(gs.e).str();
            string edge_str = tmp.substr(gs.start_pos, gs.end_pos - gs.start_pos);

            vector<int> positions;
            vector<int> scores;
            if (path_max_length_ - ed >= 0) {
                if (path_max_length_ - ed >= (int) edge_str.size()) {
                    QueueState state(gs, prev_state.i);
                    Update(state, prev_state,  ed + (int) edge_str.size(), ed + (int) edge_str.size());
                }
                if (ss_.size() - (prev_state.i + 1) > 0) {
                    SHWDistance(seq_str, edge_str, path_max_length_ - ed, positions, scores);
                    for (int i = 0; i < positions.size(); ++ i) {
                        if (positions[i] >= 0 && scores[i] >= 0) {
                            QueueState state(gs, prev_state.i + 1 + positions[i]);
                            Update(state, prev_state, ed + scores[i], ed + scores[i]);
                        }
                    }
                }
            }
        }


        virtual bool AddState(QueueState &cur_state, const EdgeId &e, int ed, int ind, utils::perf_counter &perf) = 0;

public:
        DijkstraGraphSequenceVirtual(const Graph &g, const string &ss, EdgeId start_e, const int start_p, const int path_max_length)
                  :g_(g), ss_(ss)
                   , start_e_(start_e)
                   , start_p_(start_p)
                   , path_max_length_(path_max_length){
            best_ed_.resize(ss_.size());
            for (size_t i = 0; i < best_ed_.size(); ++ i){
                best_ed_[i] = path_max_length_;
            }
            AddNewEdge(GraphState(start_e_, start_p_, g_.length(start_e_)), QueueState(), 0);
            min_score_ = -1;
        }

        void CloseGap() {
            bool found_path = false;
            int i = 0;
            utils::perf_counter perf;
            while (q_.size() > 0) {
                QueueState cur_state = q_.begin()->second;
                int ed = dist_[cur_state];
                //DEBUG("Queue edge=" << cur_state.gs.e.int_id() << " ed=" << ed);
                if (q_.size() > 1000000 || i > 1000000) {
                    DEBUG("EdgeDijkstra: queue size is too big ed=" << ed << " q_.size=" << q_.size() << " i=" << i << " s_len=" << ss_.size() << " time=" << perf.time() )
                    if (visited_.count(end_qstate_) > 0){
                        found_path = true;
                        min_score_ = dist_[end_qstate_];
                    }
                    break;
                }
                if (ed > path_max_length_) {
                    DEBUG("EdgeDijkstra: path not found ed=" << ed << " q_.size=" << q_.size() << " i=" << i << " s_len=" << ss_.size() << " time=" << perf.time()  )
                    break;
                }
                i ++;
                q_.erase(q_.begin());
                for (const EdgeId &e: g_.OutgoingEdges(g_.EdgeEnd(cur_state.gs.e))) {
                    found_path = AddState(cur_state, e, ed, i, perf);
                    if (found_path) break;
                }
                if (found_path) break;
            }

            if (found_path) {
                QueueState state = end_qstate_;
                while (!state.isempty()) {
                    gap_path_.push_back(state.gs.e);
                    if (prev_states_[state].i > state.i) {
                        INFO("start=" << prev_states_[state].i << " end=" << state.i)
                    }
                    if (state.gs.start_pos > state.gs.end_pos) {
                        INFO("start2=" << state.gs.start_pos << " end2=" << state.gs.end_pos)
                    }
                    mapping_path_.push_back(state.gs.e, omnigraph::MappingRange(Range((int) max(0, prev_states_[state].i), state.i), 
                                                                                Range(state.gs.start_pos, state.gs.end_pos) ));
                    state = prev_states_[state];
                }
                std::reverse(gap_path_.begin(), gap_path_.end());
                mapping_path_.reverse();
            }
            return;
        }

        std::vector<EdgeId> GetPath() {
           return gap_path_;
        }

        omnigraph::MappingPath<debruijn_graph::EdgeId> GetMappingPath() {
           return mapping_path_;
        }

        int GetEditDistance() {
            return min_score_;
        }       

        int GetPathEndPosition() {
            DEBUG("End position edge=" << end_qstate_.gs.e.int_id() << " end_pos=" << end_qstate_.gs.end_pos << " seq_pos=" << end_qstate_.i << " s_len=" << ss_.size())
            return end_qstate_.gs.end_pos;
        }  

        int GetSeqEndPosition() {
            DEBUG("End position edge=" << end_qstate_.gs.e.int_id() << " end_pos=" << end_qstate_.gs.end_pos << " seq_pos=" << end_qstate_.i << " s_len=" << ss_.size())
            return end_qstate_.i;
        }    

protected:
        std::set<std::pair<int, QueueState> > q_;
        std::map<QueueState, int> visited_;
        std::map<QueueState, int> dist_;
        std::map<QueueState, QueueState> prev_states_;

        std::vector<EdgeId> gap_path_;

        omnigraph::MappingPath<debruijn_graph::EdgeId> mapping_path_;

        std::vector<int> best_ed_;

        const Graph &g_;
        const string &ss_;
        EdgeId start_e_;
        const int start_p_;
        int path_max_length_;
        int min_score_;
        QueueState end_qstate_;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DijkstraGapFiller: public DijkstraGraphSequenceVirtual {  

private:

        virtual bool AddState(QueueState &cur_state, const EdgeId &e, int ed, int ind, utils::perf_counter &perf) {
            bool found_path = false;
            if (reachable_vertex_.count(g_.EdgeEnd(cur_state.gs.e)) > 0) { 
                if (reachable_vertex_.count(g_.EdgeEnd(e)) > 0) {
                    GraphState next_state(e, 0, g_.length(e));
                    AddNewEdge(next_state, cur_state, ed);
                }
                if (e == end_e_ && path_max_length_ - ed >= 0){
                    string seq_str = ss_.substr(cur_state.i + 1, ss_.size() - cur_state.i - 1);
                    string tmp = g_.EdgeNucls(e).str();
                    string edge_str = tmp.substr(0, end_p_);
                    int score = NWDistance(seq_str, edge_str, path_max_length_ - ed);
                    if (score != -1) {
                        if (ed +  score < path_max_length_) {
                            DEBUG("Update max_path_len=" << ed + score);
                        }
                        path_max_length_ = min(path_max_length_, ed + score);
                        QueueState state(GraphState(e, 0, end_p_), ss_.size() - 1);
                        Update(state, cur_state, ed + score, ed + score);
                        if (ed + score == path_max_length_) {
                             min_score_ = ed + score;
                             DEBUG("++=Final ed=" << ed + score<< " q_.size=" << q_.size() << " i=" << ind << " s_len=" << ss_.size() << " time=" << perf.time() )
                             found_path = true;
                        }
                    }
                    //DEBUG("TIME.QueueUpdate=" << perf.time());
                }
            }
            return found_path;
        }

public:
        DijkstraGapFiller(const Graph &g, const string &ss, EdgeId start_e, EdgeId end_e,
                                  const int start_p, const int end_p, const int path_max_length, 
                                  const std::map<VertexId, size_t> &reachable_vertex)
                  : DijkstraGraphSequenceVirtual(g, ss, start_e, start_p, path_max_length)
                   , end_e_(end_e) , end_p_(end_p)
                   , reachable_vertex_(reachable_vertex){
            GraphState end_gstate(end_e_, 0, end_p_);
            end_qstate_ = QueueState(end_gstate, ss_.size() - 1);
            bool found_path = false;
            if (start_e_ == end_e_ && end_p_ - start_p_ > 0){
                string seq_str = ss_;
                string tmp = g_.EdgeNucls(start_e_).str();
                string edge_str = tmp.substr(start_p_, end_p_ - start_p_);
                int score = NWDistance(seq_str, edge_str, path_max_length_);
                if (score != -1) {
                    if (score < path_max_length_) {
                        DEBUG("Update max_path_len=" << score);
                    }
                    path_max_length_ = min(path_max_length_, score);
                    QueueState state(GraphState(start_e_, start_p_, end_p_), ss_.size() - 1);
                    Update(state, QueueState(), score, score);
                    if (score == path_max_length_) {
                         min_score_ = score;
                         DEBUG("++=Final ed=" << score<< " q_.size=" << q_.size() << " s_len=" << ss_.size())
                         found_path = true;
                         end_qstate_ = state;
                    }
                }
            }
            if (found_path) {
                QueueState state = end_qstate_;
                while (!state.isempty()) {
                    gap_path_.push_back(state.gs.e);
                    mapping_path_.push_back(state.gs.e, omnigraph::MappingRange(Range( (int) max(0, prev_states_[state].i), state.i), 
                                                                                Range(state.gs.start_pos, state.gs.end_pos) ));
                    state = prev_states_[state];
                }
                std::reverse(gap_path_.begin(), gap_path_.end());
                mapping_path_.reverse();
            }
        }

    protected:
        EdgeId end_e_;
        const int end_p_;
        const std::map<VertexId, size_t> &reachable_vertex_;
};



class DijkstraEndsReconstructor: public DijkstraGraphSequenceVirtual {

private:
        virtual bool AddState(QueueState &cur_state, const EdgeId &e, int ed, int ind, utils::perf_counter &perf) {
            bool found_path = false;
            //DEBUG("end_pos=" << cur_state.gs.end_pos - (int) g_.length(cur_state.gs.e))
            GraphState next_state(e, max(0, cur_state.gs.end_pos - (int) g_.length(cur_state.gs.e) ), g_.length(e));
            AddNewEdge(next_state, cur_state, ed);
            int remaining = ss_.size() - cur_state.i;
            if (g_.length(e) + path_max_length_ > remaining && path_max_length_ - ed >= 0){
                string seq_str = ss_.substr(cur_state.i + 1, ss_.size() - cur_state.i - 1 );
                string tmp = g_.EdgeNucls(e).str();
                vector<int> positions;
                vector<int> scores;
                SHWDistance(tmp, seq_str, path_max_length_ - ed, positions, scores);
                if (scores.size() > 0) {
                    int score = scores[0];
                    if (ed +  score < path_max_length_) {
                        DEBUG("Update max_path_len=" << ed + score);
                    }
                    path_max_length_ = min(path_max_length_, ed + score);
                    QueueState state(GraphState(e, 0, positions[0] + 1), ss_.size() - 1);
                    Update(state, cur_state, ed + score, ed + score, true);
                    if (ed + score == path_max_length_) {
                         min_score_ = ed + score;
                         DEBUG("++=Final ed=" << ed + score);
                         DEBUG("EdgeDijkstra: path was found ed=" << ed + score << " q_.size=" << q_.size() << " i=" << ind << " s_len=" << ss_.size() << " time=" << perf.time() )
                         found_path = true;
                         end_qstate_ = state;
                    }
                }
            }
            return found_path;
        }

public:
        DijkstraEndsReconstructor(const Graph &g, const string &ss, EdgeId start_e, const int start_p, const int path_max_length)
                  :DijkstraGraphSequenceVirtual(g, ss, start_e, start_p, path_max_length) {
            end_qstate_ = QueueState();
            bool found_path = false;
            if (g_.length(start_e_) + g_.k() - start_p_ + path_max_length_ > ss_.size()){
                string seq_str = ss_;
                string tmp = g_.EdgeNucls(start_e_).str();
                string edge_str = tmp.substr(start_p_, g_.length(start_e_) + g_.k() - start_p_);
                vector<int> positions;
                vector<int> scores;
                SHWDistance(edge_str, seq_str, path_max_length_, positions, scores);
                if (scores.size() > 0) {
                    int score = scores[0];
                    if (score < path_max_length_) {
                        DEBUG("Update max_path_len=" << score);
                    }
                    path_max_length_ = min(path_max_length_, score);
                    QueueState state(GraphState(start_e_, start_p_, start_p_ + positions[0] + 1), ss_.size() - 1);
                    Update(state, QueueState(), score, score, true);
                    if (score == path_max_length_) {
                         min_score_ = score;
                         DEBUG("++=Final ed=" << score);
                         DEBUG("EdgeDijkstra: path was found ed=" << score << " q_.size=" << q_.size() << " s_len=" << ss_.size() )
                         found_path = true;
                         end_qstate_ = state;
                    }
                }
            }
            if (found_path) {
                QueueState state = end_qstate_;
                while (!state.isempty()) {
                    gap_path_.push_back(state.gs.e);
                    if (prev_states_[state].i > state.i) {
                        INFO("2start=" << prev_states_[state].i << " end=" << state.i)
                    }
                    if (state.gs.start_pos > state.gs.end_pos) {
                        INFO("2start2=" << state.gs.start_pos << " end2=" << state.gs.end_pos)
                    }
                    mapping_path_.push_back(state.gs.e, omnigraph::MappingRange(Range( (int) max(0, prev_states_[state].i), state.i), 
                                                                                Range(state.gs.start_pos, state.gs.end_pos) ));
                    state = prev_states_[state];
                }
                std::reverse(gap_path_.begin(), gap_path_.end());
                mapping_path_.reverse();
            }

        }    
};

}