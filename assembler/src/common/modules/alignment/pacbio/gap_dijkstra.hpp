#include "assembly_graph/index/edge_multi_index.hpp"
#include "modules/alignment/edge_index_refiller.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"
#include "utils/perf/perfcounter.hpp"
#include "string_distance.hpp"

namespace pacbio {

struct GapClosingConfig {
    bool run_dijkstra; // default: true
    bool restore_ends; // default: false
    // dijkstra
    int max_vertex_in_gap;
    int max_gap_length;
    int queue_limit;
    int iteration_limit;
    bool find_shortest_path; // default: false
    bool restore_mapping; // default: false

    GapClosingConfig()
        :run_dijkstra(false), 
         restore_ends(false),
         max_vertex_in_gap(-1),
         max_gap_length(-1),
         queue_limit(-1),
         iteration_limit(-1), 
         find_shortest_path(false), 
         restore_mapping(false) {}
};


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

    bool empty()
    {
        return *this == QueueState();
    }
};


class DijkstraGraphSequenceBase {

protected:
        bool ShouldUpdateQueue(int seq_ind, int ed)
        {
            //return true;
            if (seq_ind == -1){
                return true;
            }
            VERIFY(seq_ind < (int) ss_.size())
            VERIFY(seq_ind >= 0)
            if (best_ed_[seq_ind] + 20 >= ed) {
                if (seq_ind != ss_.size() - 1) {
                    best_ed_[seq_ind] = min(best_ed_[seq_ind], ed);
                }  
                return true; 
            } else {
                return false;
            }
        }

        void Update(const QueueState &state, const QueueState &prev_state, int score)
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

        void AddNewEdge(const GraphState &gs, const QueueState &prev_state, int ed)
        {
            utils::perf_counter perf0;
            string tmp = g_.EdgeNucls(gs.e).str();
            string edge_str = tmp.substr(gs.start_pos, gs.end_pos - gs.start_pos);
            DEBUG("TIME.Substrs=" << perf0.time())
            if ( 0 == (int) edge_str.size()) {
                utils::perf_counter perf1;
                QueueState state(gs, prev_state.i);
                Update(state, prev_state,  ed);
                DEBUG("TIME.QueueUpdate=" << perf1.time())
                return;
            }
            utils::perf_counter perf2;
            int len = min( (int) g_.length(gs.e) - gs.start_pos + path_max_length_, (int) ss_.size() - (prev_state.i + 1) );
            string seq_str = ss_.substr(prev_state.i + 1, len);
            DEBUG("TIME.Substrs=" << perf2.time())
            vector<int> positions;
            vector<int> scores;
            if (path_max_length_ - ed >= 0) {
                if (path_max_length_ - ed >= (int) edge_str.size()) {
                    utils::perf_counter perf3;
                    QueueState state(gs, prev_state.i);
                    Update(state, prev_state,  ed + (int) edge_str.size());
                    DEBUG("TIME.QueueUpdate=" << perf3.time())
                }
                if (ss_.size() - (prev_state.i + 1) > 0) {
                    utils::perf_counter perf4;
                    SHWDistance(seq_str, edge_str, path_max_length_ - ed, positions, scores);
                    DEBUG("TIME.SHWDistance=" << perf4.time())
                    utils::perf_counter perf5;
                    for (int i = 0; i < positions.size(); ++ i) {
                        if (positions[i] >= 0 && scores[i] >= 0) {
                            QueueState state(gs, prev_state.i + 1 + positions[i]);
                            Update(state, prev_state, ed + scores[i]);
                        }
                    }
                    DEBUG("TIME.QueueUpdate=" << perf5.time())
                }
            }
            //INFO("Add edge: eid=" << gs.e.int_id() << " " << gs.start_pos << " " << gs.end_pos << " " << seq_str.size() << " " << q_.size() << " " << path_max_length_)
        }


        virtual bool AddState(const QueueState &cur_state, EdgeId e, int ed, int ind, utils::perf_counter &perf) = 0;

        virtual bool IsEndPosition(const QueueState &cur_state) = 0;

public:
        DijkstraGraphSequenceBase(const Graph &g, const GapClosingConfig &gap_cfg, string ss, EdgeId start_e, int start_p, int path_max_length)
                  :g_(g)
                   , gap_cfg_(gap_cfg)    
                   , ss_(ss)
                   , start_e_(start_e)
                   , start_p_(start_p)
                   , path_max_length_(path_max_length)
                   , queue_limit_(gap_cfg_.queue_limit)
                   , iter_limit_(gap_cfg_.iteration_limit){
            best_ed_.resize(ss_.size(), path_max_length_);
            AddNewEdge(GraphState(start_e_, start_p_, (int) g_.length(start_e_)), QueueState(), 0);
            min_score_ = -1;
        }

        void CloseGap() {
            bool found_path = false;
            int i = 0;
            utils::perf_counter perf;
            while (q_.size() > 0) {
                QueueState cur_state = q_.begin()->second;
                int ed = visited_[cur_state];
                //INFO("Queue edge=" << cur_state.gs.e.int_id() << " ed=" << ed);
                if (q_.size() > queue_limit_ || i > iter_limit_) {
                    //INFO("EdgeDijkstra: queue size is too big ed=" << ed << " q_.size=" << q_.size() << " i=" << i << " s_len=" << ss_.size() << " time=" << perf.time() )
                    if (visited_.count(end_qstate_) > 0){
                        found_path = true;
                        min_score_ = visited_[end_qstate_];
                    }
                    break;
                }
                int tid = omp_get_thread_num();
                //INFO("TID=" << tid << " Queue edge=" << cur_state.gs.e.int_id() <<  " " << cur_state.gs.start_pos <<  " " << cur_state.gs.end_pos << " " << cur_state.i  << " ed=" << ed << " pml=" << path_max_length_ );
                if (IsEndPosition(cur_state)) {
                    found_path = true;
                    break;   
                }
                if (ed > path_max_length_) {
                    //INFO("EdgeDijkstra: path not found ed=" << ed << " path_max_length=" << path_max_length_ << " q_.size=" << q_.size() << " i=" << i << " s_len=" << ss_.size() << " time=" << perf.time()  )
                    break;
                }
                i ++;
                q_.erase(q_.begin());
                for (const EdgeId &e: g_.OutgoingEdges(g_.EdgeEnd(cur_state.gs.e))) {
                    found_path = AddState(cur_state, e, ed, i, perf);
                    if (!gap_cfg_.find_shortest_path && found_path) break;
                }
                if (!gap_cfg_.find_shortest_path && found_path) break;
            }
            //INFO("TIME.Dijkstra=" << perf.time())
            if (found_path) {
                QueueState state = end_qstate_;
                while (!state.empty()) {
                    gap_path_.push_back(state.gs.e);
                    int tid = omp_get_thread_num();
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

        vector<EdgeId> GetPath() {
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
        set<pair<int, QueueState> > q_;
        map<QueueState, int> visited_;
        map<QueueState, QueueState> prev_states_;

        vector<EdgeId> gap_path_;

        omnigraph::MappingPath<debruijn_graph::EdgeId> mapping_path_;

        vector<int> best_ed_;

        const Graph &g_;
        const GapClosingConfig &gap_cfg_;
        const string ss_;
        EdgeId start_e_;
        const int start_p_;
        int path_max_length_;
        int min_score_;
        QueueState end_qstate_;

        const size_t queue_limit_;
        const size_t iter_limit_;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DijkstraGapFiller: public DijkstraGraphSequenceBase {  

private:

        virtual bool AddState(const QueueState &cur_state, EdgeId e, int ed, int ind, utils::perf_counter &perf) {
            bool found_path = false;
            if (reachable_vertex_.size() == 0 || reachable_vertex_.count(g_.EdgeEnd(cur_state.gs.e)) > 0) { 
                if (reachable_vertex_.size() == 0 || reachable_vertex_.count(g_.EdgeEnd(e)) > 0) {
                    GraphState next_state(e, 0, (int) g_.length(e) );
                    AddNewEdge(next_state, cur_state, ed);
                }
                if (e == end_e_ && path_max_length_ - ed >= 0 && cur_state.i + 1 < ss_.size()){
                    utils::perf_counter perf0;
                    string seq_str = ss_.substr(cur_state.i + 1, ss_.size() - (cur_state.i + 1));
                    string tmp = g_.EdgeNucls(e).str();
                    string edge_str = tmp.substr(0, end_p_);
                    DEBUG("TIME.Substrs=" << perf0.time())
                    utils::perf_counter perf1;
                    int score = NWDistance(seq_str, edge_str, path_max_length_ - ed);
                    DEBUG("TIME.NWDistance=" << perf1.time())
                    if (score != -1) {
                        if (ed +  score < path_max_length_) {
                            //INFO("Update max_path_len2=" << ed + score);
                        }
                        utils::perf_counter perf2;
                        path_max_length_ = min(path_max_length_, ed + score);
                        QueueState state(GraphState(e, 0, end_p_), ss_.size() - 1);
                        Update(state, cur_state, ed + score);
                        DEBUG("TIME.QueueUpdate=" << perf2.time())
                        if (ed + score == path_max_length_) {
                             min_score_ = ed + score;
                             //INFO("++=Final ed2=" << ed + score<< " q_.size=" << q_.size() << " i=" << ind << " s_len=" << ss_.size() << " time=" << perf.time() )
                             found_path = true;
                        }
                    }
                }
            }
            return found_path;
        }

        virtual bool IsEndPosition(const QueueState &cur_state) {
            if (cur_state.i == end_qstate_.i && cur_state.gs.e == end_qstate_.gs.e && cur_state.gs.end_pos == end_qstate_.gs.end_pos) {
                return true;   
            }
            return false;
        }

public:
        DijkstraGapFiller(const Graph &g, const GapClosingConfig &gap_cfg, string ss, EdgeId start_e, EdgeId end_e,
                                  int start_p, int end_p, int path_max_length, 
                                  const map<VertexId, size_t> &reachable_vertex)
                  : DijkstraGraphSequenceBase(g, gap_cfg, ss, start_e, start_p, path_max_length)
                   , end_e_(end_e) , end_p_(end_p)
                   , reachable_vertex_(reachable_vertex){
            GraphState end_gstate(end_e_, 0, end_p_);
            end_qstate_ = QueueState(end_gstate, ss_.size() - 1);
            if (start_e_ == end_e_ && end_p_ - start_p_ > 0){
                string seq_str = ss_;
                string tmp = g_.EdgeNucls(start_e_).str();
                string edge_str = tmp.substr(start_p_, end_p_ - start_p_);
                int score = NWDistance(seq_str, edge_str, path_max_length_);
                if (score != -1) {
                    if (score < path_max_length_) {
                        //INFO("Update max_path_len=" << score);
                    }
                    path_max_length_ = min(path_max_length_, score);
                    QueueState state(GraphState(start_e_, start_p_, end_p_), ss_.size() - 1);
                    Update(state, QueueState(), score);
                    if (score == path_max_length_) {
                         min_score_ = score;
                         //INFO("++=Final ed=" << score<< " q_.size=" << q_.size() << " s_len=" << ss_.size())
                         //found_path = true;
                         end_qstate_ = state;
                    }
                }
            }
        }

    protected:
        EdgeId end_e_;
        const int end_p_;
        const map<VertexId, size_t> &reachable_vertex_;
};



class DijkstraEndsReconstructor: public DijkstraGraphSequenceBase {

private:
        virtual bool AddState(const QueueState &cur_state, EdgeId e, int ed, int ind, utils::perf_counter &perf) {
            bool found_path = false;
            //DEBUG("end_pos=" << cur_state.gs.end_pos - (int) g_.length(cur_state.gs.e))
            GraphState next_state(e, 0, (int) g_.length(e));
            AddNewEdge(next_state, cur_state, ed);
            int remaining = ss_.size() - cur_state.i;
            if (g_.length(e) + g_.k() + path_max_length_ - ed > remaining && path_max_length_ - ed >= 0 && cur_state.i + 1 < ss_.size()){
                string seq_str = ss_.substr(cur_state.i + 1, ss_.size() - (cur_state.i + 1) );
                string tmp = g_.EdgeNucls(e).str();
                int position = -1;
                int score = SHWDistance2(seq_str, tmp, path_max_length_ - ed, position);
                if (score >= 0) {
                    if (ed +  score < path_max_length_) {
                        DEBUG("Update max_path_len=" << ed + score);
                    }
                    path_max_length_ = min(path_max_length_, ed + score);
                    QueueState state(GraphState(e, 0, position + 1), ss_.size() - 1);
                    Update(state, cur_state, ed + score);
                    if (ed + score == path_max_length_) {
                         min_score_ = ed + score;
                         DEBUG("++=Final ed1=" << ed + score);
                         DEBUG("EdgeDijkstra1: path was found ed=" << ed + score << " q_.size=" << q_.size() << " i=" << ind << " s_len=" << ss_.size() << " time=" << perf.time() )
                         found_path = true;
                         end_qstate_ = state;
                    }
                }
            }
            return found_path;
        }

        virtual bool IsEndPosition(const QueueState &cur_state) {
            if (cur_state.i == end_qstate_.i) {
                return true;   
            }
            return false;
        }
public:
        DijkstraEndsReconstructor(const Graph &g, const GapClosingConfig &gap_cfg, string ss, EdgeId start_e, int start_p, int path_max_length)
                  :DijkstraGraphSequenceBase(g, gap_cfg, ss, start_e, start_p, path_max_length) {
            end_qstate_ = QueueState();
            if (g_.length(start_e_) + g_.k() - start_p_ + path_max_length_ > ss_.size()){
                string seq_str = ss_;
                string tmp = g_.EdgeNucls(start_e_).str();
                string edge_str = tmp.substr(start_p_, g_.length(start_e_) + g_.k() - start_p_);
                int position = -1;
                int score = SHWDistance2(seq_str, edge_str, path_max_length_, position);
                if (score != -1) {
                    if (score < path_max_length_) {
                        DEBUG("Update max_path_len=" << score);
                    }
                    path_max_length_ = min(path_max_length_, score);
                    QueueState state(GraphState(start_e_, start_p_, start_p_ + position + 1), ss_.size() - 1);
                    Update(state, QueueState(), score);
                    if (score == path_max_length_) {
                         min_score_ = score;
                         DEBUG("++=Final2 ed=" << score);
                         DEBUG("EdgeDijkstra2: path was found ed=" << score << " q_.size=" << q_.size() << " s_len=" << ss_.size() )
                         //found_path = true;
                         end_qstate_ = state;
                    }
                }
            }

        }    
};

}