#include "assembly_graph/index/edge_multi_index.hpp"
#include "modules/alignment/edge_index_refiller.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"
#include "utils/perf/perfcounter.hpp"

namespace pacbio {


inline void SHWDistance(string &target, string &query, int max_score, vector<int> &positions, vector<int> &scores) {
    if (query.size() == 0){
        for (int i = 0; i < min(max_score, (int) target.size()); ++ i) {
            positions.push_back(i);
            scores.push_back(i + 1);
        }
	   return;
    }
    VERIFY(target.size() > 0)
    edlib::EdlibEqualityPair additionalEqualities[36] = {{'U', 'T'}
                                        , {'R', 'A'}, {'R', 'G'}
                                        , {'Y', 'C'}, {'Y', 'T'}, {'Y', 'U'}
                                        , {'K', 'G'}, {'K', 'T'}, {'K', 'U'}
                                        , {'M', 'A'}, {'M', 'C'}
                                        , {'S', 'C'}, {'S', 'G'}
                                        , {'W', 'A'}, {'W', 'T'}, {'W', 'U'}
                                        , {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'B', 'U'}
                                        , {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'D', 'U'}
                                        , {'H', 'A'}, {'H', 'C'}, {'H', 'T'}, {'H', 'U'}
                                        , {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
                                        , {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'N', 'U'} };
    edlib::EdlibAlignResult result = edlib::edlibAlign(query.c_str(), (int) query.size(), target.c_str(), (int) target.size()
                                                   , edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_SHW_FULL, edlib::EDLIB_TASK_DISTANCE,
                                                                         additionalEqualities, 36));
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        positions.reserve(result.numLocations);
        scores.reserve(result.numLocations);
        for (int i = 0; i < result.numLocations; ++ i) {
            //INFO("Loc=" << result.endLocations[i] << " score=" << result.endScores[i]);
            if (result.endLocations[i] >= 0) {
                positions.push_back(result.endLocations[i]);
                scores.push_back(result.endScores[i]);
            }
        }
    }
    edlib::edlibFreeAlignResult(result);
}

inline void SHWDistance2(string &target, string &query, int max_score, vector<int> &positions, vector<int> &scores) {
    if (query.size() == 0){
        for (int i = 0; i < min(max_score, (int) target.size()); ++ i) {
            positions.push_back(i);
            scores.push_back(i + 1);
        }
       return;
    }
    VERIFY(target.size() > 0)
    edlib::EdlibEqualityPair additionalEqualities[36] = {{'U', 'T'}
                                        , {'R', 'A'}, {'R', 'G'}
                                        , {'Y', 'C'}, {'Y', 'T'}, {'Y', 'U'}
                                        , {'K', 'G'}, {'K', 'T'}, {'K', 'U'}
                                        , {'M', 'A'}, {'M', 'C'}
                                        , {'S', 'C'}, {'S', 'G'}
                                        , {'W', 'A'}, {'W', 'T'}, {'W', 'U'}
                                        , {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'B', 'U'}
                                        , {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'D', 'U'}
                                        , {'H', 'A'}, {'H', 'C'}, {'H', 'T'}, {'H', 'U'}
                                        , {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
                                        , {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'N', 'U'} };
    edlib::EdlibAlignResult result = edlib::edlibAlign(query.c_str(), (int) query.size(), target.c_str(), (int) target.size()
                                                   , edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_SHW, edlib::EDLIB_TASK_DISTANCE,
                                                                         additionalEqualities, 36));
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        positions.reserve(result.numLocations);
        scores.reserve(result.numLocations);
        for (int i = 0; i < result.numLocations; ++ i) {
            //INFO("Loc=" << result.endLocations[i] << " score=" << result.endScores[i]);
            if (result.endLocations[i] >= 0) {
                positions.push_back(result.endLocations[i]);
                scores.push_back(result.editDistance);
            }
        }
    }
    edlib::edlibFreeAlignResult(result);
}

inline int NWDistance(const string &a, const string &b, int max_score) {
    int a_len = (int) a.length();
    int b_len = (int) b.length();
    if (a_len == 0){
       if (b_len > max_score){
       	   return -1;
       } else{
           return b_len;
       }
    }
    if (b_len == 0){
       if (a_len > max_score){
           return -1;
       } else{
           return a_len;
       }
    }
    DEBUG(a_len << " " << b_len << " " << max_score);
    edlib::EdlibEqualityPair additionalEqualities[36] = {{'U', 'T'}
                                                , {'R', 'A'}, {'R', 'G'}
                                                , {'Y', 'C'}, {'Y', 'T'}, {'Y', 'U'}
                                                , {'K', 'G'}, {'K', 'T'}, {'K', 'U'}
                                                , {'M', 'A'}, {'M', 'C'}
                                                , {'S', 'C'}, {'S', 'G'}
                                                , {'W', 'A'}, {'W', 'T'}, {'W', 'U'}
                                                , {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'B', 'U'}
                                                , {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'D', 'U'}
                                                , {'H', 'A'}, {'H', 'C'}, {'H', 'T'}, {'H', 'U'}
                                                , {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
                                                , {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'N', 'U'} };

    edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len, b.c_str(), b_len
                                                   , edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_DISTANCE,
                                                                         additionalEqualities, 36));
    int score = -1;
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        score = result.editDistance;
    }
    edlib::edlibFreeAlignResult(result);
    return score;
}


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


class DijkstraGapFiller {

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

        void Update(QueueState &state, QueueState &prev_state, int score, int ed)
        {
            if (visited_.count(state) > 0) {
                if (dist_[state] >= ed){
                    q_.erase(make_pair(visited_[state], state));
                    visited_[state] = score;
                    dist_[state] = ed;
                    prev_states_[state] = prev_state;
                    if (ShouldUpdateQueue(state.i, ed)){
                        q_.insert(make_pair(score, state));
                    }
                }
            } else{
                if (ShouldUpdateQueue(state.i, ed)){
                    visited_.insert(make_pair(state, score));
                    dist_.insert(make_pair(state, ed));
                    prev_states_.insert(make_pair(state, prev_state));
                    q_.insert(make_pair(score, state));
                }
            }
        }        

        void AddNewEdge(GraphState gs, QueueState prev_state, int ed)
        {
            utils::perf_counter perf;
            int end = min( (int) g_.length(gs.e) - gs.start_pos + path_max_length_, (int) ss_.size() - (prev_state.i + 1) );
            // INFO("Add edge edge_id=" << gs.e.int_id() << " ed=" << ed << " seq_b=" << prev_state.i + 1 << " seq_e=" << prev_state.i + 1 + end  << " ss_size=" << ss_.size() 
            //                                                 << " edge_b=" << gs.start_pos << " edge_e=" << gs.end_pos)
            string seq_str = ss_.Subseq(prev_state.i + 1, prev_state.i + 1 + end).str();
            string tmp = g_.EdgeNucls(gs.e).str();
            string edge_str = tmp.substr(gs.start_pos, gs.end_pos - gs.start_pos);
            //INFO("TIME.Substrs=" << perf.time());
            
            vector<int> positions;
            vector<int> scores;
            if (path_max_length_ - ed >= 0) {
                perf.reset();
                //INFO("TIME.SHWDistance=" << perf.time());
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
                //INFO("TIME.QueueUpdate=" << perf.time());
            }
        }

    public:
        DijkstraGapFiller(const Graph &g, const Sequence &ss, EdgeId start_e, EdgeId end_e,
                                  const int start_p, const int end_p, const int path_max_length, 
                                  const std::map<VertexId, size_t> &reachable_vertex)
                  :g_(g), ss_(ss)
                   , start_e_(start_e), end_e_(end_e)
                   , start_p_(start_p), end_p_(end_p)
                   , path_max_length_(path_max_length) 
                   , reachable_vertex_(reachable_vertex){
            best_ed_.resize(ss_.size());
            for (size_t i = 0; i < best_ed_.size(); ++ i){
                best_ed_[i] = path_max_length_;
            }
            AddNewEdge(GraphState(start_e_, start_p_, g_.length(start_e_)), QueueState(), 0);
            min_score_ = -1;
        }

        void CloseGap() {
            GraphState end_gstate(end_e_, 0, end_p_);
            QueueState end_qstate(end_gstate, ss_.size() - 1);
            bool found_path = false;
            int i = 0;
            while (q_.size() > 0) {
                QueueState cur_state = q_.begin()->second;
                int ed = dist_[cur_state];
                //INFO("cur ed=" << ed << " edge=" << cur_state.gs.e.int_id())
                if (q_.size() > 1000000 || i > 1000000) {
                    INFO("queue size is too big")
                    if (visited_.count(end_qstate) > 0){
                        found_path = true;
                        min_score_ = dist_[end_qstate];
                    }
                    break;
                }
                if (ed > path_max_length_) {
                    INFO("No path found")
                    break;
                }
                if (cur_state == end_qstate && ed <= path_max_length_) {
                    min_score_ = ed;
                    INFO("++=Final ed=" << ed);
                    found_path = true;
                    break;
                }
                i ++;
                q_.erase(q_.begin());
                for (const EdgeId &e: g_.OutgoingEdges(g_.EdgeEnd(cur_state.gs.e))) {
                    if (reachable_vertex_.count(g_.EdgeEnd(cur_state.gs.e)) > 0) { 
                        if (reachable_vertex_.count(g_.EdgeEnd(e)) > 0) {
                            GraphState next_state(e, 0, g_.length(e));
                            AddNewEdge(next_state, cur_state, ed);
                        }
                        if (e == end_e_ && path_max_length_ - ed >= 0){
                            utils::perf_counter perf;
                            string seq_str = ss_.Subseq(cur_state.i + 1, ss_.size() ).str();
                            string tmp = g_.EdgeNucls(e).str();
                            string edge_str = tmp.substr(0, end_p_);
                            // INFO("End edge_id=" << e.int_id() << " ed=" << ed << " seq_b=" << cur_state.i + 1 << " seq_e=" << ss_.size() 
                            //                                 << " edge_b=" << 0 << " edge_e=" << end_p_
                            //                                 << " seq=" << seq_str << "\n edge=" << edge_str)
                            //INFO("TIME.Substrs=" << perf.time());
                            perf.reset();
                            int score = NWDistance(seq_str, edge_str, path_max_length_ - ed);
                            //INFO("TIME.NWDistance=" << perf.time());
                            perf.reset();
                            if (score != -1) {
                                if (ed +  score < path_max_length_) {
                                    INFO("Update max_path_len=" << ed + score);
                                }
                                path_max_length_ = min(path_max_length_, ed + score);
                                QueueState state(GraphState(e, 0, end_p_), ss_.size() - 1);
                                Update(state, cur_state, ed + score, ed + score);
				                if (ed + score == path_max_length_) {
                    			     min_score_ = ed + score;
                    			     INFO("++=Final ed=" << ed + score);
                    			     found_path = true;
                    			     break;
                		        }
                            }
                            //INFO("TIME.QueueUpdate=" << perf.time());
                        }
                    }
                }
		        if (found_path) break;
            }

            if (found_path) {
                QueueState state = end_qstate;
                end_qstate_ = end_qstate;
                while (!state.isempty()) {
                    gap_path_.push_back(state.gs.e);
                    state = prev_states_[state];
                }
                std::reverse(gap_path_.begin(), gap_path_.end());
            }
            return;
        }

        std::vector<EdgeId> GetPath() {
           return gap_path_;
        }

        int GetEditDistance() {
            return min_score_;
        }

        int GetEndQueueState() {
            return end_qstate_.gs.end_pos;
        }

    protected:
        std::set<std::pair<int, QueueState> > q_;
        std::map<QueueState, int> visited_;
        std::map<QueueState, int> dist_;
        std::map<QueueState, QueueState> prev_states_;

        std::vector<EdgeId> gap_path_;

        std::vector<int> best_ed_;

        const Graph &g_;
        const Sequence &ss_;
        EdgeId start_e_;
        EdgeId end_e_;
        const int start_p_;
        const int end_p_;
        int path_max_length_;
        const std::map<VertexId, size_t> &reachable_vertex_;
        int min_score_;
        QueueState end_qstate_;
};





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DijkstraEndsReconstructor {

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

        void Update(QueueState &state, QueueState &prev_state, int score, int ed, bool is_end_position=false)
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
            //INFO("Adding new state " << g_.int_id(gs.e) << " " << ed)
            utils::perf_counter perf;
            int end = min( (int) g_.length(gs.e) - gs.start_pos + path_max_length_, (int) ss_.size() - (prev_state.i + 1) );
            string seq_str = ss_.Subseq(prev_state.i + 1, prev_state.i + 1 + end).str();
            string tmp = g_.EdgeNucls(gs.e).str();
            string edge_str = tmp.substr(gs.start_pos, gs.end_pos - gs.start_pos);
            // INFO(" =>Add edge edge_id=" << gs.e.int_id() << " ed=" << ed << " seq_b=" << prev_state.i + 1 << " seq_e=" << prev_state.i + 1 + end  << " ss_size=" << ss_.size() 
            //                                                 << " edge_b=" << gs.start_pos << " edge_e=" << gs.end_pos)
            //INFO("TIME.Substrs=" << perf.time());
            vector<int> positions;
            vector<int> scores;
            if (path_max_length_ - ed >= 0) {
                perf.reset();
                //INFO("TIME.SHWDistance=" << perf.time());
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
                //INFO("TIME.QueueUpdate=" << perf.time());
            }
        }

    public:
        DijkstraEndsReconstructor(const Graph &g, const Sequence &ss, EdgeId start_e, const int start_p, const int path_max_length)
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
            if (q_.size() == 0) {
                INFO("EdgeDijkstra: queue size is too small q_.size=" << q_.size() << " start_e=" << start_e_.int_id() << " start_p=" << start_p_ << " e_len=" << g_.length(start_e_) << " s_len=" << ss_.size() << " time=" << perf.time() )
            } else {
                INFO("EdgeDijkstra: queue size begin q_.size=" << q_.size() << " start_e=" << start_e_.int_id() << " start_p=" << start_p_ << " e_len=" << g_.length(start_e_) << " s_len=" << ss_.size() << " time=" << perf.time() )   
            }
            while (q_.size() > 0) {
                QueueState cur_state = q_.begin()->second;
                int ed = dist_[cur_state];
                //INFO("Queue edge=" << cur_state.gs.e.int_id() << " ed=" << ed);
                if (q_.size() > 1000000 || i > 1000000) {
                    INFO("EdgeDijkstra: queue size is too big ed=" << ed << " q_.size=" << q_.size() << " i=" << i << " s_len=" << ss_.size() << " time=" << perf.time() )
                    if (visited_.count(end_qstate_) > 0){
                        found_path = true;
                        min_score_ = dist_[end_qstate_];
                    }
                    break;
                }
                if (ed > path_max_length_) {
                    INFO("EdgeDijkstra: path not found ed=" << ed << " q_.size=" << q_.size() << " i=" << i << " s_len=" << ss_.size() << " time=" << perf.time()  )
                    break;
                }
                i ++;
                q_.erase(q_.begin());
                for (const EdgeId &e: g_.OutgoingEdges(g_.EdgeEnd(cur_state.gs.e))) {
                    //INFO("end_pos=" << cur_state.gs.end_pos - (int) g_.length(cur_state.gs.e))
                    GraphState next_state(e, max(0, cur_state.gs.end_pos - (int) g_.length(cur_state.gs.e) ), g_.length(e));
                    AddNewEdge(next_state, cur_state, ed);
                    int remaining = ss_.size() - cur_state.i;
                    if (g_.length(e) + path_max_length_ > remaining && path_max_length_ - ed >= 0){
                        utils::perf_counter perf;
                        string seq_str = ss_.Subseq(cur_state.i + 1, ss_.size() ).str();
                        string tmp = g_.EdgeNucls(e).str();
                        //INFO("TIME.Substrs=" << perf.time());
                        perf.reset();
                        vector<int> positions;
                        vector<int> scores;
                        SHWDistance(tmp, seq_str, path_max_length_ - ed, positions, scores);
                        //INFO("TIME.NWDistance=" << perf.time());
                        perf.reset();
                        if (scores.size() > 0) {
                            int score = scores[0];
                            if (ed +  score < path_max_length_) {
                                INFO("Update max_path_len=" << ed + score);
                            }
                            path_max_length_ = min(path_max_length_, ed + score);
                            QueueState state(GraphState(e, 0, positions[0]), ss_.size() - 1);
                            Update(state, cur_state, ed + score, ed + score, true);
                            if (ed + score == path_max_length_) {
                                 min_score_ = ed + score;
                                 INFO("++=Final ed=" << ed + score);
                                 INFO("EdgeDijkstra: path was found ed=" << ed + score << " q_.size=" << q_.size() << " i=" << i << " s_len=" << ss_.size() << " time=" << perf.time() )
                                 found_path = true;
                                 end_qstate_ = state;
                                 break;
                            }
                        }
                    }
                }
                if (found_path) break;
            }
            if (q_.size() == 0) {
                INFO("EdgeDijkstra: after queue size is too small q_.size=" << q_.size() << " start_e=" << start_e_.int_id() << " start_p=" << start_p_ << " e_len=" << g_.length(start_e_) << " s_len=" << ss_.size() << " time=" << perf.time() )
            } else {
                INFO("EdgeDijkstra: after queue size begin q_.size=" << q_.size() << " start_e=" << start_e_.int_id() << " start_p=" << start_p_ << " e_len=" << g_.length(start_e_) << " s_len=" << ss_.size() << " time=" << perf.time() )   
            }

            if (found_path) {
                QueueState state = end_qstate_;
                while (!state.isempty()) {
                    gap_path_.push_back(state.gs.e);
                    state = prev_states_[state];
                }
                std::reverse(gap_path_.begin(), gap_path_.end());
                string tmp = "";
                for (int i = 1; i < gap_path_.size() - 1; ++i) {
                    size_t len = g_.length(gap_path_[i]);
                    string t = g_.EdgeNucls(gap_path_[i]).First(len).str();
                    tmp += t;
                }
                size_t len = g_.length(gap_path_[0]);
                string t = g_.EdgeNucls(gap_path_[0]).First(len).Last(len - start_p_).str();
                path_str_ += t + tmp + g_.EdgeNucls(gap_path_[gap_path_.size() - 1]).First(end_qstate_.gs.end_pos).str();
            }
            return;
        }

        std::vector<EdgeId> GetPath() {
           return gap_path_;
        }

        int GetEditDistance() {
            return min_score_;
        }

        int GetPathEndPosition() {
            INFO("End position edge=" << end_qstate_.gs.e.int_id() << " end_pos=" << end_qstate_.gs.end_pos << " seq_pos=" << end_qstate_.i << " s_len=" << ss_.size())
            return end_qstate_.gs.end_pos;
        }        

        string GetPathStr() {
            return path_str_;
        }

    protected:
        std::set<std::pair<int, QueueState> > q_;
        std::map<QueueState, int> visited_;
        std::map<QueueState, int> dist_;
        std::map<QueueState, QueueState> prev_states_;

        std::vector<EdgeId> gap_path_;

        std::vector<int> best_ed_;

        const Graph &g_;
        const Sequence &ss_;
        EdgeId start_e_;
        const int start_p_;
        int path_max_length_;
        int min_score_;
        QueueState end_qstate_;
        string path_str_;
};

//class EndsConstructor: public AbstractGapFiller {
//   virtual bool IsEnd(EdgeId e, QueueState &state){
//       int remaining = ss_.size() - state.i;
//       if (g_.length(e) - remaining > 0){
//           return true;
//       }
//       return false;
//   }
//public:
//   EndsConstructor(const Graph &g, const Sequence &ss, EdgeId start_e, EdgeId end_e,
//                                  const int start_p, const int end_p, const int path_max_length, 
//                                 const std::map<VertexId, size_t> &vs, const std::map<VertexId, size_t> &vs_b)
//           :AbstractGapFiller(g, ss, start_e, end_e, start_p, end_p, path_max_length, vs, vs_b)
//       {}
//
//};

}