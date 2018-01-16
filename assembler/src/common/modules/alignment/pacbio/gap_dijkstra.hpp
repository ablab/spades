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
    edlib::EdlibAlignResult result = edlib::edlibAlign(query.c_str(), query.size(), target.c_str(), target.size()
                                                   , edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_SHW_FULL, edlib::EDLIB_TASK_DISTANCE,
                                                                         additionalEqualities, 36));
    int score = -1;
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        score = result.editDistance;
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


class AbstractGapFiller {

        void Update(QueueState &state, QueueState &prevState, int score, int ed)
        {
            //if (g_.int_id(end_e_) == g_.int_id(state.gs.e) ) {
            //    INFO("  Add edge to q " << g_.int_id(state.gs.e) << " ed=" << score << " pos=" << state.i)
            //}
            if (visited.count(state) > 0) {
                if (dist[state] > ed && best_ed[state.i] + 10 > ed){
                    q.erase(make_pair(visited[state], state));
                    visited[state] = score;
                    dist[state] = ed;
                    best_ed[state.i] = ed;
                    prevStates[state] = prevState;
                    q.insert(make_pair(score, state));
                }
            } else{
                visited.insert(make_pair(state, score));
                best_ed[state.i] = ed;
                dist.insert(make_pair(state, ed));
                prevStates.insert(make_pair(state, prevState));
                q.insert(make_pair(score, state));
            }
        }

        void AddNewEdge(GraphState gs, QueueState prevState, int ed)
        {
            //INFO("Adding new state " << g_.int_id(gs.e) << " " << ed)
            utils::perf_counter perf;
            int end = min( (int) g_.length(gs.e) - gs.start_pos + path_max_length_, (int) ss_.size() - (prevState.i + 1) );
            string sstr = ss_.Subseq(prevState.i + 1, prevState.i + 1 + end).str();
            string tmp = g_.EdgeNucls(gs.e).str();
            string edgestr = tmp.substr(gs.start_pos, gs.end_pos - gs.start_pos);
            //INFO("TIME.Substrs=" << perf.time());
            vector<int> positions;
            vector<int> scores;
            if (path_max_length_ - ed >= 0){
                perf.reset();
                SHWDistance(sstr, edgestr, path_max_length_ - ed, positions, scores);
                //INFO("TIME.SHWDistance=" << perf.time());
                perf.reset();
                if (path_max_length_ - ed >= edgestr.size()) {
                    QueueState state(gs, prevState.i);
                    int remaining = ss_.size() - prevState.i;
                    Update(state, prevState,  ed + edgestr.size(), ed + edgestr.size());
                }
                for (int i = 0; i < positions.size(); ++ i) {
                    if (positions[i] >= 0 && scores[i] >= 0) {
                        QueueState state(gs, prevState.i + 1 + positions[i]);
                        int remaining = ss_.size() - (prevState.i + 1 + positions[i]);
                        Update(state, prevState, ed + scores[i], ed + scores[i]);
                    }
                }
                //INFO("TIME.QueueUpdate=" << perf.time());
            }
        }

        virtual bool IsEnd(EdgeId e, QueueState &state) = 0;

    public:
        AbstractGapFiller(const Graph &g, const Sequence &ss, EdgeId start_e, EdgeId end_e,
                                  const int start_p, const int end_p, const int path_max_length, 
                                  const std::map<VertexId, size_t> &vs, const std::map<VertexId, size_t> &vs_b)
                  :g_(g), ss_(ss)
                   , start_e_(start_e), end_e_(end_e)
                   , start_p_(start_p), end_p_(end_p)
                   , path_max_length_(path_max_length) 
                   , vs_(vs), vs_b_(vs_b){
            best_ed.resize(ss_.size());
            for (size_t i = 0; i < best_ed.size(); ++ i){
                best_ed[i] = path_max_length_;
            }
            AddNewEdge(GraphState(start_e_, start_p_, g_.length(start_e_)), QueueState(), 0);
            min_score = -1;
        }

        void CloseGap() {
            GraphState endgstate(end_e_, 0, end_p_);
            QueueState endqstate(endgstate, ss_.size() - 1);
            // GraphState startgstate(start_e_, start_p_, g_.length(start_e_));
            // QueueState startqstate(startgstate, 0);
            bool foundPath = false;
            int i = 0;
            while (q.size() > 0) {
                QueueState curState = q.begin()->second;
                int ed = dist[curState];
                if (q.size() > 1000000) {
                    INFO("queue size is too big")
                    if (visited.count(endqstate) > 0){
                        foundPath = true;
                        min_score = dist[endqstate];
                    }
                    break;
                }
                if (ed > path_max_length_) {
                    INFO("No path found")
                    break;
                }
                if (curState == endqstate && ed <= path_max_length_) {
                    min_score = ed;
                    INFO("+++++++++++++++++++++++++++++++++++++++++++++++=Final ed=" << ed);
                    foundPath = true;
                    break;
                }
                i ++;
                q.erase(q.begin());
                for (const auto &e: g_.OutgoingEdges(g_.EdgeEnd(curState.gs.e))) {
                    if (vs_.count(g_.EdgeEnd(curState.gs.e)) > 0) { //&& vs_.at(g_.EdgeEnd(curState.gs.e)) - (ss_.size() - curState.i) <= path_max_length_ - ed ) {
                        if (vs_.count(g_.EdgeEnd(e)) > 0) {
                            GraphState nextgs(e, 0, g_.length(e));
                            AddNewEdge(nextgs, curState, ed);
                        }
                        if (IsEnd(e, curState) && path_max_length_ - ed >= 0){
                            utils::perf_counter perf;
                            string sstr = ss_.Subseq(curState.i + 1, ss_.size() ).str();
                            string tmp = g_.EdgeNucls(e).str();
                            string edgestr = tmp.substr(0, end_p_);
                            //INFO("TIME.Substrs=" << perf.time());
                            perf.reset();
                            int score = NWDistance(sstr, edgestr, path_max_length_ - ed);
                            //INFO("TIME.NWDistance=" << perf.time());
                            perf.reset();
                            if (score != -1) {
                                if (ed +  score < path_max_length_) {
                                    INFO("Update max_path_len=" << ed + score);
                                }
                                path_max_length_ = min(path_max_length_, ed + score);
                                QueueState state(GraphState(e, 0, end_p_), ss_.size() - 1);
                                Update(state, curState, ed + score, ed + score);
                            }
                            //INFO("TIME.QueueUpdate=" << perf.time());
                        }
                    }
                }
            }

            if (foundPath) {
                QueueState state = endqstate;
                while (!state.isempty()) {
                    gapPath.push_back(state.gs.e);
                    state = prevStates[state];
                }
                std::reverse(gapPath.begin(), gapPath.end());
            }
            return;
        }

        std::vector<EdgeId> GapPath() {
            return gapPath;
        }

        int MinScore() {
            return min_score;
        }

    protected:
        std::set<std::pair<int, QueueState> > q;
        std::map<QueueState, int> visited; 
        std::map<QueueState, int> dist; 
        std::map<QueueState, QueueState> prevStates;

        std::vector<EdgeId> gapPath;

        std::vector<int> best_ed;

        const Graph &g_; 
        const Sequence &ss_;
        EdgeId start_e_;
        EdgeId end_e_;
        const int start_p_;
        const int end_p_;
        int path_max_length_;
        const std::map<VertexId, size_t> &vs_;
        const std::map<VertexId, size_t> &vs_b_;
        int min_score;
};

class GapFiller: public AbstractGapFiller {
    virtual bool IsEnd(EdgeId e, QueueState &state){
        if (e == end_e_) {
            return true;
        }
        return false;
    }

public:
    GapFiller(const Graph &g, const Sequence &ss, EdgeId start_e, EdgeId end_e,
                                  const int start_p, const int end_p, const int path_max_length, 
                                  const std::map<VertexId, size_t> &vs, const std::map<VertexId, size_t> &vs_b)
            :AbstractGapFiller(g, ss, start_e, end_e, start_p, end_p, path_max_length, vs, vs_b)
        {}

};

class EndsConstructor: public AbstractGapFiller {
    virtual bool IsEnd(EdgeId e, QueueState &state){
        int remaining = ss_.size() - state.i;
        if (g_.length(e) - remaining > 0){
            return true;
        }
        return false;
    }

public:
    EndsConstructor(const Graph &g, const Sequence &ss, EdgeId start_e, EdgeId end_e,
                                  const int start_p, const int end_p, const int path_max_length, 
                                  const std::map<VertexId, size_t> &vs, const std::map<VertexId, size_t> &vs_b)
            :AbstractGapFiller(g, ss, start_e, end_e, start_p, end_p, path_max_length, vs, vs_b)
        {}

};


}
