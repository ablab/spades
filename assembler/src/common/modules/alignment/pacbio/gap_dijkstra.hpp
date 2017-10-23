#include "assembly_graph/index/edge_multi_index.hpp"
#include "modules/alignment/edge_index_refiller.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"

namespace pacbio {

struct MetaVertex {
        EdgeId e;
        int j;
        int i;
        
        MetaVertex()
            : e(EdgeId()), j(-1), i(-1)
        {}
        
        MetaVertex(const EdgeId &_e
                   , const int &_j
                   , const int &_i)
            : e(_e), j(_j), i(_i)
        {}
        
        MetaVertex(const MetaVertex &mv)
            : e(mv.e), j(mv.j), i(mv.i)
        {}
        
        bool operator == (const MetaVertex &mv) const {
            return (this->e == mv.e  && this->i == mv.i && this->j == mv.j);
        }
        
        bool operator != (const MetaVertex &mv) const {
            return (this->e != mv.e || this->i != mv.i || this->j != mv.j);
        }
        
        bool operator < (const MetaVertex &mv) const {
            return (this->e < mv.e || ( this->e == mv.e && this->i < mv.i) || ( this->e == mv.e && this->i ==  mv.i && this->j < mv.j));
        }
    };
    

class Dijkstra {
        typedef typename boost::bimap< size_t, MetaVertex >::value_type state;
    public:
        Dijkstra(const Graph &g, const Sequence &ss, EdgeId start_e, EdgeId end_e,
                                  const int start_p, const int end_p, const int path_max_length, std::map<VertexId, size_t> &vertex_pathlen)
        : score_(-1)
        , g_(g)
        , ss_(ss)
        , start_e_(start_e)
        , end_e_(end_e)
        , start_p_(start_p)
        , end_p_(end_p)
        , path_max_length_(path_max_length)
	, vertex_pathlen_(vertex_pathlen) {
            num = 0;
            MetaVertex mv(start_e_, start_p_ - 1, -1);
            size_t newId = num ++;
            bids.insert( state(newId, mv ) );
            //INFO("Dijkstra initial e=" << start_e_.int_id() << " pos=" << start_p_-1)
            q.insert(make_pair (0, newId));
            dist[newId] = 0;
            prev[newId] = newId; 

            MetaVertex mv_b(end_e_, end_p_ + 1, (int) ss_.size());
            newId = num ++;
            bids.insert( state(newId, mv_b ) );
            //INFO("Dijkstra initial e=" << start_e_.int_id() << " pos=" << start_p_-1)
            q_b.insert(make_pair (0, newId));
            dist_b[newId] = 0;
            prev_b[newId] = newId; 

        }
        
        void AddMetaVertex(size_t &fromId, const MetaVertex &to, int score){
            auto it = bids.right.find(to);
            if (it == bids.right.end() || dist[it->second] > dist[fromId] + score){
                size_t toId = 0;
                if (it != bids.right.end()){
                    q.erase(make_pair(dist[it->second], it->second));
                    toId = it->second;
                }else{
                    toId = num ++;
                    bids.insert( state(toId, to ) );
                }
                dist[toId] = dist[fromId] + score;
                prev[toId] = fromId;
                q.insert(make_pair(dist[toId], toId));
                //INFO("Dijkstra queue: dist=" << dist[toId] << " edge=" << to.e.int_id() << " pos=" << to.j << " pos_s=" << to.i << " cur_id="<< toId);
            }
        }

        void AddMetaVertexBackward(size_t &fromId, const MetaVertex &to, int score){
            auto it = bids.right.find(to);
            if (it == bids.right.end() || dist_b[it->second] > dist_b[fromId] + score){
                size_t toId = 0;
                if (it != bids.right.end()){
                    q_b.erase(make_pair(dist_b[it->second], it->second));
                    toId = it->second;
                }else{
                    toId = num ++;
                    bids.insert( state(toId, to ) );
                }
                dist_b[toId] = dist_b[fromId] + score;
                prev_b[toId] = fromId;
                q_b.insert(make_pair(dist_b[toId], toId));
                //INFO("Dijkstra queue: dist=" << dist[toId] << " edge=" << to.e.int_id() << " pos=" << to.j << " pos_s=" << to.i << " cur_id="<< toId);
            }
        }
        
        void AddNextToForwardQueue(size_t mv_id, MetaVertex mv, int r_len) {
            int mu = 1; //subst cost
            int sigma = 1; // indel cost
            int K = int(g_.k());
            int s_len = int(ss_.size());
            if (mv.j < r_len - K - 1){
                MetaVertex to(mv.e, mv.j + 1, mv.i);
                AddMetaVertex(mv_id, to, sigma);
                if (mv.i < s_len - 1){
                    to = MetaVertex(mv.e, mv.j, mv.i + 1);
                    AddMetaVertex(mv_id, to, sigma);

                    to = MetaVertex(mv.e, mv.j + 1, mv.i + 1);
                    int score = 0;
                    if (g_.EdgeNucls(mv.e)[mv.j + 1] != ss_[mv.i + 1]){
                        score = mu;
                    } 
                    AddMetaVertex(mv_id, to, score);
                }
            } else{
                if (mv.i < s_len - 1){
                    MetaVertex to(mv.e, mv.j, mv.i + 1);
                    AddMetaVertex(mv_id, to, sigma);
                }
                int pos = r_len - K - mv.j - 1;
                VertexId endv = g_.EdgeEnd(mv.e);
		if (vertex_pathlen_.count(endv) > 0 && vertex_pathlen_[endv] < path_max_length_ - dist[mv_id] + (s_len - mv.i) + 1) {
                        vertices.insert(endv);
                        for (const auto &e: g_.OutgoingEdges(endv)) {
                                MetaVertex to(e, pos, mv.i);
                                AddMetaVertex(mv_id, to, sigma);
                                if (mv.i < s_len - 1){
                                        to = MetaVertex(e, pos, mv.i + 1);
                                        int score = 0;
                                        if (g_.EdgeNucls(e)[pos] != ss_[mv.i + 1]){
                                                score = mu;
                                        }
                                        AddMetaVertex(mv_id, to, score);
                                }
                        }
                }
            }
        }

        void AddNextToBackwardQueue(size_t mv_id, MetaVertex mv) {
            int mu = 1; //subst cost
            int sigma = 1; // indel cost
            int K = int(g_.k());
            if (mv.j > K){
                MetaVertex to(mv.e, mv.j - 1, mv.i);
                AddMetaVertexBackward(mv_id, to, sigma);
                if (mv.i > 0){
                    to = MetaVertex(mv.e, mv.j, mv.i - 1);
                    AddMetaVertexBackward(mv_id, to, sigma);

                    to = MetaVertex(mv.e, mv.j - 1, mv.i - 1);
                    int score = 0;
                    if (g_.EdgeNucls(mv.e)[mv.j - 1] != ss_[mv.i - 1]){
                        score = mu;
                    } 
                    AddMetaVertexBackward(mv_id, to, score);
                }
            } else{
                if (mv.i > 0){
                    MetaVertex to(mv.e, mv.j, mv.i - 1);
                    AddMetaVertexBackward(mv_id, to, sigma);
                }
                VertexId startv = g_.EdgeStart(mv.e);
                vertices.insert(startv);
                for (const auto &e: g_.IncomingEdges(startv)) {
                    int pos = int(g_.EdgeNucls(e).size()) - K + mv.j - 1;
                    //INFO("New edge " << g_.int_id(e) << " pos=" << pos  << " mv_i=" << mv.i ) 
                    //INFO("Nucls edge " << g_.int_id(e) << " " << g_.EdgeNucls(e) << "\n" << " Seq nucls " << ss_.str())
                    MetaVertex to(e, pos, mv.i);
                    AddMetaVertexBackward(mv_id, to, sigma);
                    if (mv.i > 0){
                        to = MetaVertex(e, pos, mv.i - 1);
                        int score = 0;
                        if (g_.EdgeNucls(e)[pos] != ss_[mv.i - 1]){
                            score = mu;
                        } 
                        AddMetaVertexBackward(mv_id, to, score);
                    }
                } 
            }
        }

        std::vector<EdgeId> RestoreForwardPath(MetaVertex final_vertex) {
            vector<EdgeId> ans;
            int tid = omp_get_thread_num();
            DEBUG("Dijkstra forwardrestore " << " mv_i=" << final_vertex.i  << " mv_j=" << final_vertex.j << " edgeId=" << g_.int_id(final_vertex.e))
            if (final_vertex != MetaVertex()){
                MetaVertex cur_mv = final_vertex;
                std::string final_str = "";
                size_t cur_id = bids.right.at(cur_mv);
                DEBUG(" Dijkstra: Final score " << dist[cur_id] ) 
                score_ = dist[cur_id];
                DEBUG("Dijkstra forwardrestore2 " << " mv_i=" << final_vertex.i  << " mv_j=" << final_vertex.j << " edgeId=" << g_.int_id(final_vertex.e))
                while (cur_id != prev[cur_id]) {
                    cur_mv = bids.left.at(cur_id);
                    MetaVertex prev_mv = bids.left.at(prev[cur_id]);
                    if (cur_mv.e != prev_mv.e){
                        ans.push_back(cur_mv.e);
                    }
                    cur_id = prev[cur_id];
                }
                cur_mv = bids.left.at(cur_id);
                ans.push_back(cur_mv.e);
                DEBUG("Dijkstra forwardrestore3 " << " mv_i=" << final_vertex.i  << " mv_j=" << final_vertex.j << " edgeId=" << g_.int_id(final_vertex.e))
                std::reverse(ans.begin(), ans.end());
                DEBUG("ThreadNum=" << tid << ". Ended Path construction size=" << ans.size());
                return ans;
            } else {
                DEBUG(" Dijkstra: final vertex is empty");
                return vector<EdgeId>(0);
            }

        }

        std::vector<EdgeId> RestoreBackwardPath(MetaVertex final_vertex_b) {
            vector<EdgeId> ans;
            int tid = omp_get_thread_num();
            DEBUG("Dijkstra  backwardrestore " << " mv_i=" << final_vertex_b.i  << " mv_j=" << final_vertex_b.j << " edgeId=" << g_.int_id(final_vertex_b.e))
            if (final_vertex_b != MetaVertex()){
                MetaVertex cur_mv = final_vertex_b;
                std::string final_str = "";
                size_t cur_id = bids.right.at(cur_mv);
                DEBUG(" Dijkstra: Final score " << dist_b[cur_id] ) 
                score_ = dist_b[cur_id];
                DEBUG("Dijkstra  backwardrestore2 " << " mv_i=" << final_vertex_b.i  << " mv_j=" << final_vertex_b.j << " edgeId=" << g_.int_id(final_vertex_b.e)
                                                    << " cur_id =" << cur_id << " prev_b="<< prev_b[cur_id])
                while (cur_id != prev_b[cur_id]) {
                    cur_mv = bids.left.at(cur_id);
                    MetaVertex prev_mv = bids.left.at(prev_b[cur_id]);
                    if (cur_mv.e != prev_mv.e){
                        ans.push_back(cur_mv.e);
                    }
                    cur_id = prev_b[cur_id];
                    //INFO("Dijkstra  backwardrestore-c " << " mv_i=" << cur_mv.i  << " mv_j=" << cur_mv.j << " edgeId=" << g_.int_id(cur_mv.e))
                }
                cur_mv = bids.left.at(cur_id);
                ans.push_back(cur_mv.e);
                DEBUG("Dijkstra  backwardrestore3 " << " mv_i=" << final_vertex_b.i  << " mv_j=" << final_vertex_b.j << " edgeId=" << g_.int_id(final_vertex_b.e))
                //std::reverse(ans.begin(), ans.end());
                DEBUG("backward ThreadNum=" << tid << ". Ended Path construction size=" << ans.size());
                return ans;
            } else {
                DEBUG(" Dijkstra: final vertex is empty");
                return vector<EdgeId>(0);
            }
        }
        
        std::vector<EdgeId> Run() {
            MetaVertex final_vertex = MetaVertex();
            int s_len = int(ss_.size());
            int tid = omp_get_thread_num();
            DEBUG("ThreadNum=" << tid << ". Run Dijkstra. Seq_len=" << s_len);
            while (!q.empty()){
                if (q.size() >= 1000) {
                    DEBUG("ThreadNum=" << tid << ". Queue size: " << q.size() << " " << vertices.size()  );
                    DEBUG("ThreadNum=" << tid << ". Stop!");
                    break;
                }
                size_t mv_id = q.begin()->second;
                MetaVertex mv = bids.left.at(mv_id);
                bids_best.insert(state(mv_id, mv));
                //INFO("Dijsktra queue cur_id="<< mv_id << " dist=" << dist[mv_id] << " mv_i=" << mv.i  << " mv_j=" << mv.j << " edgeId=" << g_.int_id(mv.e) );
                q.erase(q.begin());
                int r_len = int(g_.EdgeNucls(mv.e).size());
                if (dist[mv_id] > path_max_length_ || (mv.e == end_e_ && mv.i == s_len - 1 && mv.j == end_p_)){
                    if (mv.e == end_e_ && mv.i == s_len - 1 && mv.j == end_p_){
                        final_vertex = mv;
                    }
                    break;
                }

                AddNextToForwardQueue(mv_id, mv, r_len);
            }
            MetaVertex final_vertex_b = MetaVertex();
/*            if (final_vertex == MetaVertex()) {
                while (!q_b.empty()){
                    if (q_b.size() >= 2500 || vertices.size() > 20) {
                        DEBUG("ThreadNum=" << tid << ". Queue size: " << q_b.size() << " " << vertices.size() );
                        DEBUG("ThreadNum=" << tid << ". Stop!");
                        break;
                    }
                    size_t mv_id = q_b.begin()->second;
                    MetaVertex mv = bids.left.at(mv_id);
                    if (bids_best.left.count(mv_id) > 0){
                        bids_best_b.insert(state(mv_id, mv));
                    }
                    //INFO("Dijsktra queue cur_id="<< mv_id << " dist=" << dist_b[mv_id] << " mv_i=" << mv.i  << " mv_j=" << mv.j << " edgeId=" << g_.int_id(mv.e) );
                    q_b.erase(q_b.begin());
                    
                    if (dist_b[mv_id] > path_max_length_ || (mv.e == start_e_ && mv.i == 0 && mv.j == start_p_)){
                        if (mv.e == start_e_ && mv.i == 0 && mv.j == start_p_){
                            final_vertex_b = mv;
                        }
                        
                        break;
                    }

                    AddNextToBackwardQueue(mv_id, mv);
                }
            }
*/            MetaVertex final_vertex_fb = MetaVertex();
            int best_score = path_max_length_;
/*            if (final_vertex_b == MetaVertex() && final_vertex == MetaVertex()){
                DEBUG("Dijkstra FB working ThreadNum = " << tid << " num=" << bids.size())
                typename boost::bimap< size_t, MetaVertex >::const_iterator iter;
                for (iter = bids.begin(); iter != bids.end(); ++iter ) {
                    int i = (int) iter->left;
                    if (dist.count(i) > 0 && dist_b.count(i) > 0){
                        int cur_score = dist[i] + dist_b[i];
                        if (cur_score < best_score){
                            final_vertex_fb = iter->right;
                            best_score = cur_score;
                            DEBUG("Dijkstra FB working ThreadNum = " << tid << " score=" << cur_score << " i=" << i << " mv_i=" 
                                                                << final_vertex_fb.i  << " mv_j=" << final_vertex_fb.j << " edgeId=" << g_.int_id(final_vertex_fb.e));
                        }
                    }
                }
            }
*/            DEBUG("Dijkstra ThreadNum = " << tid << " num=" << bids.size())
            if (final_vertex != MetaVertex()) {
                DEBUG("For worked! ")
                return RestoreForwardPath(final_vertex);
            }
            if (final_vertex_b != MetaVertex()) {
                DEBUG("Back worked! ")
                return RestoreBackwardPath(final_vertex_b);
            }
            if (final_vertex_fb != MetaVertex()) {
                INFO("FB worked! " << best_score << " ThreadNum=" << tid << " s=" << g_.int_id(start_e_) << " e="<< g_.int_id(end_e_) << " mv_i=" 
                                                            << final_vertex_fb.i  << " mv_j=" << final_vertex_fb.j << " edgeId=" << g_.int_id(final_vertex_fb.e) << " score=" <<best_score)
                std::vector<EdgeId> fpath = RestoreForwardPath(final_vertex_fb);
                std::vector<EdgeId> bpath = RestoreBackwardPath(final_vertex_fb);
                std::vector<EdgeId> fbpath;
                fbpath.reserve(fpath.size() + bpath.size() - 1);
                fbpath.insert(fbpath.end(), fpath.begin(), fpath.end() - 1);
                fbpath.insert(fbpath.end(), bpath.begin(), bpath.end());
                INFO("FB worked! " << best_score << " ThreadNum=" << tid << " fb_sz = " << fbpath.size() <<  " fpath=" << fpath.size() << " bpath=" << bpath.size());
                return fbpath;
            } else {
               DEBUG(" Dijkstra: final vertex is empty");
               return vector<EdgeId>(0); 
            }
            DEBUG(" Dijkstra: final vertex is empty");
            return vector<EdgeId>(0); 
        }

        int score_;

    private:
        const Graph &g_;
        const Sequence &ss_;
        EdgeId start_e_;
        EdgeId end_e_;
        int start_p_;
        int end_p_;
        int path_max_length_;
        std::map<VertexId, size_t> &vertex_pathlen_;
        boost::bimap< size_t, MetaVertex > bids;
        boost::bimap< size_t, MetaVertex > bids_best;
        boost::bimap< size_t, MetaVertex > bids_best_b;
        size_t num;
        std::set< std::pair<int, size_t> > q;
        std::map<size_t, int> dist;
        std::map<size_t, size_t> prev;
        set<VertexId> vertices;       

        std::map<size_t, int> dist_b;
        std::map<size_t, size_t> prev_b;
        std::set< std::pair<int, size_t> > q_b;
        set<VertexId> vertices_b;
    };
}