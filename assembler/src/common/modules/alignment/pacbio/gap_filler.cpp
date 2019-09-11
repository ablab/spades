//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "modules/alignment/pacbio/gap_filler.hpp"

#include "assembly_graph/paths/path_utils.hpp"

namespace sensitive_aligner {

using namespace std;


string GapFiller::PathToString(std::vector<EdgeId> &path) const {
    Sequence merged_seq = MergeSequences(g_, path);  
    if (merged_seq.size() > 0) {
        return merged_seq.Subseq(0, merged_seq.size() - g_.k()).str(); 
    } else {
        return "";
    }
}

GapFillerResult GapFiller::BestScoredPathDijkstra(const string &s,
        const GraphPosition &start_pos,
        const GraphPosition &end_pos,
        int path_max_length, int score) const {
    GapClosingConfig gap_cfg = cfg_.gap_cfg;
    VertexId start_v = g_.EdgeEnd(start_pos.edgeid);
    VertexId end_v = g_.EdgeStart(end_pos.edgeid);
    auto path_searcher_b = omnigraph::DijkstraHelper<debruijn_graph::Graph>::
                                CreateBackwardBoundedDijkstra(g_, path_max_length);
    path_searcher_b.Run(end_v);
    auto path_searcher = omnigraph::DijkstraHelper<debruijn_graph::Graph>::
                                CreateBoundedDijkstra(g_, path_max_length);
    path_searcher.Run(start_v);
    const auto &reached_vertices_b = path_searcher_b.ProcessedVertices();
    const auto &reached_vertices = path_searcher.ProcessedVertices();

    unordered_map<VertexId, size_t> vertex_pathlen;
    for (auto v : reached_vertices_b) {
        if (reached_vertices.count(v) > 0) {
            vertex_pathlen[v] = path_searcher_b.GetDistance(v);
        }
    }
    int s_len = int(s.size());
    int ed_limit = score;
    if (score == numeric_limits<int>::max()) {
        ed_limit = min(max(gap_cfg.ed_lower_bound, s_len /
                            gap_cfg.max_ed_proportion), gap_cfg.ed_upper_bound);
    }
    DEBUG(" Dijkstra: String length " << s_len << "  "  << (size_t) s_len <<
          " max-len " << ed_limit << " vertex_num=" << vertex_pathlen.size());
    GapFillerResult dijkstra_res;
    if ((vertex_pathlen.size() == 0 ||
            ((size_t) s_len) * vertex_pathlen.size() > gap_cfg.max_gs_states) &&
            start_pos.edgeid.int_id() != end_pos.edgeid.int_id()) {
        DEBUG("Dijkstra won't run: Too big gap or too many paths " << s_len << " " << vertex_pathlen.size());
        if (vertex_pathlen.size() == 0) {
            dijkstra_res.return_code.not_connected = true;
        } else {
            dijkstra_res.return_code.long_gap = true;
            dijkstra_res.return_code.wide_gap = true;
        }
        return dijkstra_res;
    }
    DijkstraGapFiller gap_filler(g_, gap_cfg, s,
                                 start_pos.edgeid, end_pos.edgeid,
                                 (int) start_pos.position, (int) end_pos.position,
                                 ed_limit, vertex_pathlen);
    gap_filler.CloseGap();
    dijkstra_res.score = gap_filler.edit_distance();
    dijkstra_res.return_code = gap_filler.return_code();
    if (dijkstra_res.score == numeric_limits<int>::max()) {
        DEBUG("Dijkstra didn't find anything")
        return dijkstra_res;
    }
    dijkstra_res.full_intermediate_path = gap_filler.path();
    return dijkstra_res;
}

GapFillerResult GapFiller::BestScoredPathBruteForce(const string &seq_string,
        const GraphPosition &start_pos,
        const GraphPosition &end_pos,
        int path_min_length, int path_max_length) const {
    VertexId start_v = g_.EdgeEnd(start_pos.edgeid);
    VertexId end_v = g_.EdgeStart(end_pos.edgeid);
    TRACE(" Traversing tangled region. Start and end vertices resp: "
                    << g_.int_id(start_v) << " " << g_.int_id(end_v));
    omnigraph::PathStorageCallback<debruijn_graph::Graph> callback(g_);
    GapFillerResult bf_res;
    bf_res.return_code.no_path = true;
    int return_code = ProcessPaths(g_,
                                   path_min_length, path_max_length,
                                   start_v, end_v,
                                   callback);
    DEBUG("PathProcessor result: " << return_code << " limits " << path_min_length << " " << path_max_length);
    vector<vector<EdgeId> > paths = callback.paths();
    size_t best_path_ind = paths.size();
    int best_score = numeric_limits<int>::max();
    if (paths.size() == 0) {
        DEBUG("need to find best scored path between " << paths.size() << " , seq_len " << seq_string.length());
        DEBUG ("no paths");
        return bf_res;
    }
    if (seq_string.length() > cfg_.pb.max_contigs_gap_length) {
        DEBUG("need to find best scored path between " << paths.size() << " , seq_len " << seq_string.length());
        DEBUG("Gap is too large");
        return bf_res;
    }
    string s_add = g_.EdgeNucls(start_pos.edgeid).
                      Subseq(start_pos.position, g_.length(start_pos.edgeid)).str();
    string e_add = g_.EdgeNucls(end_pos.edgeid).Subseq(0, end_pos.position).str();
    bool additional_debug = (paths.size() > 1 && paths.size() < 10);
    for (size_t i = 0; i < paths.size(); i++) {
        DEBUG("path len " << paths[i].size());
        if (paths[i].size() == 0) {
            DEBUG ("Pathprocessor returns path with size = 0")
        }
        
        string cur_string = s_add + PathToString(paths[i]) + e_add;
        TRACE("cur_string: " << cur_string << "\n seq_string " << seq_string);
        int cur_score = StringDistance(cur_string, seq_string);
        //DEBUG only
        if (additional_debug) {
            TRACE("candidate path number " << i << " , len " << cur_string.length());
            TRACE("graph candidate: " << cur_string);
            TRACE("in pacbio read: " << seq_string);
            for (auto j_iter = paths[i].begin(); j_iter != paths[i].end();
                    ++j_iter) {
                DEBUG(g_.int_id(*j_iter));
            }
            DEBUG("score: " << cur_score);
        }
        if (cur_score < best_score) {
            best_score = cur_score;
            best_path_ind = i;
        }
    }
    TRACE(best_score);
    bf_res.score = best_score;
    if (best_score == numeric_limits<int>::max()) {
        if (paths.size() < 10) {
            for (size_t i = 0; i < paths.size(); i++) {
                DEBUG ("failed with strings " << seq_string <<
                        " " << s_add + PathToString(paths[i]) + e_add);
            }
        }
        DEBUG (paths.size() << " paths available");
        return bf_res;
    }
    if (additional_debug) {
        TRACE("best score found! Path " << best_path_ind << " score " << best_score);
    }
    bf_res.full_intermediate_path = paths[best_path_ind];
    bf_res.full_intermediate_path.insert(bf_res.full_intermediate_path.begin(), start_pos.edgeid);
    bf_res.full_intermediate_path.push_back(end_pos.edgeid);
    if (return_code == 0) {
        bf_res.return_code.status = 0;
    }
    return bf_res;
}

GapFillerResult GapFiller::Run(const string &s,
                               const GraphPosition &start_pos,
                               const GraphPosition &end_pos,
                               int path_min_length, int path_max_length) const {
    utils::perf_counter pc;
    GapClosingConfig gap_cfg = cfg_.gap_cfg;
    auto bf_res = BestScoredPathBruteForce(s, start_pos, end_pos, path_min_length, path_max_length);
    double bf_time = pc.time();
    pc.reset();
    if (gap_cfg.run_dijkstra && bf_res.return_code.status != 0) {
        auto dijkstra_res = BestScoredPathDijkstra(s, start_pos, end_pos, path_max_length, bf_res.score);
        DEBUG("BruteForce run: return_code=" << bf_res.return_code.status
              << " score=" << bf_res.score << " time_bf=" << bf_time
              << " Dijkstra run: return_code=" << dijkstra_res.return_code.status
              << " score=" <<  dijkstra_res.score << " time_d=" << pc.time() << " len=" << s.size() << "\n")
        if (dijkstra_res.return_code.status == 0) {
            return dijkstra_res;
        }
    }
    if (bf_res.score != numeric_limits<int>::max() && !gap_cfg.find_shortest_path) {
        bf_res.return_code.status = 0;
    }
    return bf_res;

}

//////////////////////////////////////////////////////////////////

GraphPosition GapFiller::ConjugatePosition(const GraphPosition &start_pos) const {
    GraphPosition cstart_pos;
    cstart_pos.edgeid = g_.conjugate(start_pos.edgeid);
    cstart_pos.position = min((int) g_.length(start_pos.edgeid),
                             (int) g_.length(start_pos.edgeid) + (int) g_.k() - (int) start_pos.position);
    DEBUG("Backward e=" << cstart_pos.edgeid.int_id() << " sp=" << cstart_pos.position);
    return cstart_pos;
}

void GapFiller::UpdatePath(vector<debruijn_graph::EdgeId> &path,
                           vector<EdgeId> &ans,
                           MappingPoint p, PathRange &range, bool forward, GraphPosition &old_start_pos) const {
    if (forward) {
        size_t end_pos = p.edge_pos;
        size_t end_pos_seq = p.seq_pos;
        while (end_pos < g_.k() && ans.size() > 0) {
            ans.pop_back();
            end_pos += g_.length(ans[ans.size() - 1]);
        }
        for (int i = 1; i < (int) ans.size(); ++i) {
            path.push_back(ans[i]);
        }
        range.path_end.seq_pos = end_pos_seq;
        range.path_end.edge_pos = end_pos;
    } else {
        vector<debruijn_graph::EdgeId> cur_sorted;
        size_t end_pos = p.edge_pos;
        int start = (int) g_.length(ans[ans.size() - 1]) + (int) g_.k() - (int) end_pos;
        int cur_ind = (int) ans.size() - 1;
        while (cur_ind >= 0 && start - (int) g_.length(ans[cur_ind]) > 0) {
            start -= (int) g_.length(ans[cur_ind]);
            cur_ind --;
        }
        if (cur_ind > 0) {
            cur_sorted.push_back(g_.conjugate(ans[cur_ind]) );
        }
        for (int i = cur_ind - 1; i > 0; --i) {
            cur_sorted.push_back(g_.conjugate(ans[i]));
        }
        for (size_t i = 0; i < path.size(); ++i) {
            cur_sorted.push_back(path[i]);
        }
        if (path.size() == cur_sorted.size() && start > (int) old_start_pos.position) {
            return;
        }
        path = cur_sorted;
        range.path_start.seq_pos = p.seq_pos;
        range.path_start.edge_pos = start;
    }
}

GapFillerResult GapFiller::Run(Sequence &s,
                               GraphPosition &start_pos,
                               bool forward,
                               vector<debruijn_graph::EdgeId> &path,
                               PathRange &range) const {
    VERIFY(path.size() > 0);
    EndsClosingConfig ends_cfg = cfg_.ends_cfg;
    GraphPosition old_start_pos = start_pos;
    if (!forward) {
        s = !s;
        start_pos = ConjugatePosition(start_pos);
    }
    GapFillerResult res;
    size_t s_len = int(s.size());
    int score =  min(min(
                        max(ends_cfg.ed_lower_bound, (int) s_len / ends_cfg.max_ed_proportion),
                        ends_cfg.ed_upper_bound),
                    (int) s_len);
    if ((int) s_len >  ends_cfg.max_restorable_length && s_len >  g_.length(start_pos.edgeid) - start_pos.position + g_.k()) {
        DEBUG("EdgeDijkstra: sequence is too long " << s_len)
        res.return_code.long_gap = true;
        return res;
    }
    if (s_len < 1) {
        DEBUG("EdgeDijkstra: sequence is too small " << s_len)
        return res;
    }
    utils::perf_counter pc;
    DijkstraEndsReconstructor algo(g_, ends_cfg, s.str(), start_pos.edgeid, (int) start_pos.position, score);
    algo.CloseGap();
    score = algo.edit_distance();
    res.return_code = algo.return_code();
    if (score == numeric_limits<int>::max()) {
        DEBUG("EdgeDijkstra didn't find anything edge=" << start_pos.edgeid.int_id()
              << " s_start=" << start_pos.position << " seq_len=" << s.size())
        return res;
    }
    vector<EdgeId> ans = algo.path();
    MappingPoint p(forward ? algo.seq_end_position() + range.path_end.seq_pos : range.path_start.seq_pos - algo.seq_end_position(), algo.path_end_position());
    UpdatePath(path, ans, p, range, forward, old_start_pos);
    return res;
}

} // namespace sensitive_aligner