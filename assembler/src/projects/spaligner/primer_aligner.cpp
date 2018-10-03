#include "primer_aligner.hpp"

#include "sequence/sequence_tools.hpp"

namespace rna16S_mapping {

int EditDistance(const string &a, const string &b, int &start_pos, int &end_pos) {
    //INFO("Before ed")
    int a_len = (int) a.length();
    int b_len = (int) b.length();
    VERIFY(a_len > 0);
    VERIFY(b_len > 0);
    edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len, b.c_str(), b_len
                                     , edlib::edlibNewAlignConfig(a_len, edlib::EDLIB_MODE_HW, edlib::EDLIB_TASK_LOC,
                                             additionalEqualities, 36));
    int score = -1;
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        if (result.numLocations > 0) {
            score = result.editDistance;
            start_pos = result.startLocations[0];
            end_pos = result.endLocations[0];
        } else {
            WARN("EditDistance: something wrong with edlib result");
        }
    }
    edlib::edlibFreeAlignResult(result);
    return score;
}

void PrimerAligner::AlignPrimer(const io::SingleRead &read, int &v, std::vector<EdgeId> &e, std::vector<MappingRange> &range) {
    v = cfg_.graph_max_ed;
    int k = 10;
    for (auto it = gp_.g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId eid = *it;
        std::string edge_str = gp_.g.EdgeNucls(eid).str();
        int start_pos = 0;
        int end_pos = 1;
        int dist = EditDistance(read.GetSequenceString(), edge_str, start_pos, end_pos);
        if (dist == -1) continue;
        if (v >= dist && start_pos < gp_.g.length(eid)) {
            e.push_back(eid);
            range.push_back(MappingRange(Range(0, read.size()), Range(start_pos, end_pos + 1)) );
        }
    }

    std::string ans = " dist=" +  std::to_string(v) + " num=" +  std::to_string(e.size());
    // for (EdgeId eid: e) {
    //     ans += " " + std::to_string(eid.int_id()) + ",";
    // }
    // ans += "\n";
    // for (MappingRange r: range) {
    //     ans += " " + std::to_string(r.mapped_range.start_pos) + "-" + std::to_string(r.mapped_range.end_pos) + ",";
    // }
    INFO("Primer name=" << read.name() << " seq=" << read.GetSequenceString() << " " << ans);
}

void PrimerAligner::FormDistanceMatrix() {
    for (size_t i = 0; i < primers_.size(); ++ i) {
        for (auto e : primers_[i].e_) {
            VertexId end_v = gp_.g.EdgeEnd(e);
            if (dist_.count(end_v) == 0) {
                auto path_searcher = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(gp_.g, 1500);
                path_searcher.Run(end_v);
                auto reached_vertices = path_searcher.ProcessedVertices();
                dist_.insert(std::pair<VertexId, std::set<VertexId> > (end_v, std::set<VertexId>()));
                for (auto j_iter = reached_vertices.begin(); j_iter != reached_vertices.end(); ++j_iter) {
                    dist_[end_v].insert(*j_iter);
                }
            }
        }
    }
}

void PrimerAligner::PreparePrimers(std::vector<io::SingleRead> &wrappedprimers, int threads) {
    #pragma omp parallel num_threads(threads)
    #pragma omp for
    for (size_t i = 0 ; i < wrappedprimers.size(); ++i) {
        int v = 0;
        std::vector<EdgeId> e;
        std::vector<MappingRange> range;
        AlignPrimer(wrappedprimers[i], v, e, range);
        const std::string &name = wrappedprimers[i].name();
        const std::string &seq = wrappedprimers[i].GetSequenceString();
        INFO("i=" << i)
        #pragma omp critical
        {
            if (e.size() > 0) {
                primers_.push_back(ReadMapping(name, seq, v, e, range));
            }
        }
    }
    std::sort(primers_.begin(), primers_.end());
    FormDistanceMatrix();
}


} //namespace rna16S_mapping