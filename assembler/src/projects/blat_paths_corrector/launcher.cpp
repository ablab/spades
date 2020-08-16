//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "common.hpp"

#include "utils/logger/logger.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "common/assembly_graph/core/basic_graph_stats.hpp"
#include "common/assembly_graph/core/graph.hpp"
#include "common/assembly_graph/paths/path_processor.hpp"
#include "common/assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include "common/modules/path_extend/pe_utils.hpp"

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <boost/optional.hpp>

using namespace debruijn_graph;
using namespace path_extend;
using namespace std;

namespace {

struct PathWithBorderEdgesIndexies {
    GappedPath path;
    size_t first_edge_index;
    size_t last_edge_index;
};

class SequenceCorrector {
    Graph const & graph;
    PathThreadingParams params;
    GraphCoverageMap const & scaffold_coverage_map;
    PathWithEdgePostions const & path_info;
    std::string const & seq;
public:
    SequenceCorrector(Graph const & graph, PathThreadingParams params, GraphCoverageMap const & scaffold_coverage_map, PathWithEdgePostions const & path_info, std::string const & seq)
        : graph(graph)
        , params(params)
        , scaffold_coverage_map(scaffold_coverage_map)
        , path_info(path_info)
        , seq(seq)
    {}

    pair<std::string, std::vector<PathWithBorderEdgesIndexies>> GetBestSequence() const;
private:
    pair<GappedPath, size_t> GetNextPath(size_t start_pos) const;
    
    boost::optional<pair<GappedPath, size_t>> FindFiller(size_t start_pos) const;

    GappedPath GetBestMachedPath(EdgeId start_edge, EdgeId end_edge, size_t distance) const;

    GappedPath ConnectWithScaffolds(EdgeId start, EdgeId end) const;

    size_t EdgeLen(EdgeId edge) const { return graph.length(edge); }

    EdgeId GetEdge(size_t pos) const { 
        VERIFY_MSG(path_info.edge_set[pos].size() == 1, "alternatives are not supported");
        return path_info.edge_set[pos].front(); 
    }

    int GetPos(size_t pos) const { return path_info.positions[pos]; }

    size_t Size() const noexcept { return path_info.edge_set.size(); }
    
    size_t GetDistance(size_t start_pos, size_t end_pos) const { 
        size_t start_after_edge = GetPos(start_pos) + EdgeLen(GetEdge(start_pos));
        if (start_after_edge <= GetPos(end_pos))
            return GetPos(end_pos) - start_after_edge; 
        DEBUG("Negative distance between " << graph.int_id(GetEdge(start_pos)) << " and " << graph.int_id(GetEdge(end_pos)));
        return -1ull;
        // return 0;
    }

    std::string BackMapping(vector<PathWithBorderEdgesIndexies> const & paths) const;
};

std::string SequenceCorrector::BackMapping(vector<PathWithBorderEdgesIndexies> const & paths) const {
    stringstream ss;
    ScaffoldSequenceMaker seq_maker(graph);
    auto current_seq_pos = std::numeric_limits<long long>::min();
    for (size_t i = 0; i < paths.size(); ++i) {
        auto const & path = paths[i];
        size_t should_be_dropped = 0;
        BidirectionalPath path_seq(graph, path.path);
        auto new_seq = seq_maker.MakeSequence(path_seq);
        VERIFY(new_seq.size() == path_seq.Length() + graph.k());

        auto right_shift = graph.length(path_seq.Back()) + graph.k();
        auto start_pos = GetPos(path.first_edge_index);
        auto end_pos = GetPos(path.last_edge_index) + right_shift;
        if (current_seq_pos < start_pos) {
            if (0 <= current_seq_pos && current_seq_pos < seq.size())
                ss << seq.substr(current_seq_pos, start_pos - current_seq_pos);
        } else {
            should_be_dropped = current_seq_pos - start_pos;
            WARN("The new path prefix with len=" << should_be_dropped << " of " << new_seq.size() << " would be dropped");
        }
        if (new_seq.size() != end_pos - start_pos) {
            DEBUG("Replasing contig part with len=" << end_pos - start_pos << " using path with len=" << new_seq.size());
        }
        if (should_be_dropped >= new_seq.size()) {
            WARN("Oh no! Full path would be dropped");
            continue;
        }

        ss << new_seq.substr(should_be_dropped);
        current_seq_pos = end_pos;
    }
    current_seq_pos = std::max(0ll, current_seq_pos);
    if (current_seq_pos < seq.size())
        ss << seq.substr(current_seq_pos);
    return ss.str();
}

pair<std::string, std::vector<PathWithBorderEdgesIndexies>> SequenceCorrector::GetBestSequence() const {
    std::vector<PathWithBorderEdgesIndexies> paths;
    size_t start_pos = 0;
    while (start_pos < Size()) {
        auto path = GetNextPath(start_pos);
        if (path.first.Size() > 1) {
            paths.push_back({std::move(path.first), start_pos, path.second});
        }
        start_pos = path.second + 1;
    }
    return {BackMapping(paths), std::move(paths)};
}

pair<GappedPath, size_t> SequenceCorrector::GetNextPath(size_t start_pos) const {
    GappedPath current_path;
    current_path.PushBack(GetEdge(start_pos));
    while (start_pos + 1 < Size()) {
        auto path = FindFiller(start_pos);
        if (!path.is_initialized())
            break;
        current_path.PushBack(std::move(path->first));
        start_pos = path->second;
        if (current_path.Back() != GetEdge(start_pos))
            current_path.PushBack(GetEdge(start_pos));
    }
    return {std::move(current_path), start_pos};
}

boost::optional<pair<GappedPath, size_t>> SequenceCorrector::FindFiller(size_t start_pos) const {
    size_t end_pos = start_pos + 1;
    auto start_edge = GetEdge(start_pos);
    for (; end_pos - start_pos <= params.max_steps_forward && end_pos < Size(); ++end_pos) {
        GappedPath current_path;
        auto end_edge  = GetEdge(end_pos);

        if (graph.EdgeStart(end_edge) == graph.EdgeEnd(start_edge)) {
            return {{std::move(current_path), end_pos}};
        }
 
        size_t distance = GetDistance(start_pos, end_pos);
        if (distance == -1ull)
            continue;
        current_path = GetBestMachedPath(start_edge, end_edge, distance);
        if (params.use_scaffolds && current_path.Empty())
            current_path = ConnectWithScaffolds(start_edge, end_edge);
 
        if (!current_path.Empty()) {
            return {{std::move(current_path), end_pos}};
        }
    }
    return {};
}

GappedPath SequenceCorrector::GetBestMachedPath(EdgeId start_edge, EdgeId end_edge, size_t distance) const {
    VertexId target_vertex = graph.EdgeStart(end_edge);
    VertexId start_vertex = graph.EdgeEnd(start_edge);
    omnigraph::PathStorageCallback<Graph> path_storage(graph);
    auto min_len = static_cast<size_t>(std::floor((double)distance*(1-params.good_distance_coeff)));
    auto max_len = params.max_distance;
    cout << "<1>-----------------------\n";
    omnigraph::ProcessPaths(graph, min_len, max_len, start_vertex, target_vertex, path_storage);
    const auto& detected_paths = path_storage.paths();
    cout << "<cut here>----------------\n";

    auto max_len2 = static_cast<size_t>(std::ceil ((double)distance*(1+params.good_distance_coeff)));
    omnigraph::PathStorageCallback<Graph> path_storage2(graph);
    cout << "<2>-----------------------\n";
    omnigraph::ProcessPaths(graph, min_len, max_len2, start_vertex, target_vertex, path_storage2);
    const auto& detected_paths2 = path_storage2.paths();
    cout << "<cut here>----------------\n";
    
    GappedPath answer;
    vector<pair<size_t, size_t>> scores;
    vector<pair<size_t, size_t>> scores2;
    if (distance > 0) {
        for (size_t i = 0; i < detected_paths.size(); ++i) {
            auto path_len = CumulativeLength(graph, detected_paths[i]);
            size_t difference = (size_t) abs((long long)(path_len) - (long long)(distance));
            if (math::le(double(difference), double(distance) * params.good_distance_coeff))
                scores.emplace_back(difference, i);
        }
        for (size_t i = 0; i < detected_paths2.size(); ++i) {
            auto path_len = CumulativeLength(graph, detected_paths2[i]);
            size_t difference = (size_t) abs((long long)(path_len) - (long long)(distance));
            if (math::le(double(difference), double(distance) * params.good_distance_coeff))
                scores2.emplace_back(difference, i);
        }
        if (scores.size() != scores2.size()) {
            cout << "wtf, scores1 = " << scores.size() << ", scores2 = " << scores2.size() << "\n";
            cout << "scores1:\n";
            for (auto const & x : scores) {
                cout << "diff = " << x.first << ", len = " << CumulativeLength(graph, detected_paths[x.second]) << "\n";
            }
            cout << "scores2:\n";
            for (auto const & x : scores2) {
                cout << "diff = " << x.first << ", len = " << CumulativeLength(graph, detected_paths2[x.second]) << '\n';
            }
            cout << "max_len = " << max_len << "\n"
            "max_len2 = " << max_len2 << "\n"
            "distance = " << distance << "\n"
            "double(distance) * params.good_distance_coeff = " << double(distance) * params.good_distance_coeff << "\n";
            cout << "paths1:\n";
            for (auto const & x : scores) {
                BidirectionalPath p(graph, detected_paths[x.second]);
                p.PrintINFO();
                cout << "====================\n";
            }
            cout << "paths2:\n";
            for (auto const & x : scores2) {
                BidirectionalPath p(graph, detected_paths2[x.second]);
                p.PrintINFO();
                cout << "====================\n";
            }
            exit(2);
        }
        auto by_dist = [](const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b) { return a.first < b.first; };
        std::sort(scores.begin(), scores.end(), by_dist);
        std::sort(scores2.begin(), scores2.end(), by_dist);
        for (size_t i = 0; i < scores.size(); ++i) {
            VERIFY_MSG(scores[i].first == scores2[i].first, "d1 = " << scores[i].first << ", d2 = " << scores2[i].first);
        }
        if (scores.size() > 1) {
            if (math::le(double(scores[0].first), double(scores[1].first) * params.best_of_good_coeff))
                scores.erase(scores.begin() + 1, scores.end());
        }
    }

    if (scores.size() == 1) {
        size_t index = scores.front().second;
        answer.PushBack(detected_paths[index]);
    }
    return answer;
}

GappedPath SequenceCorrector::ConnectWithScaffolds(EdgeId start, EdgeId end) const {
    auto cover_start = scaffold_coverage_map.GetCoveringPaths(start);
    auto cover_end = scaffold_coverage_map.GetCoveringPaths(end);
    BidirectionalPathSet intersection;
    for (auto p : cover_start) {
        if (cover_end.count(p) > 0)
            intersection.insert(p);
    }
    GappedPath answer;
    if (intersection.size() == 1) {
        auto connecting_scaffold = *intersection.begin();
        auto start_pos = connecting_scaffold->FindAll(start);
        auto end_pos = connecting_scaffold->FindAll(end);
        if (start_pos.size() == 1 && end_pos.size() == 1 && start_pos.front() < end_pos.front()) {
            DEBUG("Found unique connecting scaffold " << start_pos.front() + 1 << " " << end_pos.front() + 1);
            answer.PushBack(connecting_scaffold->SubPath(start_pos.front() + 1, end_pos.front() + 1));
        }
    }
    return answer;
}

} // namespase

path_extend::PathContainer Launch(debruijn_graph::GraphPack const & gp,
                                  PathThreadingParams params,
                                  PathWithEdgePostionsContainer const & input_paths,
                                  std::vector<SeqString> & contigs,
                                  std::vector<std::string> const & paths_names,
                                  PathContainer const & scaffolds,
                                  size_t nthreads)
{
    params.use_scaffolds = !input_paths.empty();
    auto const & graph = gp.get<Graph>();
    GraphCoverageMap cover_map(graph);
    cover_map.AddPaths(scaffolds);

    PathContainer total_paths;
    for (size_t i = 0; i < input_paths.size(); ++i) {
        INFO("Processing path# " << i << " with " << input_paths[i].edge_set.size() << " edges");

        const auto& path_name = paths_names[i];
        auto & contig = *find_if(contigs.begin(), contigs.end(), [&path_name](SeqString const & contig){return contig.name == path_name;});
        SequenceCorrector corrector(graph, params, cover_map, input_paths[i], contig.seq);

        auto data = corrector.GetBestSequence();
        contig.seq = data.first;
        for (auto const & path : data.second) {
            auto p = make_unique<BidirectionalPath>(graph, std::move(path.path));
            auto cp = make_unique<BidirectionalPath>(p->Conjugate());
            total_paths.AddPair(p.release(), cp.release());
        }
    }
    INFO("DONE");
    return total_paths;
}