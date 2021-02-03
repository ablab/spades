//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "common.hpp"
#include "replacer.hpp"

#include "utils/logger/logger.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "common/assembly_graph/core/basic_graph_stats.hpp"
#include "common/assembly_graph/core/graph.hpp"
#include "common/assembly_graph/paths/path_processor.hpp"
#include "common/assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include "common/modules/path_extend/pe_utils.hpp"
#include "ssw/ssw_cpp.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <list>
#include <fstream>
#include <boost/optional.hpp>

using namespace debruijn_graph;
using namespace path_extend;
using namespace std;

namespace {

struct PathWithBorderEdgesIndexies {
    SimpleBidirectionalPath path;
    size_t first_edge_index;
    size_t last_edge_index;
};

enum class PathFiller {
    UseDistance,
    UseAligner
};

template<PathFiller filler_mode>
class SequenceCorrector {
    Graph const & graph;
    PathThreadingParams params;
    GraphCoverageMap const & scaffold_coverage_map;
    PathWithEdgePostions path_info;
    std::string const & seq;
    size_t nthreads;
public:
    SequenceCorrector(Graph const & graph, PathThreadingParams params, GraphCoverageMap const & scaffold_coverage_map, PathWithEdgePostions const & path_info, std::string const & seq, size_t nthreads)
        : graph(graph)
        , params(params)
        , scaffold_coverage_map(scaffold_coverage_map)
        , path_info(path_info)
        , seq(seq)
        , nthreads(nthreads)
    {}

    pair<std::string, std::vector<PathWithBorderEdgesIndexies>> GetBestSequence() const;
private:
    pair<SimpleBidirectionalPath, size_t> GetNextPath(size_t start_pos, size_t & adjacent_edges_connected) const;

    boost::optional<pair<SimpleBidirectionalPath, size_t>> FindFiller(size_t start_pos, size_t & adjacent_edges_connected) const;

    SimpleBidirectionalPath GetBestMachedPathUsingDistance(EdgeId start_edge, EdgeId end_edge, size_t distance, string const & ref) const;
    SimpleBidirectionalPath GetBestMachedPathUsingAligner(EdgeId start_edge, EdgeId end_edge, size_t distance, string const & ref) const;

    SimpleBidirectionalPath ConnectWithScaffolds(EdgeId start, EdgeId end) const;

    EdgeId GetEdge(size_t pos) const {
        return path_info.edges[pos];
    }

    long long GetStartPos(size_t pos) const { return path_info.start_positions[pos]; }
    long long GetEndPos(size_t pos) const { return path_info.end_positions[pos]; }

    size_t Size() const noexcept { return path_info.edges.size(); }

    size_t GetDistance(size_t start_pos, size_t end_pos) const { 
        long long dist = GetStartPos(end_pos) - GetStartPos(start_pos) + graph.k();
        if (dist >= 0)
            return dist;
        DEBUG("Negative distance between " << graph.int_id(GetEdge(start_pos)) << " and " << graph.int_id(GetEdge(end_pos)));
        return -1ull;
        // return 0;
    }

    std::list<ReplaceInfo> ConvertToReplaceInfo(vector<PathWithBorderEdgesIndexies> const & paths) const;
};

template<PathFiller filler_mode>
std::list<ReplaceInfo> SequenceCorrector<filler_mode>::ConvertToReplaceInfo(vector<PathWithBorderEdgesIndexies> const & paths) const {
    list<ReplaceInfo> mapping_info;
    ScaffoldSequenceMaker seq_maker(graph);
    for (size_t i = 0; i < paths.size(); ++i) {
        ReplaceInfo info(seq_maker.MakeSequence(*BidirectionalPath::create(graph, paths[i].path)));
        auto start_pos = GetStartPos(paths[i].first_edge_index);
        auto end_pos = GetEndPos(paths[i].last_edge_index);
        if (start_pos < 0)
            info.drop_from_head = -start_pos;
        if (seq.size() < end_pos)
            info.drop_from_tail = end_pos - seq.size();
        if (info.ShouldBeDropped())
            continue;
        info.contig_start_pos = start_pos + info.drop_from_head;
        info.contig_end_pos = end_pos - info.drop_from_tail;
        mapping_info.push_back(std::move(info));
    }
    return mapping_info;
}

template<PathFiller filler_mode>
pair<std::string, std::vector<PathWithBorderEdgesIndexies>> SequenceCorrector<filler_mode>::GetBestSequence() const {
    std::vector<PathWithBorderEdgesIndexies> paths;
    size_t adjacent_edges_connected = 0;
    size_t start_pos = 0;
    while (start_pos < Size()) {
        auto path = GetNextPath(start_pos, adjacent_edges_connected);
        if (path.first.Size() > 0) {
            paths.push_back({std::move(path.first), start_pos, path.second});
        }
        start_pos = path.second + 1;
    }
    INFO("Made " << paths.size() << " path" << (paths.size() != 1 ? "s" : ""));
    INFO("Connected " << adjacent_edges_connected << " adjacent edges pair" << (adjacent_edges_connected == 1 ? "" : "s"));

    return {Replace(seq, ConvertToReplaceInfo(paths)), std::move(paths)};
}

template<PathFiller filler_mode>
pair<SimpleBidirectionalPath, size_t> SequenceCorrector<filler_mode>::GetNextPath(size_t start_pos, size_t & adjacent_edges_connected) const {
    SimpleBidirectionalPath current_path;
    current_path.PushBack(GetEdge(start_pos));
    while (start_pos + 1 < Size()) {
        auto path = FindFiller(start_pos, adjacent_edges_connected);
        if (!path.is_initialized())
            break;
        current_path.PushBack(std::move(path->first));
        start_pos = path->second;
        if (current_path.Back() != GetEdge(start_pos))
            current_path.PushBack(GetEdge(start_pos));
    }
    return {std::move(current_path), start_pos};
}

template<PathFiller filler_mode>
boost::optional<pair<SimpleBidirectionalPath, size_t>> SequenceCorrector<filler_mode>::FindFiller(size_t start_pos, size_t & adjacent_edges_connected) const {
    size_t end_pos = start_pos + 1;
    auto start_edge = GetEdge(start_pos);
    for (; end_pos - start_pos <= params.max_steps_forward && end_pos < Size(); ++end_pos) {
        SimpleBidirectionalPath current_path;
        auto end_edge  = GetEdge(end_pos);

        if (graph.EdgeStart(end_edge) == graph.EdgeEnd(start_edge) && GetStartPos(end_pos) + graph.k() == GetEndPos(start_pos)) {
            ++adjacent_edges_connected;
            return {{std::move(current_path), end_pos}};
        }

        if (params.use_agressive_filling) {
            size_t distance = GetDistance(start_pos, end_pos);
            if (distance == -1ull)
                return {};
            VERIFY(GetStartPos(end_pos) >= (long long)distance);
            size_t start_after_edge = GetStartPos(end_pos) - distance;
            if (filler_mode == PathFiller::UseDistance)
                current_path = GetBestMachedPathUsingDistance(start_edge, end_edge, distance, seq.substr(start_after_edge, distance));
            else
                current_path = GetBestMachedPathUsingAligner(start_edge, end_edge, distance, seq.substr(start_after_edge, distance));
        }

        if (params.use_scaffolds && current_path.Empty())
            current_path = ConnectWithScaffolds(start_edge, end_edge);
 
        if (!current_path.Empty()) {
            return {{std::move(current_path), end_pos}};
        }
    }
    return {};
}

template<PathFiller filler_mode>
SimpleBidirectionalPath SequenceCorrector<filler_mode>::GetBestMachedPathUsingDistance(EdgeId start_edge, EdgeId end_edge, size_t distance, string const &) const {
    ScaffoldSequenceMaker seq_maker(graph);
    VertexId target_vertex = graph.EdgeStart(end_edge);
    VertexId start_vertex = graph.EdgeEnd(start_edge);
    omnigraph::PathStorageCallback<Graph> path_storage(graph);
    auto min_len = static_cast<size_t>(std::floor((double)distance*(1-params.good_distance_coeff)));
    auto max_len = static_cast<size_t>(std::ceil ((double)distance*(1+params.good_distance_coeff)));
    omnigraph::ProcessPaths(graph, min_len, max_len, start_vertex, target_vertex, path_storage);
    const auto& detected_paths = path_storage.paths();
    
    SimpleBidirectionalPath answer;
    vector<pair<size_t, size_t>> scores;
    if (distance > 0) {
        for (size_t i = 0; i < detected_paths.size(); ++i) {
            auto path_len = CumulativeLength(graph, detected_paths[i]);
            size_t difference = (size_t) abs((long long)(path_len) - (long long)(distance));
            if (math::le(double(difference), double(distance) * params.good_distance_coeff))
                scores.emplace_back(difference, i);
        }
        auto by_dist = [](const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b) { return a.first < b.first; };
        std::sort(scores.begin(), scores.end(), by_dist);
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

template<PathFiller filler_mode>
SimpleBidirectionalPath SequenceCorrector<filler_mode>::GetBestMachedPathUsingAligner(EdgeId start_edge, EdgeId end_edge, size_t distance, string const & ref) const {
    ScaffoldSequenceMaker seq_maker(graph);
    StripedSmithWaterman::Aligner aligner(1, 1, 1, 1);
    aligner.SetReferenceSequence(ref.data(), (int)ref.size());
    VertexId target_vertex = graph.EdgeStart(end_edge);
    VertexId start_vertex = graph.EdgeEnd(start_edge);
    omnigraph::PathStorageCallback<Graph> path_storage(graph);
    auto min_len = static_cast<size_t>(std::floor((double)distance*(1-params.good_distance_coeff)));
    auto max_len = static_cast<size_t>(std::ceil ((double)distance*(1+params.good_distance_coeff)));
    omnigraph::ProcessPaths(graph, min_len, max_len, start_vertex, target_vertex, path_storage);
    const auto& detected_paths = path_storage.paths();
    
    SimpleBidirectionalPath answer;
    vector<pair<size_t, size_t>> scores;
    if (distance > 0) {
        // #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
        #pragma omp parallel for schedule(runtime)
        for (size_t i = 0; i < detected_paths.size(); ++i) {
            StripedSmithWaterman::Alignment alignment;
            auto path = BidirectionalPath::create(graph, detected_paths[i]);
            auto query_seq = seq_maker.MakeSequence(*path);
            StripedSmithWaterman::Filter filter(false, false, 0, uint8_t(1 + std::max(ref.size(), query_seq.size())));
            
            if (aligner.Align(query_seq.data(), filter, &alignment)) {
                #pragma omp critical
                {
                    scores.emplace_back(alignment.sw_score, i);
                }
            }
        }

        if (scores.size() > 1) {
            auto by_score = [](const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b) { return a.first > b.first; };
            std::sort(scores.begin(), scores.end(), by_score);
        }
    }

    if (!scores.empty() && scores.front().first > 0) {
        size_t index = scores.front().second;
        answer.PushBack(detected_paths[index]);
    }
    return answer;
}

template<PathFiller filler_mode>
SimpleBidirectionalPath SequenceCorrector<filler_mode>::ConnectWithScaffolds(EdgeId start, EdgeId end) const {
    auto cover_start = scaffold_coverage_map.GetCoveringPaths(start);
    auto cover_end = scaffold_coverage_map.GetCoveringPaths(end);
    BidirectionalPathSet intersection;
    for (auto p : cover_start) {
        if (cover_end.count(p) > 0)
            intersection.insert(p);
    }
    SimpleBidirectionalPath answer;
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
                                  PathContainer const & scaffolds,
                                  size_t nthreads)
{
    params.use_scaffolds &= !input_paths.empty();
    auto const & graph = gp.get<Graph>();
    GraphCoverageMap cover_map(graph);
    cover_map.AddPaths(scaffolds);

    #ifdef GOOD_NAME
    for (auto & contig : contigs) {
        if (contig.name != GOOD_NAME) {
            contig.seq = "";
        }
    }
    #endif

    PathContainer total_paths;
    // #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    for (size_t i = 0; i < input_paths.size(); ++i) {
        const auto& path_name = input_paths[i].path_name;
        #ifdef GOOD_NAME
        if (path_name != GOOD_NAME) {
            continue;
        }
        #endif
        auto & contig = *find_if(contigs.begin(), contigs.end(), [&path_name](SeqString const & contig){return contig.name == path_name;});
        INFO("Processing path [" << contig.name << "] # " << i + 1 << " (of " << input_paths.size() << ") with " << input_paths[i].edges.size() << " edges");
        SequenceCorrector<PathFiller::UseDistance> corrector(graph, params, cover_map, input_paths[i], contig.seq, nthreads);

        auto data = corrector.GetBestSequence();
        contig.seq = data.first;
        PathContainer result;
        for (auto const & path : data.second) {
            auto p = BidirectionalPath::create(graph, std::move(path.path));
            auto cp = BidirectionalPath::clone_conjugate(p);
            result.AddPair(move(p), move(cp));
        }

        // #pragma omp critical
        {
            total_paths.AddContainer(std::move(result));
        }
    }

    INFO("DONE");
    return total_paths;
}