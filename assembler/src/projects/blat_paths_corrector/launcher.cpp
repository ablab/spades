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
    GappedPath path;
    size_t first_edge_index;
    size_t last_edge_index;
};

struct BackMappingInfo {
    BidirectionalPath path;
    size_t contig_start_pos; // inclusive
    size_t contig_end_pos;   // exclusive
    size_t drop_from_head;
    size_t drop_from_tail;

    BackMappingInfo(BidirectionalPath path)
        : path(std::move(path))
        , contig_start_pos(0)
        , contig_end_pos(0)
        , drop_from_head(0)
        , drop_from_tail(0)
    {}

    size_t LengthInNucls() const {
        return path.Length() + path.g().k();
    }

    bool ShouldBeDropped() const {
        return LengthInNucls() <= drop_from_head + drop_from_tail;
    }

    size_t ResultLengthInNucls() const {
        VERIFY_MSG(LengthInNucls() > drop_from_head + drop_from_tail, "LoL? LengthInNucls() = " << LengthInNucls() << ", drop_from_head = " << drop_from_head << ", drop_from_tail = " << drop_from_tail);
        return LengthInNucls() - drop_from_head - drop_from_tail;
    }
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
    {
        DropOverlappedEdges();
    }

    pair<std::string, std::vector<PathWithBorderEdgesIndexies>> GetBestSequence() const;
private:
    pair<GappedPath, size_t> GetNextPath(size_t start_pos) const;
    
    boost::optional<pair<GappedPath, size_t>> FindFiller(size_t start_pos) const;

    GappedPath GetBestMachedPathUsingDistance(EdgeId start_edge, EdgeId end_edge, size_t distance, string const & ref) const;
    GappedPath GetBestMachedPathUsingAligner(EdgeId start_edge, EdgeId end_edge, size_t distance, string const & ref) const;

    GappedPath ConnectWithScaffolds(EdgeId start, EdgeId end) const;

    size_t EdgeLen(EdgeId edge) const { return graph.length(edge); }
    size_t EdgeLenInNucl(EdgeId edge) const { return EdgeLen(edge) + graph.k(); }

    EdgeId GetEdge(size_t pos) const { 
        VERIFY_MSG(path_info.edge_set[pos].size() == 1, "alternatives are not supported");
        return path_info.edge_set[pos].front(); 
    }

    long long GetPos(size_t pos) const { return path_info.positions[pos]; }

    size_t Size() const noexcept { return path_info.edge_set.size(); }
    
    size_t GetDistance(size_t start_pos, size_t end_pos) const { 
        long long start_after_edge = GetPos(start_pos) + EdgeLen(GetEdge(start_pos));
        if (start_after_edge <= GetPos(end_pos))
            return GetPos(end_pos) - start_after_edge;
        DEBUG("Negative distance between " << graph.int_id(GetEdge(start_pos)) << " and " << graph.int_id(GetEdge(end_pos)));
        return -1ull;
        // return 0;
    }

    void DropOverlappedEdges() {
        vector<bool> is_bad_edge(Size(), false);
        for (size_t i = 0; i + 1 < Size(); ++i) {
            if (GetPos(i + 1) < GetPos(i) + (long long)EdgeLenInNucl(GetEdge(i))) {
                is_bad_edge[i] = true;
                is_bad_edge[i + 1] = true;
            }
        }
        PathWithEdgePostions new_path_info;
        for (size_t i = 0; i < Size(); ++i) {
            if (!is_bad_edge[i]) {
                new_path_info.edge_set.push_back(std::move(path_info.edge_set[i]));
                new_path_info.positions.push_back(path_info.positions[i]);
            }
        }
        auto dropped_pats_cnt = Size();
        path_info = std::move(new_path_info);
        dropped_pats_cnt -= Size();
        if (dropped_pats_cnt > 0)
            WARN(dropped_pats_cnt << " edge" << (dropped_pats_cnt != 1 ? "s" : "") << " would be dropped");
    }

    std::string BackMappingDropCurrentSuffix(vector<PathWithBorderEdgesIndexies> const & paths) const;
    std::string BackMappingDropNextPreffix(vector<PathWithBorderEdgesIndexies> const & paths) const;
    std::list<BackMappingInfo> ConvertToBackMappingInfo(vector<PathWithBorderEdgesIndexies> const & paths) const;
    std::string MakeSeq(list<BackMappingInfo> const & mapping_info) const;
};

template<PathFiller filler_mode>
std::list<BackMappingInfo> SequenceCorrector<filler_mode>::ConvertToBackMappingInfo(vector<PathWithBorderEdgesIndexies> const & paths) const {
    list<BackMappingInfo> mapping_info;
    for (size_t i = 0; i < paths.size(); ++i) {
        BackMappingInfo info(BidirectionalPath(graph, paths[i].path));
        auto right_shift = static_cast<long long>(EdgeLenInNucl(info.path.Back()));
        auto start_pos = GetPos(paths[i].first_edge_index);
        auto end_pos = GetPos(paths[i].last_edge_index) + right_shift;
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
std::string SequenceCorrector<filler_mode>::MakeSeq(list<BackMappingInfo> const & mapping_info) const {
    size_t current_pos = 0;
    stringstream ss;
    ScaffoldSequenceMaker seq_maker(graph);
    for (auto const & path : mapping_info) {
        if (path.drop_from_head > 0)
            WARN("The new path prefix with len=" << path.drop_from_head << " of " << path.LengthInNucls() << " would be dropped");
        if (path.drop_from_tail > 0)
            WARN("The new path suffix with len=" << path.drop_from_tail << " of " << path.LengthInNucls() << " would be dropped");
        ss << seq.substr(current_pos, path.contig_start_pos - current_pos);
        ss << seq_maker.MakeSequence(path.path).substr(path.drop_from_head, path.ResultLengthInNucls());
        current_pos = path.contig_end_pos;
    }
    ss << seq.substr(current_pos);
    return ss.str();
}

template<PathFiller filler_mode>
std::string SequenceCorrector<filler_mode>::BackMappingDropNextPreffix(vector<PathWithBorderEdgesIndexies> const & paths) const {
    auto mapping_info = ConvertToBackMappingInfo(paths);
    if (mapping_info.empty())
        return seq;

    auto current_path = mapping_info.begin();
    while (true) {
        auto next_path = std::next(current_path);
        if (next_path == mapping_info.end())
            break;
        if (next_path->contig_start_pos < current_path->contig_end_pos) {
            next_path->drop_from_head += current_path->contig_end_pos - next_path->contig_start_pos;
            next_path->contig_start_pos = current_path->contig_end_pos;
        }
        if (next_path->ShouldBeDropped()) {
            mapping_info.erase(next_path);
            continue;
        }
        ++current_path;
    }
    if (mapping_info.size() != paths.size()) {
        size_t dropped_pats_cnt = paths.size() - mapping_info.size();
        WARN("Oh no! Full " << dropped_pats_cnt << " path" << (dropped_pats_cnt != 1 ? "s" : "") << " would be dropped");
    }

    return MakeSeq(mapping_info);
}

template<PathFiller filler_mode>
std::string SequenceCorrector<filler_mode>::BackMappingDropCurrentSuffix(vector<PathWithBorderEdgesIndexies> const & paths) const {
    auto mapping_info = ConvertToBackMappingInfo(paths);
    if (mapping_info.empty())
        return seq;

    auto current_path = mapping_info.rbegin();
    while (true) {
        auto prev_path = std::next(current_path);
        if (prev_path == mapping_info.rend())
            break;
        if (current_path->contig_start_pos < prev_path->contig_end_pos) {
            prev_path->drop_from_tail += prev_path->contig_end_pos - current_path->contig_start_pos;
            prev_path->contig_end_pos = current_path->contig_start_pos;
        }
        auto const & base = std::prev(prev_path.base());
        VERIFY(base->path.GetId() == prev_path->path.GetId());
        if (prev_path->ShouldBeDropped()) {
            mapping_info.erase(std::prev(prev_path.base()));
            continue;
        }
        ++current_path;
    }
    if (mapping_info.size() != paths.size()) {
        size_t dropped_pats_cnt = paths.size() - mapping_info.size();
        WARN("Oh no! Full " << dropped_pats_cnt << " path" << (dropped_pats_cnt != 1 ? "s" : "") << " would be dropped");
    }

    return MakeSeq(mapping_info);
}

template<PathFiller filler_mode>
pair<std::string, std::vector<PathWithBorderEdgesIndexies>> SequenceCorrector<filler_mode>::GetBestSequence() const {
    std::vector<PathWithBorderEdgesIndexies> paths;
    size_t start_pos = 0;
    while (start_pos < Size()) {
        auto path = GetNextPath(start_pos);
        if (path.first.Size() > 0) {
            paths.push_back({std::move(path.first), start_pos, path.second});
        }
        start_pos = path.second + 1;
    }
    INFO("Made " << paths.size() << " path" << (paths.size() != 1 ? "s" : ""));
    return {BackMappingDropCurrentSuffix(paths), std::move(paths)};
    // return {BackMappingDropNextPreffix(paths), std::move(paths)};
}

template<PathFiller filler_mode>
pair<GappedPath, size_t> SequenceCorrector<filler_mode>::GetNextPath(size_t start_pos) const {
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

template<PathFiller filler_mode>
boost::optional<pair<GappedPath, size_t>> SequenceCorrector<filler_mode>::FindFiller(size_t start_pos) const {
    size_t end_pos = start_pos + 1;
    auto start_edge = GetEdge(start_pos);
    for (; end_pos - start_pos <= params.max_steps_forward && end_pos < Size(); ++end_pos) {
        GappedPath current_path;
        auto end_edge  = GetEdge(end_pos);

        if (graph.EdgeStart(end_edge) == graph.EdgeEnd(start_edge)) {
            return {{std::move(current_path), end_pos}};
        }

        if (params.use_agressive_filling) {
            size_t distance = GetDistance(start_pos, end_pos);
            if (distance == -1ull)
                return {};
            VERIFY(GetPos(end_pos) >= (long long)distance);
            size_t start_after_edge = GetPos(end_pos) - distance;
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
GappedPath SequenceCorrector<filler_mode>::GetBestMachedPathUsingDistance(EdgeId start_edge, EdgeId end_edge, size_t distance, string const &) const {
    ScaffoldSequenceMaker seq_maker(graph);
    VertexId target_vertex = graph.EdgeStart(end_edge);
    VertexId start_vertex = graph.EdgeEnd(start_edge);
    omnigraph::PathStorageCallback<Graph> path_storage(graph);
    auto min_len = static_cast<size_t>(std::floor((double)distance*(1-params.good_distance_coeff)));
    auto max_len = static_cast<size_t>(std::ceil ((double)distance*(1+params.good_distance_coeff)));
    omnigraph::ProcessPaths(graph, min_len, max_len, start_vertex, target_vertex, path_storage);
    const auto& detected_paths = path_storage.paths();
    
    GappedPath answer;
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
GappedPath SequenceCorrector<filler_mode>::GetBestMachedPathUsingAligner(EdgeId start_edge, EdgeId end_edge, size_t distance, string const & ref) const {
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
    
    GappedPath answer;
    vector<pair<size_t, size_t>> scores;
    if (distance > 0) {
        // #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
        #pragma omp parallel for schedule(runtime)
        for (size_t i = 0; i < detected_paths.size(); ++i) {
            StripedSmithWaterman::Alignment alignment;
            BidirectionalPath path(graph, detected_paths[i]);
            auto query_seq = seq_maker.MakeSequence(path);
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
GappedPath SequenceCorrector<filler_mode>::ConnectWithScaffolds(EdgeId start, EdgeId end) const {
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
        const auto& path_name = paths_names[i];
        #ifdef GOOD_NAME
        if (path_name != GOOD_NAME) {
            continue;
        }
        #endif
        INFO("Processing path# " << i + 1 << " (of " << input_paths.size() << ") with " << input_paths[i].edge_set.size() << " edges");
        auto & contig = *find_if(contigs.begin(), contigs.end(), [&path_name](SeqString const & contig){return contig.name == path_name;});
        SequenceCorrector<PathFiller::UseDistance> corrector(graph, params, cover_map, input_paths[i], contig.seq, nthreads);

        auto data = corrector.GetBestSequence();
        contig.seq = data.first;
        PathContainer result;
        for (auto const & path : data.second) {
            auto p = make_unique<BidirectionalPath>(graph, std::move(path.path));
            auto cp = make_unique<BidirectionalPath>(p->Conjugate());
            result.AddPair(p.release(), cp.release());
        }

        // #pragma omp critical
        {
            total_paths.AddContainer(std::move(result));
        }
    }

    INFO("DONE");
    return total_paths;
}