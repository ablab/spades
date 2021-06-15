//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "helpers/common.hpp"
#include "helpers/replacer.hpp"
#include "filler_chooser.hpp"

#include "utils/logger/logger.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "common/assembly_graph/paths/path_processor.hpp"
#include "common/assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include "common/modules/path_extend/pe_utils.hpp"

#include <string>
#include <unordered_map>
#include <vector>
#include <list>
#include <fstream>

using namespace debruijn_graph;
using namespace path_extend;
using namespace std;

namespace {

struct PathWithBorderEdgesIndexies {
    SimpleBidirectionalPath path;
    size_t first_edge_index;
    size_t last_edge_index;
};
class SequenceCorrector {
    Graph const & graph;
    FillerChooser const & filler_chooser;
    PathThreadingParams params;
    GraphCoverageMap const & scaffold_coverage_map;
    PathWithEdgePostions path_info;
    std::string const & seq;
    size_t nthreads;
public:
    SequenceCorrector(Graph const & graph, FillerChooser const & filler_chooser, PathThreadingParams params, GraphCoverageMap const & scaffold_coverage_map, PathWithEdgePostions const & path_info, std::string const & seq, size_t nthreads)
        : graph(graph)
        , filler_chooser(filler_chooser)
        , params(params)
        , scaffold_coverage_map(scaffold_coverage_map)
        , path_info(path_info)
        , seq(seq)
        , nthreads(nthreads)
    {}

    pair<std::list<ReplaceInfo>, std::vector<PathWithBorderEdgesIndexies>> GetBestSequence() const;
private:
    pair<SimpleBidirectionalPath, size_t> GetNextPath(size_t start_pos, size_t & adjacent_edges_connected) const;

    boost::optional<pair<Gap, size_t>> FindFiller(size_t start_pos, size_t & adjacent_edges_connected) const;

    boost::optional<string> GetBestMachedPath(EdgeId start_edge, EdgeId end_edge, string const & ref) const;

    SimpleBidirectionalPath ConnectWithScaffolds(EdgeId start, EdgeId end) const;

    EdgeId GetEdge(size_t pos) const {
        return path_info.edges[pos];
    }

    long long GetStartPos(size_t pos) const { return path_info.start_positions[pos]; }
    long long GetEndPos(size_t pos) const { return path_info.end_positions[pos]; }

    size_t Size() const noexcept { return path_info.edges.size(); }

    /// @returns amount of nucls from the end of the 'start_pos' edge to the start of the 'end_pos' edge.
    /// @returns -1ull when the distance is negative.
    size_t GetDistance(size_t start_pos, size_t end_pos) const { 
        long long dist = GetStartPos(end_pos) - GetEndPos(start_pos);
        if (dist >= 0)
            return dist;
        DEBUG("Negative distance between " << graph.int_id(GetEdge(start_pos)) << " and " << graph.int_id(GetEdge(end_pos)));
        return -1ull;
        // return 0;
    }

    std::list<ReplaceInfo> ConvertToReplaceInfo(vector<PathWithBorderEdgesIndexies> const & paths) const;
};

std::list<ReplaceInfo> SequenceCorrector::ConvertToReplaceInfo(vector<PathWithBorderEdgesIndexies> const & paths) const {
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

pair<std::list<ReplaceInfo>, std::vector<PathWithBorderEdgesIndexies>> SequenceCorrector::GetBestSequence() const {
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

    return {ConvertToReplaceInfo(paths), std::move(paths)};
}

pair<SimpleBidirectionalPath, size_t> SequenceCorrector::GetNextPath(size_t start_pos, size_t & adjacent_edges_connected) const {
    SimpleBidirectionalPath current_path;
    current_path.PushBack(GetEdge(start_pos));
    while (start_pos + 1 < Size()) {
        auto filler = FindFiller(start_pos, adjacent_edges_connected);
        if (!filler.is_initialized())
            break;
        start_pos = filler->second;
        current_path.PushBack(GetEdge(start_pos), std::move(filler->first));
    }
    return {std::move(current_path), start_pos};
}

boost::optional<pair<Gap, size_t>> SequenceCorrector::FindFiller(size_t start_pos, size_t & adjacent_edges_connected) const {
    size_t end_pos = start_pos + 1;
    auto start_edge = GetEdge(start_pos);
    for (; end_pos - start_pos <= params.max_steps_forward && end_pos < Size(); ++end_pos) {
        auto end_edge  = GetEdge(end_pos);

        if (graph.EdgeStart(end_edge) == graph.EdgeEnd(start_edge) && GetStartPos(end_pos) + graph.k() == GetEndPos(start_pos)) {
            ++adjacent_edges_connected;
            return {{Gap(), end_pos}};
        }

        if (params.use_agressive_filling) {
            size_t distance = GetDistance(start_pos, end_pos);
            if (distance == -1ull)
                return {};
            distance += 2 * graph.k(); // add overlaping with start_edge and end_edge
            VERIFY(GetEndPos(start_pos) >= (long long)graph.k());
            size_t start_after_edge = GetEndPos(start_pos) - (long long)graph.k();
            auto gap_str = GetBestMachedPath(start_edge, end_edge, seq.substr(start_after_edge, distance));

            if (gap_str.is_initialized()) {
                assert(gap_str->size() >= graph.k());
                auto gap_size_in_nucl = gap_str->size() - 2 * graph.k(); // cut off overlaping with start_edge and end_edge
                gap_str = (gap_str->size() > 2 * graph.k() ? gap_str->substr(graph.k(), gap_size_in_nucl) : "");
                return {{Gap(std::move(*gap_str), gap_size_in_nucl+graph.k()), end_pos}};
            }
        }

        /// TODO: think how to correct ConnectWithScaffolds
        if (params.use_scaffolds) {
            assert(false);
            // ConnectWithScaffolds(start_edge, end_edge);
        }
    }
    return {};
}

boost::optional<string> SequenceCorrector::GetBestMachedPath(EdgeId start_edge, EdgeId end_edge, string const & ref) const {
    if (ref.empty())
        return {};

    auto distance = ref.size();
    VertexId target_vertex = graph.EdgeStart(end_edge);
    VertexId start_vertex = graph.EdgeEnd(start_edge);
    omnigraph::PathStorageCallback<Graph> path_storage(graph);
    auto min_len = static_cast<size_t>(std::floor((double)distance*(1-params.good_distance_coeff)));
    auto max_len = static_cast<size_t>(std::ceil ((double)distance*(1+params.good_distance_coeff)));
    auto error = omnigraph::ProcessPaths(graph, min_len, max_len, start_vertex, target_vertex, path_storage);

    if (error) {
        WARN("WOW! omnigraph::ProcessPaths failed!, error: " << error);
        return {};
    }

    const auto& detected_paths = path_storage.paths();
    return filler_chooser(detected_paths, ref);
}

SimpleBidirectionalPath SequenceCorrector::ConnectWithScaffolds(EdgeId start, EdgeId end) const {
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
        
        AlignerFiller filler_chooser(graph, params.score_domination_coeff);
        // DistanceFiller filler_chooser(graph, params.good_distance_coeff, params.best_of_good_coeff);
        SequenceCorrector corrector(graph, filler_chooser, params, cover_map, input_paths[i], contig.seq, nthreads);

        auto data = corrector.GetBestSequence();
        contig.seq = ReplaceAndDump(contig.seq, std::move(data.first), contig.name);
        PathContainer result;
        for (auto const & path : data.second)
            result.Add(BidirectionalPath::create(graph, std::move(path.path)));

        // #pragma omp critical
        {
            total_paths.AddContainer(std::move(result));
        }
    }

    INFO("DONE");
    return total_paths;
}
