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

using namespace debruijn_graph;
using namespace path_extend;
using namespace std;

namespace {

class PathsCorrector {
    Graph const & graph;
    PathThreadingParams params;
public:

    PathsCorrector(Graph const & graph)
        : graph(graph)
    {}

    void InsertLCSP(EdgeId start, EdgeId end, const std::vector<std::vector<EdgeId>>& detected_paths, BidirectionalPath& res) const;

    void ProcessEdgePairUniqueOnly(EdgeId start, EdgeId end, size_t distance, BidirectionalPath& res) const;

    void ProcessPathUniqueOnly(const PathWithEdgePostions& path_info, const GraphCoverageMap& scaffold_coverage_map, vector<deque<const int*>>& refs_to_pos, PathContainer& result) const;

    void ExtendUniqueForward(BidirectionalPath& path) const;
};

std::vector<EdgeId> LCPS(const std::vector<std::vector<EdgeId>>& paths, bool prefix, size_t edges_used) {
    size_t min_size = paths.front().size();
    for (const auto& p : paths)
        min_size = std::min(min_size, p.size());
    min_size -= edges_used;

    size_t index = 0;
    bool all_equal = true;
    while (index < min_size && all_equal) {
        for (size_t i = 1; i < paths.size(); ++i) {
            auto e = paths.front()[prefix ? index : paths.front().size() - 1 - index];
            if (index >= paths[i].size() or e != paths[i][prefix ? index : paths[i].size() - 1 - index]) {
                all_equal = false;
                break;
            }
        }
        if (all_equal)
            ++index;
    }

    std::vector<EdgeId> result;
    for (size_t i = 0; i < index; ++i) {
        result.push_back(paths.front()[prefix ? i : paths.front().size() - 1 - i]);
    }
    return prefix ? result : std::vector<EdgeId>(result.rbegin(), result.rend());
}

void PathsCorrector::InsertLCSP(EdgeId start, EdgeId end, const std::vector<std::vector<EdgeId>>& detected_paths,
                                BidirectionalPath& res) const
{
    bool self_looped = true;
    for (size_t i = 1; i < detected_paths.size(); ++i) {
        if (std::find(detected_paths[i].begin(), detected_paths[i].end(), start) == detected_paths[i].end()
            and std::find(detected_paths[i].begin(), detected_paths[i].end(), end) == detected_paths[i].end()) {
            self_looped = false;
            break;
        }
    }
    if (self_looped) {
        DEBUG("Unique shortest path found between " << graph.int_id(start) << " and " << graph.int_id(end));
        res.PushBack(detected_paths.front());
    } else {
        auto lcp = LCPS(detected_paths, true, 0);
        auto lcs = LCPS(detected_paths, false, lcp.size());

        res.PushBack(lcp);
        if (!lcs.empty()) {
            auto dummy = lcs.front();
            res.PushBack(dummy, Gap(- int(graph.length(dummy))));
            res.PushBack(lcs);
        } else if (!res.Empty()) {
            auto dummy = res.Back();
            res.PushBack(dummy, Gap(- int(graph.length(dummy))));
        }
        DEBUG("Path with loop found between " << graph.int_id(start) << " and " << graph.int_id(end));
    }
}

bool in_bounds(double middle, double border, double border_factor) {
    VERIFY(math::le(0.0, border_factor) && math::le(border_factor, 1.0));
    return math::le(border * border_factor, middle) && math::le(middle, border * (2 - border_factor));
}

void PathsCorrector::ProcessEdgePairUniqueOnly(EdgeId start, EdgeId end, size_t distance, BidirectionalPath& res) const {
    VertexId target_vertex = graph.EdgeStart(end);
    VertexId start_vertex = graph.EdgeEnd(start);
    omnigraph::PathStorageCallback<Graph> path_storage(graph);
    omnigraph::ProcessPaths(graph, 0, params.max_distance, start_vertex, target_vertex, path_storage);

    res.Clear();
    const auto& detected_paths = path_storage.paths();

    if (path_storage.size() == 0) {
        DEBUG("No paths found between " << graph.int_id(start) << " and " << graph.int_id(end));
        return;
    }
    if (path_storage.size() == 1) {
        DEBUG("Unique path found between " << graph.int_id(start) << " and " << graph.int_id(end));
        res.PushBack(detected_paths.front());
        return;
    }
    std::vector<std::pair<size_t, size_t>> good_paths;
    if (distance > 0 && params.use_distances) {
        for (size_t i = 0; i < detected_paths.size(); ++i) {
            size_t difference = (size_t) abs((long long)(CumulativeLength(graph, detected_paths[i])) - (long long)(distance));
            if (math::le(double(difference), double(distance) * params.good_distance_coeff)) {
                good_paths.emplace_back(i, difference);
            }
        }
        if (good_paths.size() > 1) {
            std::sort(good_paths.begin(), good_paths.end(), 
                [](const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b) {
                    return a.second < b.second;
                });
            if (math::le(double(good_paths[0].second), double(good_paths[1].second) * params.best_of_good_coeff)) {
                good_paths.erase(good_paths.begin() + 1, good_paths.end());
            }
        }
    }

    if (good_paths.size() == 1) {
        size_t index = good_paths.front().first;
        DEBUG("Unique path with proper distance found between " << graph.int_id(start) << " and " << graph.int_id(end));
        res.PushBack(path_storage.paths()[index]);
        return;
    } 
    // if (good_paths.size() > 1) {
    //     std::vector<std::vector<EdgeId>> good_distance_paths;
    //     for (const auto& gp : good_paths) {
    //         good_distance_paths.push_back(std::move(detected_paths[gp.first]));
    //     }
    //     InsertLCSP(start, end, good_distance_paths, res);
    // } else {
    //     // InsertLCSP(start, end, detected_paths, res);
    // }
}


void ConnectWithScaffolds(EdgeId start, EdgeId end, const GraphCoverageMap& scaffold_coverage_map, BidirectionalPath& res) {
    res.Clear();
    DEBUG("Looking connecting scaffold");
    auto cover_start = scaffold_coverage_map.GetCoveringPaths(start);
    auto cover_end = scaffold_coverage_map.GetCoveringPaths(end);
    BidirectionalPathSet intersection;
    for (auto p : cover_start) {
        if (cover_end.count(p) > 0) {
            intersection.insert(p);
        }
    }
    DEBUG("Intersection found");
    if (intersection.size() == 1) {
        auto connecting_scaffold = *intersection.begin();
        auto start_pos = connecting_scaffold->FindAll(start);
        auto end_pos = connecting_scaffold->FindAll(end);
        if (start_pos.size() == 1 && end_pos.size() == 1 && start_pos.front() < end_pos.front()) {
            DEBUG("Found unique connecting scaffold " << start_pos.front() + 1 << " " << end_pos.front() + 1);
            res.PushBack(connecting_scaffold->SubPath(start_pos.front() + 1, end_pos.front() + 1));
        }
    }
}

void PathsCorrector::ProcessPathUniqueOnly(const PathWithEdgePostions& path_info,
                                           const GraphCoverageMap& scaffold_coverage_map,
                                           vector<deque<const int*>>& refs_to_pos,
                                           PathContainer& result) const
{
    auto current_path = make_unique<BidirectionalPath>(graph);
    deque<const int*> current_refs;
    const auto& path = path_info.edge_set;
    size_t i = 0;
    size_t j = 1;
    while (i < path.size() - 1) {
        const auto& alternative_starts = path[i];
        const auto& alternative_ends = path[j];
        DEBUG("Processing edges set #" << i << ", " << alternative_starts.size() << " starts, " << alternative_ends.size() << " ends");
        if (alternative_starts.empty() or alternative_ends.empty()) {
            ERROR("Zero alternatives given");
            ERROR("Path size " << path.size() << ", pos " << i << ", j" << j);
            ++i;
            j = i + 1;
            continue;
        }
        if (alternative_starts.size() != 1 or alternative_ends.size() != 1) {
            WARN("Multiple alternatives, expecting unique edges");
        }

        auto edge_start = alternative_starts.front();
        auto edge_end  = alternative_ends.front();
        DEBUG("Processing edges " << graph.int_id(edge_start) << " and " << graph.int_id(edge_end));

        if (current_path->Empty()) {
            current_path->PushBack(edge_start);
            current_refs.push_back(&path_info.positions[i]);
        }
        if (graph.EdgeStart(edge_end) == graph.EdgeEnd(edge_start)) {
            DEBUG("Adjacent edges, no need to look for paths");
            current_path->PushBack(edge_end);
            current_refs.push_back(&path_info.positions[j]);
            i = j;
            ++j;
            continue;
        }

        int distance = path_info.positions.empty() ? 0 :
            path_info.positions[j] - path_info.positions[i] - (int) graph.length(edge_start);
        if (distance < 0) {
            DEBUG("Negative distance between " << graph.int_id(edge_start) << " and " << graph.int_id(edge_end));
            distance = 0;
        }
        BidirectionalPath resulting_path(graph);

        ProcessEdgePairUniqueOnly(edge_start, edge_end, size_t(distance), resulting_path);

        if (params.use_scaffolds && resulting_path.Empty()) {
            ConnectWithScaffolds(edge_start, edge_end, scaffold_coverage_map, resulting_path);
            DEBUG("Resulting path size " << resulting_path.Size())
        }


        if (resulting_path.Empty()) {
            if (params.make_transitive_connections &&
                j - i < params.max_steps_forward && j < path.size() - 1) {
                DEBUG("No paths found, trying next " << j)
                ++j;
                continue;
            }
            if (1 < current_path->Size()) {
                auto cp = new BidirectionalPath(current_path->Conjugate());
                result.AddPair(current_path.release(), cp);
                refs_to_pos.push_back(std::move(current_refs));
            }
            current_refs.clear();
            j = i + 1;
            current_path = make_unique<BidirectionalPath>(graph, path[j].front());
            current_refs.push_back(&path_info.positions[j]);
        } else {
            for (size_t index = 0; index < resulting_path.Size(); ++index) {
                auto e = resulting_path[index];
                //Dirty fix for dummy edges that sign broken loop (see LCPS)
                if (resulting_path.GapAt(index).gap == - int(graph.length(e))) {
                    auto cp = new BidirectionalPath(current_path->Conjugate());
                    result.AddPair(current_path.release(), cp);
                    refs_to_pos.push_back(std::move(current_refs));
                    current_path = make_unique<BidirectionalPath>(graph);
                }
                else {
                    current_path->PushBack(e, resulting_path.GapAt(index));
                    current_refs.push_back(nullptr);
                }
            }
            if (current_path->Empty() || current_path->Back() != edge_end) {
                current_path->PushBack(edge_end);
                current_refs.push_back(&path_info.positions[j]);
            } else {
                current_refs.back() = &path_info.positions[j];
            }
        }
        i = j;
        ++j;
    }
    if (1 < current_path->Size()) {
        auto cp = new BidirectionalPath(current_path->Conjugate());
        result.AddPair(current_path.release(), cp);
        refs_to_pos.push_back(std::move(current_refs));
    }
}

void PathsCorrector::ExtendUniqueForward(BidirectionalPath& path) const {
    while (graph.OutgoingEdgeCount(graph.EdgeEnd(path.Back())) == 1) {
        auto e = *graph.OutgoingEdges(graph.EdgeEnd(path.Back())).begin();
        if (path.FindFirst(e) != -1)
            break;
        path.PushBack(e);
    }
}

void ExntendPathsEnds(PathsCorrector const & corrector, vector<deque<const int*>>& refs_to_pos, PathContainer& current_paths) {
    if (current_paths.size() == 0)
        return;

    corrector.ExtendUniqueForward(*current_paths.Get(0));
    while (refs_to_pos[0].size() != current_paths.Get(0)->Size())
        refs_to_pos[0].push_back(nullptr);
    
    corrector.ExtendUniqueForward(*current_paths.GetConjugate(0));
    while (refs_to_pos[0].size() != current_paths.Get(0)->Size())
        refs_to_pos[0].push_front(nullptr);

    size_t last = current_paths.size() - 1;

    corrector.ExtendUniqueForward(*current_paths.Get(last));
    while (refs_to_pos[last].size() != current_paths.Get(last)->Size())
        refs_to_pos[last].push_back(nullptr);

    corrector.ExtendUniqueForward(*current_paths.GetConjugate(last));
        while (refs_to_pos[last].size() != current_paths.Get(last)->Size())
        refs_to_pos[last].push_front(nullptr);

}

pair<long long, long long> CalcReplasedBounds(BidirectionalPath const & path, deque<const int*> const & refs) {
    auto left_point = find_if(refs.begin(), refs.end(), [](const int* p) {return p != nullptr;}) - refs.begin();
    auto right_point = refs.size() - (find_if(refs.rbegin(), refs.rend(), [](const int* p) {return p != nullptr;}) - refs.rbegin()) - 1;
    auto left_shift = path.Length() - path.LengthAt(left_point) + path.GapAt(left_point).gap;
    auto right_shift = path.LengthAt(right_point) - path.GapAt(right_point).gap + path.g().k();
    VERIFY(left_point == 0);
    VERIFY(right_point + 1 == path.Size());
    // VERIFY(left_point != right_point);
    // if (left_point == right_point) 
    //     VERIFY(left_shift + right_shift == path.Length() + path.g().k());
    return {*refs[left_point] - left_shift, *refs[right_point] + right_shift};
}

void BackMapping(PathContainer const & paths, vector<deque<const int*>> const & refs_to_pos, SeqString & contig) {
    stringstream ss;
    auto current_pos = std::numeric_limits<long long>::min();
    for (size_t i = 0; i < paths.size(); ++i) {
        auto const & path = paths[i];
        auto const & refs = refs_to_pos[i];
        ScaffoldSequenceMaker seq_maker(path.g());
        auto bounds = CalcReplasedBounds(path, refs);
        size_t should_be_dropped = 0;
        if (current_pos < bounds.first) {
            if (0 <= current_pos && current_pos < contig.seq.size())
                ss << contig.seq.substr(current_pos, bounds.first - current_pos);
        } else {
            should_be_dropped = current_pos - bounds.first;
            WARN("The new path prefix with len=" << should_be_dropped << " of " << path.Length() + path.g().k() << " would be dropped");
        }
        if (path.Length() + path.g().k() != bounds.second - bounds.first) {
            WARN("Replasing contig part with len=" << bounds.second - bounds.first << " using path with len=" << path.Length() + path.g().k());
        }
        if (should_be_dropped >= path.Length() + path.g().k()) {
            WARN("Oh no! Full path would be dropped");
            continue;
        }

        ss << seq_maker.MakeSequence(path).substr(should_be_dropped);
        current_pos = bounds.second;
    }
    current_pos = std::max(0ll, current_pos);
    if (current_pos < contig.seq.size())
        ss << contig.seq.substr(current_pos);
    contig.seq = ss.str();
}

} // namespase

path_extend::PathContainer Launch(debruijn_graph::GraphPack const & gp,
                                  PathThreadingParams params,
                                  PathWithEdgePostionsContainer const & input_paths,
                                  std::vector<SeqString> & contigs,
                                  std::vector<std::string> const & paths_names,
                                  PathContainer const & scaffolds)
{
    INFO("ExSPAnder repeat resolving tool started");
    params.use_scaffolds = !input_paths.empty();
    auto const & graph = gp.get<Graph>();
    GraphCoverageMap cover_map(graph);
    cover_map.AddPaths(scaffolds);

    PathsCorrector corrector(graph);
    PathContainer total_paths;
    size_t counter = 0;
    for (size_t i = 0; i < input_paths.size(); ++i) {
        const auto& path = input_paths[i];
        const auto& path_name = paths_names[i];
        vector<deque<const int*>> refs_to_pos;
        PathContainer result;
        INFO("Processing path# " << counter << " with " << path.edge_set.size() << " edges");
        if (path.edge_set.size() == 0) {
            INFO("Skipping empty path #" << counter++);
            continue;
        }

        if (path.edge_set.size() == 1) {
            INFO("Skipping almost empty path #" << counter++);
            continue;
            // auto p = new BidirectionalPath(graph, path.edge_set.front());
            // auto cp = new BidirectionalPath(p->Conjugate());
            // result.AddPair(p, cp);
            // refs_to_pos.emplace_back(1, &path.positions.front());
        } else {
            corrector.ProcessPathUniqueOnly(path, cover_map, refs_to_pos, result);
        }

        VERIFY(refs_to_pos.size() == result.size());
        for (size_t i = 0; i < refs_to_pos.size(); ++i) {
            VERIFY(refs_to_pos[i].size() == result[i].Size());
        }

        if (params.extend_unique_paths) {
            INFO("Extending tips" )
            ExntendPathsEnds(corrector, refs_to_pos, result);
        }

        VERIFY(refs_to_pos.size() == result.size());
        for (size_t i = 0; i < refs_to_pos.size(); ++i) {
            VERIFY(refs_to_pos[i].size() == result[i].Size());
        }

        for (auto const & p_ref : refs_to_pos) {
            bool hs = false;
            for (auto ref : p_ref)
                hs |= (ref != nullptr);
            VERIFY(hs);
        }
        
        INFO("Starting back mapping")
        auto & contig = *find_if(contigs.begin(), contigs.end(), [&path_name](SeqString const & contig){return contig.name == path_name;});
        BackMapping(result, refs_to_pos, contig);
        INFO("Back mapping was finished")

        total_paths.AddContainer(std::move(result));
    }
    INFO("DONE");
    return total_paths;
}
