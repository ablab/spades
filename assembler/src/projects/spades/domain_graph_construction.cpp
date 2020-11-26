//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "domain_graph_construction.hpp"
#include "assembly_graph/core/construction_helper.hpp"
#include "assembly_graph/paths/path_processor.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "modules/path_extend/pe_utils.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "domain_graph.hpp"
#include "domain_matcher.hpp"
#include "utils/filesystem/path_helper.hpp"

namespace debruijn_graph {

template<class Graph>
class SetOfForbiddenEdgesPathChooser : public omnigraph::PathProcessor<Graph>::Callback {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph &g_;
    std::set<std::vector<EdgeId>> forbidden_edges_;
    std::vector<EdgeId> answer_path_;

    bool CheckCoverageDiff(const std::vector<EdgeId> &path) const {
        double min_coverage = std::numeric_limits<double>::max();
        double max_coverage = std::numeric_limits<double>::min();

        for (size_t i = 0; i < path.size(); ++i) {
            min_coverage = std::min(min_coverage, g_.coverage(path[i]));
            max_coverage = std::max(max_coverage, g_.coverage(path[i]));
        }
        return math::ge(50.0, max_coverage/min_coverage);
    }

    bool IsNewPathBetter(const path_extend::BidirectionalPath &current, const path_extend::BidirectionalPath &candidate) const {
        int current_length = (int)current.Length();
        int candidate_length = (int)candidate.Length();
        return candidate_length < current_length;
    }

public:
    SetOfForbiddenEdgesPathChooser(const Graph &g, const std::set<std::vector<EdgeId>> &forbidden_edges)
            : g_(g), forbidden_edges_(forbidden_edges) {}

    void HandleReversedPath(const std::vector<EdgeId> &reversed_path) override {
        std::vector<EdgeId> forward_path = reversed_path;
        std::reverse(forward_path.begin(), forward_path.end());

        for (const auto& forbidden_path : forbidden_edges_) {
            if (forbidden_path.empty())
                continue;

            bool cross_start = std::find(std::begin(forward_path), std::end(forward_path), forbidden_path.front()) != forward_path.end();
            bool cross_end = find(std::begin(forward_path), std::end(forward_path), forbidden_path.back()) != forward_path.end();

            if (cross_end && cross_start)
                return;
        }

        if (answer_path_.empty()) {
            if (!CheckCoverageDiff(forward_path))
               return;

            answer_path_ = forward_path;
            return;
        }

        if (!CheckCoverageDiff(forward_path))
            return;

        auto current = path_extend::BidirectionalPath::create(g_, answer_path_);
        auto candidate = path_extend::BidirectionalPath::create(g_, forward_path);
        if (IsNewPathBetter(*current, *candidate)) {
            answer_path_ = forward_path;
        }
    }

    void reset() {
        answer_path_.clear();
    }

    const std::vector<EdgeId>& answer() {
        return answer_path_;
    }
};

class DomainGraphConstructor {
public:
    DomainGraphConstructor(debruijn_graph::GraphPack &gp)
            : gp_(gp), domain_graph_(gp_.get<Graph>()) {}

    nrps::DomainGraph& ConstructGraph(const nrps::ContigAlnInfo &info) {
        INFO("Domain graph construction started");
        INFO("Constructing nodes...");
        ConstructNodes(info);
        INFO("Constructing strong edges...");
        ConstructStrongEdges();
        INFO("Constructing weak edges...");
        ConstructWeakEdges();
        INFO("Domain graph construction constructed, total vertices: " << domain_graph_.size() << ", edges: " << domain_graph_.e_size());
        return domain_graph_;
    }

private:

    void ConnectWithWeakEdge(VertexId v1, VertexId v2,
                             SetOfForbiddenEdgesPathChooser<Graph> &chooser) {
        const auto &g = gp_.get<Graph>();
        DEBUG("Trying to connect " << domain_graph_.GetVertexName(v1) << " and " << domain_graph_.GetVertexName(v2) << " with weak edge");
        int last_mapping = (int)(g.length(domain_graph_.mapping_path(v1).back().first) - domain_graph_.mapping_path(v1).end_pos());
        int first_mapping = (int)domain_graph_.mapping_path(v2).start_pos();

        if (5000 < last_mapping + first_mapping)
            return;

        int min_len = 0;
        if (g.EdgeEnd(domain_graph_.domain_edges(v1).back()) != g.EdgeStart(domain_graph_.domain_edges(v2).front())) {
            DEBUG("Trying to find paths from " << g.EdgeEnd(domain_graph_.domain_edges(v1).back()) << " to " << g.EdgeStart(domain_graph_.domain_edges(v2).front()));
            ProcessPaths(g, min_len, 4000 - last_mapping - first_mapping, g.EdgeEnd(domain_graph_.domain_edges(v1).back()), g.EdgeStart(domain_graph_.domain_edges(v2).front()), chooser);
            if (!chooser.answer().empty()) {
                DEBUG("Path was found");
                auto p = path_extend::BidirectionalPath::create(g, chooser.answer());
                DEBUG("Start vertex: " << g.EdgeStart(p->Front()).int_id());
                DEBUG("End vertex: " << g.EdgeEnd(p->Back()).int_id());
                DEBUG("Path:");
                p->PrintDEBUG();
                domain_graph_.AddEdge(v1, v2, false, chooser.answer(), p->Length() + last_mapping + first_mapping);
            }
            else {
                DEBUG("Path was not found");
            }
        } else {
            domain_graph_.AddEdge(v1, v2, false, std::vector<EdgeId>(), last_mapping + first_mapping);
        }
        chooser.reset();
    }

    std::vector<VertexId> VerticesReachedFrom(VertexId start_vertex) {
        auto bounded_dijkstra = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(gp_.get<Graph>(),
                                                                                        4000, 10000);
        bounded_dijkstra.Run(start_vertex);
        TRACE("Reached vertices size - " << bounded_dijkstra.ReachedVertices());
        return bounded_dijkstra.ReachedVertices();
    }

    void ConstructWeakEdges() {
        const auto &g = gp_.get<Graph>();

        std::set<std::vector<EdgeId>> forbidden_edges;
        for (VertexId v : domain_graph_.vertices())
            forbidden_edges.insert(domain_graph_.mapping_path(v).simple_path());

        SetOfForbiddenEdgesPathChooser<Graph> chooser(g, forbidden_edges);
        size_t current_index = 0;
        for (VertexId v1 : domain_graph_.vertices()) {
            if (++current_index % 100 == 0)
                INFO(current_index << " of " << domain_graph_.size() << " processed.");

            if (domain_graph_.HasStrongEdge(v1) || !domain_graph_.NearContigEnd(v1))
                continue;

            auto reached_vertices = VerticesReachedFrom(g.EdgeEnd(domain_graph_.domain_edges(v1).back()));
            std::set<VertexId> reached_vertices_set(reached_vertices.begin(), reached_vertices.end());

            for (VertexId v2 : domain_graph_.vertices()) {
                if (reached_vertices_set.count(g.EdgeStart(domain_graph_.domain_edges(v2).front())) &&
                    v1 != v2 &&
                    domain_graph_.conjugate(v1) != v2 &&  domain_graph_.GetEdgesBetween(v1, v2).size() == 0 &&
                    !domain_graph_.HasStrongIncomingEdge(v2) && domain_graph_.NearContigStart(v2)) {
                    ConnectWithWeakEdge(v1, v2, chooser);
                }
            }
        }
    }

    class PairComparator {
      public:
        int operator()(const std::pair<int,int> &lhs, const std::pair<int,int> &rhs) const {
            if (lhs.first == rhs.first)
                return lhs.second < rhs.second;
            return lhs.first < rhs.first;
        }
    };

    std::pair<int,int> SearchForSubvector(const path_extend::BidirectionalPath &scaffold, const MappingPath<EdgeId> &domain) const {
        if (domain.size() > scaffold.Size())
            return { -1, -1 };
        scaffold.PrintDEBUG();
        DEBUG(domain.simple_path());
        for (size_t i = 0; i < scaffold.Size() - domain.size() + 1; ++i) {
            for (size_t j = 0; j < domain.size(); ++j) {
                if (scaffold[i + j] != domain[j].first)
                    break;

                if (j == domain.size() - 1) {
                    DEBUG( int(i) << " " << int(i+j) );
                    return { int(i), int(i+j) };
                }
            }
        }
        DEBUG( -1 << " " << -1 );
        return { -1,-1 };
    }

    std::pair<int,int> FindMappingToPath(const path_extend::BidirectionalPath &scaffold, const MappingPath<EdgeId> &domain, std::vector<EdgeId> &edges) const {
        auto res = SearchForSubvector(scaffold, domain);
        if (res.first == -1) {
            return std::make_pair<int,int>(-1,-1);
        }
        size_t start = 0;
        size_t index = 0;
        const auto &g = gp_.get<Graph>();
        for (;index != res.first; ++index) {
            start += g.length(scaffold[index]) + scaffold.GapAt(index).gap;
        }

        start += domain[0].second.mapped_range.start_pos;
        size_t end_index = res.second;
        size_t end = 0;
        for (size_t i = 0; i < end_index; ++i) {
            end += g.length(scaffold[i]) + scaffold.GapAt(i).gap;
        }
        end += domain.back().second.mapped_range.end_pos;
        for (size_t i = index; i <= end_index; ++i) {
            edges.push_back(scaffold[i]);
        }
        DEBUG("Positions " << start << " " << end);
        return std::make_pair(start, end);
    }

    void ConstructStrongEdgesInternal(VertexId current_vertex,
                                      const path_extend::BidirectionalPath &path,
                                      std::map<size_t, std::map<std::pair<int, int>, std::vector<std::pair<VertexId, std::vector<EdgeId>>>, PairComparator>> &mappings_for_path) {
        std::vector<EdgeId> edges;
        std::pair<int, int> coords = FindMappingToPath(path, domain_graph_.mapping_path(current_vertex), edges);
        if (coords.first == -1)
            return;
        mappings_for_path[path.GetId()][coords].push_back(std::make_pair(current_vertex, edges));

        if (coords.second + 5000 > path.Length()) {
            DEBUG("set coord");
            DEBUG(current_vertex);
            DEBUG(domain_graph_.conjugate(current_vertex));
            domain_graph_.SetContigNearEnd(current_vertex);
            DEBUG(domain_graph_.NearContigEnd(current_vertex));
        }
        if (coords.first < 5000) {
            DEBUG("set coord");
            DEBUG(current_vertex);
            DEBUG(domain_graph_.conjugate(current_vertex));
            domain_graph_.SetContigNearEnd(domain_graph_.conjugate(current_vertex));
        }
    }

    size_t GetIndexFromPosition(size_t position, const path_extend::BidirectionalPath &path) const {
        size_t index = 0;
        for (index = 0; index < path.Size(); ++index) {
            size_t current_pos = path.Length() - path.LengthAt(index) + path.graph().length(path[index]);
            if (current_pos > position)
                break;
        }
        return index;
    }

    std::vector<EdgeId> FindEdgesBetweenMappings(int first_mapping_end_coord, int second_mapping_start_coord,
                                                 const path_extend::BidirectionalPath &path) const {
        if (first_mapping_end_coord < 0 || second_mapping_start_coord < 0)
            return {};

        size_t first_mapping_end = GetIndexFromPosition(first_mapping_end_coord, path);
        size_t second_mapping_start = GetIndexFromPosition(second_mapping_start_coord, path);
        if (first_mapping_end > second_mapping_start)
            return {};
    
        std::vector<EdgeId> answer;
        const auto &g = gp_.get<Graph>();
        for (size_t i = first_mapping_end + 1; i < second_mapping_start; ++i) {
            if (answer.size() != 0 && g.EdgeEnd(answer.back()) != g.EdgeStart(path[i])) {
                auto dijkstra = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(g, 500,
                                                                                        30, true /* collect traceback */);
                dijkstra.Run(g.EdgeEnd(answer.back()));
                auto shortest_path = dijkstra.GetShortestPathTo(g.EdgeStart(path[i]));
                DEBUG("Shortest path");
                DEBUG(shortest_path);
                for (auto e : shortest_path)
                    answer.push_back(e);
            }
            answer.push_back(path[i]);
        }

        return answer;
    }

    void ConstructStrongEdges() {
        const auto &g = gp_.get<Graph>();
        path_extend::GraphCoverageMap coverage_map(g, gp_.get<path_extend::PathContainer>("exSPAnder paths"));
        std::map<size_t, std::map<std::pair<int, int>, std::vector<std::pair<VertexId, std::vector<EdgeId>>>, PairComparator>> mappings_for_path;
        std::unordered_map<size_t, path_extend::BidirectionalPath*> from_id_to_path;
        for (VertexId current_vertex : domain_graph_.vertices()) {
            DEBUG("Processing vertex " << domain_graph_.GetVertexName(current_vertex));
            const auto &mapping_path = domain_graph_.mapping_path(current_vertex);
            if (mapping_path.empty())
                continue;

            EdgeId first = mapping_path.front().first;
            auto path_container = coverage_map.GetCoveringPaths(first);
            for (const auto& path_pair : path_container) {
                from_id_to_path[path_pair->GetId()] = path_pair;
                from_id_to_path[path_pair->GetConjPath()->GetId()] = path_pair->GetConjPath();
                ConstructStrongEdgesInternal(current_vertex, *path_pair, mappings_for_path);
            }
        }

        std::set<VertexId> removed_vertices;

        for (auto p : mappings_for_path) {
            DEBUG("Processing path " << p.first);
            if (from_id_to_path[p.first]->IsCanonical())
                continue;

            std::pair<std::pair<int, int>, std::pair<VertexId, std::vector<EdgeId>>> prev(std::make_pair(-1, -1), std::make_pair(VertexId(0), std::vector<EdgeId>()));

            for (auto &external_p : p.second) { // std::pair<std::pair<int, int>, std::vector<std::pair<VertexId, std::vector<EdgeId>>
                for (const auto& imaps : external_p.second) {
                    auto maps = std::make_pair(external_p.first, imaps);
                    DEBUG("Processing mapping " << maps.second.first);
                    DEBUG("Mapping start: " << maps.first.first << ". Mapping end: " << maps.first.second);

                    if (removed_vertices.count(maps.second.first)) {
                        DEBUG("Already deleted");
                        continue;
                    }

                    if (prev.first.first == -1) {
                        prev = maps;
                        continue;
                    }

                    if (prev.first.second > maps.first.first) {
                        DEBUG("Mapping intersects with other, skipping");
                        if (!removed_vertices.count(maps.second.first)) {
                            removed_vertices.insert(domain_graph_.conjugate(maps.second.first));
                            removed_vertices.insert(maps.second.first);
                            if (!domain_graph_.IncomingEdgeCount(maps.second.first) && !domain_graph_.OutgoingEdgeCount(maps.second.first))
                                domain_graph_.DeleteVertex(maps.second.first);
                        }
                        continue;
                    }

                    if (prev.first.second < maps.first.first && maps.first.first - prev.first.second < 20000 &&
                        !removed_vertices.count(maps.second.first) && !removed_vertices.count(prev.second.first) &&
                        domain_graph_.GetEdgesBetween(prev.second.first, maps.second.first).size() == 0) {
                        DEBUG("Connecting " << prev.second << " and " << maps.second);
                        domain_graph_.AddEdge(prev.second.first, maps.second.first, true,
                                              FindEdgesBetweenMappings(prev.first.second, maps.first.first, *from_id_to_path[p.first]), maps.first.first - prev.first.second);
                    }
                    prev = maps;
                }

            }

        }
    }


    //TODO: try some good coverage strategy
    bool IsInsideRepeat(VertexId v) const {
        if (domain_graph_.domain_edges(v).size() > 1)
            return false;

        const auto &g = gp_.get<Graph>();
        EdgeId e = domain_graph_.domain_edges(v).front();
        if (g.IncomingEdgeCount(g.EdgeStart(e)) > 1 ||
            g.OutgoingEdgeCount(g.EdgeEnd(e)) > 1)
            return true;

        return false;
    }

    void ConstructNodes(const nrps::ContigAlnInfo &info) {
        auto mapper = MapperInstance(gp_);
        std::vector<MappingPath<EdgeId>> graph_edges(info.size());
        INFO("Reconstructing alignment paths");
#       pragma omp parallel for
        for (size_t i = 0; i < info.size(); ++i) {
            Sequence sequence(info[i].seq);
            DEBUG(sequence.str());
            graph_edges[i] = mapper->MapSequence(sequence);
        }

        INFO("Creating vertices")
        unsigned id = 1;
        for (size_t i = 0; i < info.size(); ++i) {
            const auto &aln = info[i];
            auto &edges = graph_edges[i];
            if (edges.simple_path().size() == 0)
                continue;

            std::string name = aln.name + "_" + std::to_string(id);
            DEBUG("Adding vertex " << name);

            domain_graph_.AddVertex(name, edges,
                                    edges.front().second.mapped_range.start_pos,
                                    edges.back().second.mapped_range.end_pos,
                                    aln.type, aln.desc);
            id++;
        }
        if (cfg::get().hm->set_copynumber)
            for (VertexId v : domain_graph_.vertices())
                domain_graph_.SetMaxVisited(v, IsInsideRepeat(v) ? 2 : 1);
    }

    GraphPack &gp_;
    nrps::DomainGraph domain_graph_;
    DECL_LOGGER("DomainGraphConstruction");
};

void DomainGraphConstruction::run(GraphPack &gp, const char*) {
    auto res = nrps::DomainMatcher().MatchDomains(gp, cfg::get().hm->hmm_set, cfg::get().output_dir);
    DomainGraphConstructor constructor(gp);
    auto &domain_graph = constructor.ConstructGraph(res);
    const auto &hm = *cfg::get().hm;
    domain_graph.FindDomainOrderings(gp,
                                     hm.component_size_part, 1, hm.start_only_from_tips,
                                     "gene_clusters.fasta", cfg::get().output_dir);
    domain_graph.ExportToDot(fs::append_path(cfg::get().output_dir, "domain_graph.dot"));
}

}
