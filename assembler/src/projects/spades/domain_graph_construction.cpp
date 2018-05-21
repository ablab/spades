#include "domain_graph_construction.hpp"

#include "assembly_graph/paths/path_processor.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"

#include "modules/path_extend/pe_utils.hpp"

#include "domain_graph.hpp"
#include "domain_matcher.hpp"

namespace debruijn_graph {

template<class Graph>
class SetOfForbiddenEdgesPathChooser : public PathProcessor<Graph>::Callback {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph &g_;
    std::set<std::vector<EdgeId>> forbidden_edges_;
    std::vector<EdgeId> answer_path_;
    std::vector<size_t> peaks;

    bool CheckCoverageDiff(const path_extend::BidirectionalPath &path) const {
        double min_coverage = std::numeric_limits<double>::max();
        double max_coverage = std::numeric_limits<double>::min();

        for (size_t i = 0; i < path.Size(); ++i) {
            min_coverage = std::min(min_coverage, path.graph().coverage(path[i]));
            max_coverage = std::max(max_coverage, path.graph().coverage(path[i]));
        }
        return math::ge(50.0, max_coverage/min_coverage);
    }

    bool IsNewPathBetter(path_extend::BidirectionalPath &current, path_extend::BidirectionalPath &candidate) const {
        int current_length = (int)current.Length();
        int candidate_length = (int)candidate.Length();
        int diff_current = std::numeric_limits<int>::max();
        int diff_cand = std::numeric_limits<int>::max();
        for (auto peak : peaks) {
            diff_current = std::min(diff_current, abs((int)peak - current_length));
            diff_cand = std::min(diff_cand, abs((int)peak - candidate_length));
        }
        return diff_cand < diff_current;
    }

public:
    SetOfForbiddenEdgesPathChooser(const Graph &g, const std::set<std::vector<EdgeId>> &forbidden_edges) :
            g_(g), forbidden_edges_(forbidden_edges) {
    }

    void HandleReversedPath(const vector<EdgeId> &reversed_path) override {
        vector<EdgeId> forward_path = reversed_path;
        std::reverse(forward_path.begin(), forward_path.end());
        bool cross_start = false;
        bool cross_end = false;

        for (auto forbidden_path : forbidden_edges_) {
            if (forbidden_path.size() != 0) {
                if (find(std::begin(forward_path), std::end(forward_path), forbidden_path.front()) != forward_path.end()) {
                    cross_start = true;
                }
                if (find(std::begin(forward_path), std::end(forward_path), forbidden_path.back()) != forward_path.end()) {
                    cross_end = true;
                }
                if (cross_end && cross_start) {
                    return;
                }
            }
        }
        if (answer_path_.empty()) {
            if (!CheckCoverageDiff(path_extend::BidirectionalPath(g_, forward_path))) {
               return;
            }
            answer_path_ = forward_path;
            return;
        }

        if (!CheckCoverageDiff(path_extend::BidirectionalPath(g_,forward_path))) {
            return;
        }

        path_extend::BidirectionalPath current(g_, answer_path_);
        path_extend::BidirectionalPath candidate(g_, forward_path);
        if (IsNewPathBetter(current, candidate)) {
            answer_path_ = forward_path;
        }
    }

    void reset() {
        answer_path_.clear();
    }

    const vector<EdgeId>& answer() {
        return answer_path_;
    }

    void set_peaks(size_t first, size_t second) {
        peaks.clear();
        peaks.push_back(first);
        peaks.push_back(second);
    }
};

class DomainGraphConstructor {
public:
    DomainGraphConstructor(conj_graph_pack &gp)
            : gp_(gp) {}

    nrps::DomainGraph& ConstructGraph(const nrps::ContigAlnInfo &info) {
        gp_.EnsureIndex();
        INFO("A-domain graph construction started");
        ConstructNodes(info);
        ConstructStrongEdges();
        ConstructWeakEdges();
        INFO("A-domain graph construction ended");
        return graph;
    }

private:

    void ConnectWithWeakEdge(std::shared_ptr<nrps::Vertex> v1, std::shared_ptr<nrps::Vertex> v2,
                             SetOfForbiddenEdgesPathChooser<Graph> &chooser) {
        DEBUG("Trying to connect " << v1->name_ << " and " << v2->name_ << " with weak edge");
        int last_mapping = (int)gp_.g.length(mappings[v1->name_].back().first) - mappings[v1->name_].end_pos();
        int first_mapping = (int)mappings[v2->name_].start_pos();

        if (5000 < last_mapping + first_mapping) {
            return;
        }
        int min_len = 0;
        chooser.set_peaks(std::max(0, 500 - last_mapping - first_mapping), std::max(0, 1000 - last_mapping - first_mapping));
        if (gp_.g.EdgeEnd(v1->domain_edges_in_row_.back()) != gp_.g.EdgeStart(v2->domain_edges_in_row_[0])) {
            DEBUG("Trying to find paths from " << gp_.g.EdgeEnd(v1->domain_edges_in_row_.back()) << " to " << gp_.g.EdgeStart(v2->domain_edges_in_row_[0]));
            ProcessPaths(gp_.g, min_len, 4000 - last_mapping - first_mapping, gp_.g.EdgeEnd(v1->domain_edges_in_row_.back()), gp_.g.EdgeStart(v2->domain_edges_in_row_[0]), chooser);
            if(!chooser.answer().empty()) {
                DEBUG("Path was found");
                path_extend::BidirectionalPath p(gp_.g, chooser.answer());
                DEBUG("Start vertex: " << gp_.g.EdgeStart(p.Front()).int_id());
                DEBUG("End vertex: " << gp_.g.EdgeEnd(p.Back()).int_id());
                DEBUG("Path:");
                p.PrintDEBUG();
                graph.addEdge(v1, v2, false, p.Length() + last_mapping + first_mapping, chooser.answer());
            }
            else {
                DEBUG("Path was not found");
            }
        } else {
            graph.addEdge(v1, v2, false, last_mapping + first_mapping, vector<EdgeId>());
        }
        chooser.reset();
    }

    void ConstructWeakEdges() {
        std::set<std::vector<EdgeId>> forbidden_edges;
        for (const auto &mapping : mappings)
            forbidden_edges.insert(mapping.second.simple_path());

        SetOfForbiddenEdgesPathChooser<Graph> chooser(gp_.g, forbidden_edges);
        for (auto v1 : graph.getNodeSet()) {
            if (!graph.HasStrongEdge(v1) && (v1->near_to_the_end_of_contig_)) {
                for (auto v2 : graph.getNodeSet()) {
                    if (v1 != v2 && v1->rc_ != v2 && !graph.HasStrongIncomingEdge(v2) && v2->near_to_the_start_of_contig_) {
                        ConnectWithWeakEdge(v1, v2, chooser);
                    }
                }
            }
        }
    }

    class PairComparator {
      public:
        int operator()(const std::pair<int,int> &lhs, const std::pair<int,int> &rhs) const {
            return lhs.first < rhs.first;
        }
    };

    std::pair<int,int> SearchForSubvector(path_extend::BidirectionalPath *scaffold, MappingPath<EdgeId> &domain) {
        if (domain.size() > scaffold->Size())
            return { -1, -1 };

        for (size_t i = 0; i < scaffold->Size() - domain.size() + 1; ++i) {
            for (size_t j = 0; j < domain.size(); ++j) {
                if ((*scaffold)[i + j] != domain[j].first)
                    break;

                if (j == domain.size() - 1)
                    return { int(i), int(i+j) };
            }
        }

        return { -1,-1 };
    }

    std::pair<int,int> FindMappingToPath(path_extend::BidirectionalPath *scaffold, MappingPath<EdgeId> &domain, std::vector<EdgeId> &edges) {
        auto res = SearchForSubvector(scaffold, domain);
        if (res.first == -1) {
            return std::make_pair<int,int>(-1,-1);
        }
        int start = 0;
        size_t index = 0;
        for (;index != res.first; ++index) {
            start += gp_.g.length((*scaffold)[index]);
        }

        start += domain[0].second.mapped_range.start_pos;
        size_t end_index = index + domain.size() - 1;
        int end = 0;
        for (size_t i = 0; i < end_index; ++i) {
            end += gp_.g.length((*scaffold)[i]);
        }
        end += domain[domain.size() - 1].second.mapped_range.end_pos;
        for (int i = index; i <= end_index; ++i) {
            edges.push_back((*scaffold)[i]);
        }
        return std::make_pair(start, end);
    }

    void ConstructStrongEdgesInternal(std::pair<const std::string, MappingPath<EdgeId>> &mappings, path_extend::BidirectionalPath *path,
            std::map<size_t, std::map<std::pair<int, int>, std::pair<std::string, std::vector<EdgeId>>, PairComparator>> &mappings_for_path) {
        std::vector<EdgeId> edges;
        std::pair<int, int> coords = FindMappingToPath(path, mappings.second, edges);
        if (coords.first == -1) {
            return;
        }
        mappings_for_path[path->GetId()][coords] = std::make_pair(mappings.first, edges);

        if (coords.second + 5000 < path->Length()) {
           auto v = graph.getNode(mappings.first);
           v->near_to_the_end_of_contig_ = false;
        }
        if (coords.first > 5000) {
            auto v = graph.getNode(mappings.first);
            v->near_to_the_start_of_contig_ = false;
        }
    }

    size_t GetIndexFromPosition(size_t position, path_extend::BidirectionalPath *path) {
        size_t index = 0;
        for (index = 0; index < path->Size(); ++index) {
            size_t current_pos = path->LengthAt(index);
            if (current_pos > position) {
                break;
            }
        }
        return index;
    }

    vector<EdgeId> FindEdgesBetweenMappings(int first_mapping_end_coord, int second_mapping_start_coord, path_extend::BidirectionalPath *path) {
        vector<EdgeId> answer;
        if (first_mapping_end_coord < 0 || second_mapping_start_coord < 0) {
            return answer;
        }
        size_t first_mapping_end = GetIndexFromPosition(first_mapping_end_coord, path);
        size_t second_mapping_start = GetIndexFromPosition(second_mapping_start_coord, path);

        if (first_mapping_end > second_mapping_start) {
            return answer;
        }

        for (size_t i = first_mapping_end + 1; i < second_mapping_start; ++i) {
            answer.push_back((*path)[i]);
        }
        return answer;
    }

    void ConstructStrongEdges() {
        path_extend::GraphCoverageMap coverage_map(gp_.g, gp_.contig_paths);
        std::map<size_t, std::map<std::pair<int, int>, std::pair<std::string, std::vector<EdgeId>>, PairComparator>> mappings_for_path;
        std::map<size_t, path_extend::BidirectionalPath*> from_id_to_path;
        for (auto domain : mappings) {
            DEBUG("Processing mapping " << domain.first);
            auto mapping_path = domain.second;
            if (!mapping_path.empty()) {
                EdgeId first = mapping_path.front().first;
                auto path_container = coverage_map.GetCoveringPaths(first);
                for (auto path_pair : path_container) {
                    from_id_to_path[path_pair->GetId()] = path_pair;
                    from_id_to_path[path_pair->GetConjPath()->GetId()] = path_pair->GetConjPath();
                    ConstructStrongEdgesInternal(domain, path_pair, mappings_for_path);
                    ConstructStrongEdgesInternal(domain, path_pair->GetConjPath(), mappings_for_path);
                }
            }
        }
        std::set<std::string> removed_vertices;
        for (auto p : mappings_for_path) {
            DEBUG("Processing path " << p.first);
            std::pair<std::pair<int, int>, std::pair<std::string, vector<EdgeId>>> prev(make_pair(-1, -1), make_pair("", std::vector<EdgeId>()));
            for (const auto& maps : p.second) {
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
                        DEBUG("Inserting " << graph.getVertex(maps.second.first)->rc_->name_);
                        DEBUG("Inserting " << maps.second.first);
                        removed_vertices.insert(graph.getVertex(maps.second.first)->rc_->name_);
                        removed_vertices.insert(maps.second.first);
                        graph.removeVertex(maps.second.first);
                    }
                    continue;
                }

                if (prev.first.second < maps.first.first && maps.first.first - prev.first.second < 20000 && !removed_vertices.count(maps.second.first) && !removed_vertices.count(prev.second.first)) {
                    DEBUG("Connecting " << prev.second << " and " << maps.second);
                    graph.addEdge(prev.second.first, maps.second.first, true, maps.first.first - prev.first.second, FindEdgesBetweenMappings(prev.first.second, maps.first.first, from_id_to_path[p.first]));
                }
                prev = maps;
            }
        }
    }

    //TODO: try some good coverage strategy
    bool IsInsideRepeat(std::shared_ptr<nrps::Vertex> v) {
        if (v->domain_edges_in_row_.size() > 1)
            return false;

        EdgeId e = v->domain_edges_in_row_[0];
        if (gp_.g.IncomingEdgeCount(gp_.g.EdgeStart(e)) > 1 ||
            gp_.g.OutgoingEdgeCount(gp_.g.EdgeEnd(e)) > 1)
            return true;

        return false;
    }

    void ConstructNodes(const nrps::ContigAlnInfo &info) {
        auto mapper = MapperInstance(gp_);

        unsigned id = 1;
        for (const auto &aln : info) {
            Sequence sequence = Sequence(aln.seq);
            std::string name = aln.name + "_" + std::to_string(id), name_rc = name + "rc";

            auto edges = mapper->MapSequence(sequence);
            if (edges.simple_path().size() == 0)
                continue;

            DEBUG("Adding vertex " << name);
            graph.addVertex(name,
                            edges.simple_path(),
                            edges.front().second.mapped_range.start_pos, edges.back().second.mapped_range.end_pos,
                            aln.type);
            auto rc_edges = mapper->MapSequence(!sequence);
            DEBUG("Adding vertex " << name_rc);
            graph.addVertex(name_rc,
                            rc_edges.simple_path(),
                            rc_edges.front().second.mapped_range.start_pos, rc_edges.back().second.mapped_range.end_pos,
                            aln.type);
            mappings[name] = edges;
            mappings[name_rc] = rc_edges;
            graph.makeRC(name, name_rc);
            id++;
        }

        for (auto v : graph.getVertexSet()) {
            v->max_visited_ = IsInsideRepeat(v) ? 2 : 1;
        }
    }

    conj_graph_pack &gp_;
    nrps::DomainGraph graph;
    std::map<std::string, MappingPath<EdgeId>> mappings;
    DECL_LOGGER("AGraph");
};

void DomainGraphConstruction::run(conj_graph_pack &gp, const char*) {
    auto res = nrps::DomainMatcher().MatchDomains(gp, cfg::get().hmm_set, cfg::get().output_dir);

    DomainGraphConstructor constructor(gp);
    auto &graph = constructor.ConstructGraph(res);

    graph.FindDomainOrderings(gp, "gene_clusters.fasta", cfg::get().output_dir);
    graph.ExportToDot(cfg::get().output_dir + "/domain_graph.dot");
    INFO("Export paths");
    graph.ExportPaths(gp, cfg::get().output_dir);
}

}
