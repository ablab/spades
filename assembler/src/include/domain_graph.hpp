//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/graph_support/contig_output.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

#include "io/reads/osequencestream.hpp"
#include "pipeline/graph_pack.hpp"

#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace nrps {
typedef debruijn_graph::EdgeId EdgeId;
struct Edge;

struct Vertex {
public:
    std::string name_;
    std::set<Edge *> edges_;

    std::set<debruijn_graph::EdgeId> unique_domain_edges_;
    std::vector<debruijn_graph::EdgeId> domain_edges_in_row_;
    size_t start_coord_;
    size_t end_coord_;
    std::shared_ptr<Vertex> rc_;

    std::string v_type_;
    bool near_to_the_end_of_contig_;
    bool near_to_the_start_of_contig_;

    bool visited_;
    size_t visited_times_;
    size_t max_visited_;

    /*
     * Constructs a vertex with the given name.
     */
    Vertex(const std::set<EdgeId> &unique_a_edges,
           const std::vector<EdgeId> &domain_edges, size_t start_coord,
           size_t end_coord, std::string type, const std::string &name = "")
            : name_(name), unique_domain_edges_(unique_a_edges),
              domain_edges_in_row_(domain_edges), start_coord_(start_coord),
              end_coord_(end_coord), v_type_(type), near_to_the_end_of_contig_(true),
              near_to_the_start_of_contig_(true), visited_(false), visited_times_(0),
              max_visited_(1) {}

    Vertex(const std::string &name = "") : name_(name) {}
    ~Vertex() {}

    Vertex(const Vertex &other) = default;
    Vertex &operator=(const Vertex &other) = default;
    Vertex &operator=(Vertex &&other) = default;
};

struct Edge {
public:
    std::shared_ptr<Vertex> start_;
    std::shared_ptr<Vertex> end_;
    bool visited_;
    bool strong_;
    size_t length_;
    std::vector<EdgeId> edges_;

    Edge(std::shared_ptr<Vertex> start = nullptr,
         std::shared_ptr<Vertex> end = nullptr,
         const std::vector<EdgeId> &edges_on_graph = std::vector<EdgeId>(),
         bool isstrong = false)
            : start_(start), end_(end), visited_(false), strong_(isstrong),
              length_(0), edges_(edges_on_graph) {}
    ~Edge() {}

    Edge &operator=(const Edge &other) = default;
    Edge &operator=(Edge &&other) = default;
};

class DomainGraph {
public:
    DomainGraph() { iteration_number_ = 0; }

    std::shared_ptr<Vertex> getNode(const std::string &name) const {
        if (node_map_.find(name) != node_map_.end())
            return node_map_.at(name);

        return nullptr;
    }

    Edge *getArc(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2) const {
        for (Edge *edge : this->getEdgeSet(v1)) {
            if (edge->end_ == v2) {
                return edge;
            }
        }

        return nullptr;
    }

    Edge *getEdge(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2) const {
        return this->getArc(v1, v2);
    }

    bool HasStrongEdge(std::shared_ptr<nrps::Vertex> v) const {
        return StrongEdgeCount(v) > 0;
    }

    bool HasStrongIncomingEdge(std::shared_ptr<nrps::Vertex> v) const {
        auto rc = v->rc_;
        return StrongEdgeCount(rc) > 0;
    }

    size_t StrongEdgeCount(std::shared_ptr<nrps::Vertex> v) const {
        return std::count_if(v->edges_.begin(), v->edges_.end(),
                             [](Edge *e) { return e->strong_; });
    }

    size_t WeakEdgeCount(std::shared_ptr<nrps::Vertex> v) const {
        return std::count_if(v->edges_.begin(), v->edges_.end(),
                             [](Edge *e) { return !e->strong_; });
    }

    Edge *addEdge(const std::string &v1, const std::string &v2, bool strong,
                  size_t length = 0,
                  const vector<EdgeId> &edges = vector<EdgeId>()) {
        return this->addEdge(getVertex(v1), getVertex(v2), strong, length, edges);
    }

    Edge *addEdge(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2,
                  bool strong, size_t length = 0,
                  const vector<EdgeId> &edges = vector<EdgeId>()) {
        Edge *e = new Edge(v1, v2, edges, strong);
        e->length_ = length;
        return addEdge(e);
    }

    void removeVertex(const std::string &vertex_name) {
        for (auto v : nodes_) {
            if (v->name_ == vertex_name) {
                nodes_.erase(v);
                break;
            }
        }
        auto it = node_map_.find(vertex_name);
        for (auto edge : it->second->edges_) {
            delete edge;
        }
        node_map_.erase(it);
    }

    bool isExistingNode(std::shared_ptr<Vertex> v) const {
        return this->getNodeSet().find(v) != this->getNodeSet().end();
    }

    Edge *addEdge(Edge *e) {
        Edge *result = this->addEdgeInternal(e);
        return result;
    }

    std::shared_ptr<Vertex>
    addVertex(const std::string &name,
              const std::vector<debruijn_graph::EdgeId> &a_edges,
              size_t start_coord, size_t end_coord, std::string type) {
        std::set<EdgeId> unique_edges;
        for (auto edge : a_edges) {
            if (a_edges_map_.find(edge) == a_edges_map_.end()) {
                unique_edges.insert(edge);
            }
            if (std::find(a_edges_map_[edge].begin(), a_edges_map_[edge].end(),
                          name) == a_edges_map_[edge].end()) {
                a_edges_map_[edge].push_back(name);
            }
        }
        return this->addVertexInternal(name, unique_edges, a_edges, start_coord,
                                       end_coord, type);
  }

    std::shared_ptr<Vertex> addVertex(std::shared_ptr<Vertex> v) {
        return this->addVertexInternal(v);
    }

    const std::set<Edge *> &getEdgeSet(std::shared_ptr<Vertex> v) const {
        return v->edges_;
    }

    std::shared_ptr<Vertex> getVertex(const std::string &name) const {
        return node_map_.at(name);
    }

    const std::set<std::shared_ptr<Vertex>> &getVertexSet() const {
        return nodes_;
    }

    std::set<std::shared_ptr<Vertex>> const &getNodeSet() const { return nodes_; }

    void ExportToDot(const std::string &output_path) const {
        std::ofstream out(output_path);
        out << "digraph domain_graph {" << std::endl;
        for (auto v : nodes_) {
            out << "\"" << v->name_ << "\""
                << " [label=\"" << v->v_type_ << " " << v->name_ << " "
                << v->max_visited_ << "\"];" << std::endl;
        }

        for (auto e : arcs_) {
            out << "\"" << e->start_->name_ << "\""
                << " -> "
                << "\"" << e->end_->name_ << "\"";
            out << " [label=";
            out << e->length_;
            if (e->strong_) {
                out << " style=bold];" << std::endl;
            } else {
                out << " style=dotted];" << std::endl;
            }
        }
        out << "}";
    }

    void makeRC(const std::string &first_vertex,
                const std::string &second_vertex) {
        node_map_[first_vertex]->rc_ = node_map_[second_vertex];
        node_map_[second_vertex]->rc_ = node_map_[first_vertex];
    }

    void FindBasicStatistic(std::ofstream &stat_stream) {
        stat_stream << "A - Adenylation domain" << std::endl;
        stat_stream << "AT - Acyltransferase domain" << std::endl;
        stat_stream << "C - Condensation domain" << std::endl;
        stat_stream << "KR - Keto-reductase domain" << std::endl;
        stat_stream << "KS - Keto-synthase domain" << std::endl;
        stat_stream << "TE - Termination domain" << std::endl;
        stat_stream << std::endl;

        std::map<std::string, int> domain_count;

        for (auto v : nodes_)
            domain_count[v->v_type_]++;

        for (auto domain_type : domain_count) {
            stat_stream << "# " << domain_type.first << "-domains - "
                        << domain_type.second << std::endl;
        }
        stat_stream << std::endl;
    }

    void FindDomainOrderings(debruijn_graph::conj_graph_pack &gp,
                             std::string output_filename) {
        std::ofstream stat_stream(cfg::get().output_dir + "/bgc_statistics.txt");
        FindBasicStatistic(stat_stream);
        path_extend::ContigWriter writer(
            gp.g, path_extend::MakeContigNameGenerator(cfg::get().mode, gp));
        std::set<std::shared_ptr<Vertex>> nodes_with_incoming;
        for (auto e : arcs_) {
            nodes_with_incoming.insert(e->end_);
        }

        std::set<std::shared_ptr<Vertex>> start_nodes;
        for (auto v : nodes_) {
            if (nodes_with_incoming.find(v) == nodes_with_incoming.end() &&
                v->edges_.size()) {
                start_nodes.insert(v);
            }
        }
        io::osequencestream_bgc oss(output_filename);
        path_extend::ScaffoldSequenceMaker seq_maker(gp.g);
        std::vector<std::vector<std::shared_ptr<Vertex>>> answer;
        int ordering_id = 1;
        int component_id = 1;

        for (auto v : start_nodes) {
            if (!v->visited_) {
                stat_stream << "BGC subgraph " << component_id << std::endl;
                FindAllPossibleArrangements(gp.g, v, answer, stat_stream);
                ordering_id = 1;
                for (auto vec : answer) {
                    OutputStatArrangement(gp.g, vec, ordering_id, stat_stream);
                    path_extend::BidirectionalPath *p =
                            new path_extend::BidirectionalPath(gp.g);
                    std::string outputstring = PathToSequence(p, vec);
                    path_extend::BidirectionalPath *conjugate =
                            new path_extend::BidirectionalPath(p->Conjugate());
                    contig_paths_.AddPair(p, conjugate);
                    OutputComponent(gp, p, component_id, ordering_id);
                    outputstring = seq_maker.MakeSequence(*p);
                    oss.SetCluster(component_id, ordering_id);
                    oss << outputstring;
                    ordering_id++;
                }
                answer.clear();
                component_id++;
            }
        }
    }

    void ExportPaths(debruijn_graph::conj_graph_pack &gp,
                     std::string output_dir) {
        auto name_generator =
                path_extend::MakeContigNameGenerator(cfg::get().mode, gp);
        path_extend::ContigWriter writer(gp.g, name_generator);
        writer.OutputPaths(contig_paths_, output_dir + "/orderings2");
    }

private:
    std::shared_ptr<Vertex>
    addVertexInternal(const std::string &name,
                      const std::set<debruijn_graph::EdgeId> &unique_edges,
                      const std::vector<debruijn_graph::EdgeId> &edges,
                      size_t start_coord, size_t end_coord, std::string type) {
        std::shared_ptr<Vertex> node = std::make_shared<Vertex>(
            unique_edges, edges, start_coord, end_coord, type);
        node->name_ = name;
        return addVertexInternal(node);
    }

    std::shared_ptr<Vertex> addVertexInternal(std::shared_ptr<Vertex> node) {
        nodes_.insert(node);
        DEBUG("Node " << node->name_ << " is inserted into the map");
        node_map_[node->name_] = node;
        return node;
    }

    Edge *addEdgeInternal(Edge *arc) {
        if (!isExistingNode(arc->start_))
            addVertex(arc->start_);
        if (!isExistingNode(arc->end_))
            addVertex(arc->end_);
        arc->start_->edges_.insert(arc);
        arcs_.insert(arc);
        return arc;
    }

    string ToString(int integer) {
        stringstream ss;
        ss << integer;
        return ss.str();
    }

    void OutputComponent(debruijn_graph::conj_graph_pack &gp,
                         path_extend::BidirectionalPath *p, int component_id,
                         int ordering_id) {
        auto edges = CollectEdges(p);
        auto comp = GraphComponent<Graph>::FromEdges(gp.g, edges.begin(),
                                                     edges.end(), true);
        std::ofstream os(cfg::get().output_dir + "/bgc_in_gfa/" +
                         ToString(component_id) + "_" + ToString(ordering_id) +
                         ".gfa");
        path_extend::GFAWriter<Graph> writer(gp.g, os);
        writer.Write(comp);
    }

    std::set<EdgeId> CollectEdges(path_extend::BidirectionalPath *p) {
        std::set<EdgeId> edge_set;
        for (size_t i = 0; i < p->Size(); ++i) {
            edge_set.insert(p->At(i));
        }
        return edge_set;
    }

    std::string PathToSequence(path_extend::BidirectionalPath *p,
                               std::vector<std::shared_ptr<Vertex>> &answer) {
        std::stringstream ss;
        DEBUG("Translating " << p->GetId() << " to sequnece");
        for (size_t i = 0; i < answer.size(); ++i) {

            auto v = answer[i];
            DEBUG("Translating vertex " << v->name_);

            for (auto e : v->domain_edges_in_row_) {
                if (p->Size() == 0 || p->Back() != e) {
                    int gap = 0;
                    if (p->Size() != 0 &&
                        p->graph().EdgeEnd(p->Back()) != p->graph().EdgeStart(e)) {
                        gap = 100;
                    }
                    p->PushBack(e, path_extend::Gap(gap));
                }
            }
            if (i != answer.size() - 1) {
                Edge *next_edge = this->getEdge(v, answer[i + 1]);
                DEBUG("Translating edge between " <<
                      v->name_ << " and " << answer[i + 1]->name_);
                DEBUG("Is edge strong " << next_edge->strong_);
                if (next_edge == nullptr) {
                    DEBUG("Something strange!");
                    continue;
                }
                for (auto e : next_edge->edges_) {
                    if (p->Size() == 0 || p->Back() != e) {
                        int gap = 0;
                        if (p->Size() != 0 &&
                            p->graph().EdgeEnd(p->Back()) != p->graph().EdgeStart(e)) {
                            gap = 100;
                        }
                        p->PushBack(e, path_extend::Gap(gap));
                    }
                }
            }
        }
        return ss.str();
    }

    void
    OutputStatArrangement(const Graph &g,
                          std::vector<std::shared_ptr<Vertex>> single_candidate,
                          int id, std::ofstream &stat_file) {
        stat_file << "BGC candidate " << id << std::endl;
        std::string delimeter = "";
        bool is_nrps = false;
        bool is_pks = false;
        for (auto v : single_candidate) {
            stat_file << delimeter;
            delimeter = "-";
            if (v->v_type_ == "AMP") {
                stat_file << "A";
                is_nrps = true;
            } else if (v->v_type_ == "CStart") {
                stat_file << "C";
                is_nrps = true;
            } else if (v->v_type_ == "AT") {
                stat_file << "AT";
                is_pks = true;
            } else if (v->v_type_ == "TE") {
                stat_file << "TE";
            } else if (v->v_type_ == "KR") {
                stat_file << "KR";
                is_pks = true;
            } else if (v->v_type_ == "KS") {
                stat_file << "KS";
                is_pks = true;
            } else {
                stat_file << v->v_type_;
            }
        }
        stat_file << std::endl;
        stat_file << "Predicted type: ";
        if (is_nrps && is_pks) {
            stat_file << "NRPS/PKS";
        } else if (is_nrps) {
            stat_file << "NRPS";
        } else if (is_pks) {
            stat_file << "PKS";
        } else {
            stat_file << "None";
        }
        stat_file << std::endl;
        stat_file << "Domain cordinates:" << std::endl;
        size_t current_coord = 0;
        path_extend::BidirectionalPath *p = new path_extend::BidirectionalPath(g);
        for (size_t i = 0; i < single_candidate.size(); ++i) {
            auto v = single_candidate[i];
            DEBUG("Translating vertex " << v->name_);
            for (auto e : v->domain_edges_in_row_) {
                if (p->Size() == 0 || p->Back() != e) {
                    int gap = 0;
                    if (p->Size() != 0 &&
                        p->graph().EdgeEnd(p->Back()) != p->graph().EdgeStart(e)) {
                        gap = 100;
                    }
                    p->PushBack(e, path_extend::Gap(gap));
                    current_coord += g.length(e) + gap;
                }
            }

            size_t sum = 0;
            for (auto e2 : v->domain_edges_in_row_) {
                sum += g.length(e2);
            }

            stat_file << v->start_coord_ + current_coord - sum << " ";

            stat_file << v->end_coord_ + current_coord - g.length(p->Back())
                << std::endl;

            if (i != single_candidate.size() - 1) {
                Edge *next_edge = this->getEdge(v, single_candidate[i + 1]);
                if (next_edge == nullptr) {
                    continue;
                }
                for (auto e : next_edge->edges_) {
                    if (p->Size() == 0 || p->Back() != e) {
                        int gap = 0;
                        if (p->Size() != 0 &&
                            p->graph().EdgeEnd(p->Back()) != p->graph().EdgeStart(e)) {
                            gap = 100;
                        }
                        p->PushBack(e, path_extend::Gap(gap));
                        current_coord += g.length(e) + gap;
                    }
                }
            }
        }
    }

    void OutputStat(std::set<std::shared_ptr<Vertex>> &preliminary_visited,
                    std::ofstream &stat_file) {
        stat_file << "# domains in the component - ";
        stat_file << preliminary_visited.size() << std::endl;
        int strong_edge_count =
                std::accumulate(preliminary_visited.begin(), preliminary_visited.end(),
                                0, [&](int prev, std::shared_ptr<Vertex> v) {
                                    return prev + StrongEdgeCount(v);
                                });
        int weak_edge_count =
                std::accumulate(preliminary_visited.begin(), preliminary_visited.end(),
                                0, [&](int prev, std::shared_ptr<Vertex> v) {
                                    return prev + (WeakEdgeCount(v));
                                });

        stat_file << "# strong/weak edges in the component - ";
        stat_file << strong_edge_count << "/" << weak_edge_count << std::endl;

        if (preliminary_visited.size() <= 5) {
            //                stat_file << "Warning: small number of domains, possible
            //                fragmentation" << std::endl;
        }

        if (strong_edge_count < (int)preliminary_visited.size() / 2) {
            //                stat_file << "Warning: small number of strong edges,
            //                unlikely to recover complete BGC" << std::endl;
        }
    }

    size_t GetMaxVisited(const Graph &g, std::shared_ptr<Vertex> v,
                         double base_coverage) {
        double low_coverage = std::numeric_limits<double>::max();
        for (auto e : v->domain_edges_in_row_)
            low_coverage = std::min(low_coverage, g.coverage(e));

        return size_t(round(low_coverage / base_coverage));
  }

    void SetCopynumber(const Graph &g,
                       std::set<std::shared_ptr<Vertex>> &preliminary_visited) {
        double base_coverage = std::numeric_limits<double>::max();
        for (auto v : preliminary_visited) {
            for (auto e : v->domain_edges_in_row_) {
                if (g.length(e) > 500) {
                    base_coverage = std::min(base_coverage, g.coverage(e));
                }
            }
        }

        if (math::eq(base_coverage, std::numeric_limits<double>::max()))
            return;

        for (auto v : preliminary_visited) {
            size_t value = std::max(size_t(1), GetMaxVisited(g, v, base_coverage));
            if (v->max_visited_ != value) {
                DEBUG(v->name_ << " copynumber has changed from " << v->max_visited_
                      << " to " << value);
      }
            v->max_visited_ = value;
        }
    }

    void FindAllPossibleArrangements(const Graph &g, std::shared_ptr<Vertex> v,
                                     std::vector<std::vector<std::shared_ptr<Vertex>>> &answer,
                                     std::ofstream &stat_file) {
        std::set<std::shared_ptr<Vertex>> preliminary_visited;
        preliminary_visited.insert(v);
        PrelimDFS(v, preliminary_visited);
        size_t component_size = preliminary_visited.size();
        if (component_size == 1)
            return;

        SetCopynumber(g, preliminary_visited);
        OutputStat(preliminary_visited, stat_file);
        std::vector<std::shared_ptr<Vertex>> current;
        for (auto v : preliminary_visited) {
            v->visited_ = true;
            v->rc_->visited_ = true;
        }
        preliminary_visited.clear();
        iteration_number_ = 0;
        FinalDFS(v, current, preliminary_visited, answer, component_size);
    }

    void PrelimDFS(std::shared_ptr<Vertex> v,
                   std::set<std::shared_ptr<Vertex>> &preliminary_visited) {
        for (auto e : v->edges_) {
            std::shared_ptr<Vertex> to = e->end_;
            if (preliminary_visited.find(to) == preliminary_visited.end()) {
                preliminary_visited.insert(to);
                PrelimDFS(to, preliminary_visited);
            }
        }
    }

    void FinalDFS(std::shared_ptr<Vertex> v,
                  std::vector<std::shared_ptr<Vertex>> &current,
                  std::set<std::shared_ptr<Vertex>> preliminary_visited,
                  std::vector<std::vector<std::shared_ptr<Vertex>>> &answer,
                  size_t component_size) {
        iteration_number_++;

        current.push_back(v);
        v->visited_times_++;
        if (answer.size() > 50 || iteration_number_ > 10000000) {
            return;
        }

        bool was_extended = false;
        if (this->HasStrongEdge(v)) {
            for (auto arc : v->edges_) {
                if (arc->strong_) {
                    std::shared_ptr<Vertex> to = arc->end_;
                    if (preliminary_visited.find(to) != preliminary_visited.end()) {
                        continue;
                    }
                    was_extended = true;
                    FinalDFS(to, current, preliminary_visited, answer, component_size);
                }
            }
        } else {
            for (auto arc : v->edges_) {
                std::shared_ptr<Vertex> to = arc->end_;
                if (to->visited_times_ >= to->max_visited_) {
                    continue;
                }
                if (preliminary_visited.find(to) != preliminary_visited.end()) {
                    continue;
                }
                was_extended = true;
                FinalDFS(to, current, preliminary_visited, answer, component_size);
            }
        }
        if (current.size() >= component_size && was_extended == false) {
            answer.push_back(current);
        }
        v->visited_times_--;
        current.pop_back();
    }

    std::set<std::shared_ptr<Vertex>> nodes_; /* The set of nodes in the graph */
    std::set<Edge *> arcs_;                   /* The set of arcs in the graph  */
    std::map<std::string, std::shared_ptr<Vertex>> node_map_; /* A map from names to nodes     */
    std::map<EdgeId, std::vector<std::string>> a_edges_map_; /* A map from edges to nodenames */
    path_extend::PathContainer contig_paths_;
    size_t iteration_number_;

    DECL_LOGGER("AGraph2");
};

}
