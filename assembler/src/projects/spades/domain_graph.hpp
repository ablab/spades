//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/graph_support/contig_output.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

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

struct cmpshared_ptrs {
    bool operator()(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2) const {
        return v1->name_ < v2->name_;
    }
};

public:
    DomainGraph() { iteration_number_ = 0; }

    std::shared_ptr<Vertex> getNode(const std::string &name) const {
        auto entry = node_map_.find(name);
        if (entry != node_map_.end())
            return entry->second;

        return nullptr;
    }

    Edge *getArc(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2) const {
        for (Edge *edge : getEdgeSet(v1)) {
            if (edge->end_ == v2) {
                return edge;
            }
        }

        return nullptr;
    }

    Edge *getEdge(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2) const {
        return getArc(v1, v2);
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
                  const std::vector<EdgeId> &edges = std::vector<EdgeId>()) {
        return addEdge(getVertex(v1), getVertex(v2), strong, length, edges);
    }

    Edge *addEdge(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2,
                  bool strong, size_t length = 0,
                  const std::vector<EdgeId> &edges = std::vector<EdgeId>()) {
        Edge *e = new Edge(v1, v2, edges, strong);
        e->length_ = length;
        return addEdge(e);
    }

    void removeVertex(const std::string &vertex_name) {
        auto vertex = getVertex(vertex_name);
        auto vertex_rc = getVertex(vertex_name)->rc_;
        for (auto v : nodes_) {
            if (v->name_ == vertex_name || v->name_ == vertex_rc->name_) {
                nodes_.erase(v);
            }
        }

        auto it = node_map_.find(vertex_name);
        for (auto edge : it->second->edges_) {
            removeEdge(edge);
        }

        auto it2 = node_map_.find(vertex_rc->name_);
        for (auto edge : it2->second->edges_) {
            removeEdge(edge);
        }
        node_map_.erase(it);
        node_map_.erase(it2);
    }

    void removeEdge(Edge *edge) {
        for (auto e : arcs_) {
            if (e->start_ == edge->start_ && e->end_ == edge->end_ ) {
                auto v_rc = edge->end_->rc_;
                for (auto e2 : v_rc->edges_) {
                    if (e2->end_ == e->start_->rc_) {
                        v_rc->edges_.erase(e2);
                        break;
                    }
                }
                v_rc->edges_.erase(e);
                arcs_.erase(e);
            }
        }
        arcs_.erase(edge);
        delete edge;
    }

    bool isExistingNode(std::shared_ptr<Vertex> v) const {
        return nodes_.count(v);
    }

    std::shared_ptr<Vertex> addVertex(const std::string &name,
                                      const std::vector<debruijn_graph::EdgeId> &a_edges,
                                      size_t start_coord, size_t end_coord, std::string type) {
        std::set<EdgeId> unique_edges;
        for (auto edge : a_edges) {
            if (!a_edges_map_.count(edge))
                unique_edges.insert(edge);

            auto &nodes = a_edges_map_[edge];
            if (std::find(nodes.begin(), nodes.end(), name) == nodes.end())
                nodes.push_back(name);
        }
        return addVertexInternal(name, unique_edges, a_edges, start_coord,
                                 end_coord, type);
  }

    std::shared_ptr<Vertex> addVertex(std::shared_ptr<Vertex> v) { return addVertexInternal(v); }
    Edge *addEdge(Edge *e) { return addEdgeInternal(e); }

    const std::set<Edge *> &getEdgeSet(std::shared_ptr<Vertex> v) const { return v->edges_; }
    std::shared_ptr<Vertex> getVertex(const std::string &name) const { return node_map_.at(name); }
    const std::set<std::shared_ptr<Vertex>, cmpshared_ptrs> &getVertexSet() const { return nodes_; }
    std::set<std::shared_ptr<Vertex>, cmpshared_ptrs> const &getNodeSet() const { return nodes_; }

    void ExportToDot(const std::string &output_path) const;

    void makeRC(const std::string &first_vertex,
                const std::string &second_vertex) {
        node_map_[first_vertex]->rc_ = node_map_[second_vertex];
        node_map_[second_vertex]->rc_ = node_map_[first_vertex];
    }

    void FindBasicStatistic(std::ofstream &stat_stream);
    void FindDomainOrderings(debruijn_graph::GraphPack &gp,
                             const std::string &output_filename, const std::string &output_dir);
    void ExportPaths(debruijn_graph::GraphPack &gp,
                     const std::string &output_dir);

private:
    std::shared_ptr<Vertex>
    addVertexInternal(const std::string &name,
                      const std::set<debruijn_graph::EdgeId> &unique_edges,
                      const std::vector<debruijn_graph::EdgeId> &edges,
                      size_t start_coord, size_t end_coord, std::string type) {
        std::shared_ptr<Vertex> node =
                std::make_shared<Vertex>(unique_edges, edges, start_coord, end_coord, type);
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

    void OutputComponent(debruijn_graph::GraphPack &gp,
                         path_extend::BidirectionalPath *p, int component_id,
                         int ordering_id);

    std::string PathToSequence(path_extend::BidirectionalPath *p,
                               std::vector<std::shared_ptr<Vertex>> &answer);

    void OutputStatArrangement(const debruijn_graph::Graph &g,
                               std::vector<std::shared_ptr<Vertex>> single_candidate,
                               int id, std::ofstream &stat_file);
    void OutputStat(std::set<std::shared_ptr<Vertex>> &preliminary_visited,
                    std::ofstream &stat_file);

    size_t GetMaxVisited(const debruijn_graph::Graph &g, std::shared_ptr<Vertex> v,
                         double base_coverage) {
        double low_coverage = std::numeric_limits<double>::max();
        for (auto e : v->domain_edges_in_row_)
            low_coverage = std::min(low_coverage, g.coverage(e));

        return size_t(round(low_coverage / base_coverage));
  }

    void SetCopynumber(const debruijn_graph::Graph &g,
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

    void FindAllPossibleArrangements(const debruijn_graph::Graph &g, std::shared_ptr<Vertex> v,
                                     std::vector<std::vector<std::shared_ptr<Vertex>>> &answer,
                                     std::ofstream &stat_file);

    void PrelimDFS(std::shared_ptr<Vertex> v,
                   std::set<std::shared_ptr<Vertex>> &preliminary_visited);

    void FinalDFS(std::shared_ptr<Vertex> v,
                  std::vector<std::shared_ptr<Vertex>> &current,
                  std::set<std::shared_ptr<Vertex>> preliminary_visited,
                  std::vector<std::vector<std::shared_ptr<Vertex>>> &answer,
                  size_t component_size);



    std::set<std::shared_ptr<Vertex>, cmpshared_ptrs> nodes_; /* The set of nodes in the graph */
    std::set<Edge *> arcs_;                   /* The set of arcs in the graph  */
    std::map<std::string, std::shared_ptr<Vertex>> node_map_; /* A map from names to nodes     */
    std::map<EdgeId, std::vector<std::string>> a_edges_map_; /* A map from edges to nodenames */
    path_extend::PathContainer contig_paths_;
    size_t iteration_number_;

    DECL_LOGGER("AGraph2");
};

}
