//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/assembly_graph/paths/bidirectional_path.hpp"
#include "common/assembly_graph/core/graph.hpp"

#include "boost/optional/optional.hpp"

#include <unordered_set>

namespace path_extend {

namespace scaffold_graph {

enum ScaffoldVertexT { Edge = 0, Path = 1 };

class InnerScaffoldVertex {
 public:
    virtual ~InnerScaffoldVertex() = default;

    virtual size_t GetId() const = 0;
    virtual ScaffoldVertexT GetType() const = 0;
    virtual size_t GetLengthFromGraph(const debruijn_graph::Graph &g) const = 0;
    virtual double GetCoverageFromGraph(const debruijn_graph::Graph &g) const = 0;
    virtual shared_ptr<InnerScaffoldVertex> GetConjugateFromGraph(const debruijn_graph::Graph &g) const = 0;
    virtual debruijn_graph::VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const = 0;
    virtual debruijn_graph::VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const = 0;
    virtual boost::optional<debruijn_graph::EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const = 0;
    virtual boost::optional<debruijn_graph::EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const = 0;
    virtual string GetSequence(const debruijn_graph::Graph &g) const = 0;
    virtual size_t GetSize() const = 0;

    virtual debruijn_graph::EdgeId GetLastEdge() const = 0;
    virtual EdgeId GetFirstEdge() const = 0;
    virtual std::unordered_set<EdgeId> GetAllEdges() const = 0;

    virtual BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const = 0;
    virtual BidirectionalPath GetPath(const debruijn_graph::Graph &g) const = 0;
    virtual string str(const debruijn_graph::Graph& g) const = 0;
};

class EdgeIdVertex : public InnerScaffoldVertex {
 private:
    debruijn_graph::EdgeId edge_;

 public:
    explicit EdgeIdVertex(EdgeId edge_);

    size_t GetId() const override;
    ScaffoldVertexT GetType() const override;
    size_t GetLengthFromGraph(const debruijn_graph::Graph &g) const override;
    double GetCoverageFromGraph(const debruijn_graph::Graph &g) const override;
    shared_ptr<InnerScaffoldVertex> GetConjugateFromGraph(const debruijn_graph::Graph &g) const override;
    debruijn_graph::VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const override;
    debruijn_graph::VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const override;
    optional<EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    optional<EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    std::string GetSequence(const debruijn_graph::Graph &g) const override;
    size_t GetSize() const override;

    EdgeId GetLastEdge() const override;
    EdgeId GetFirstEdge() const override;
    std::unordered_set<EdgeId> GetAllEdges() const override;

    string str(const debruijn_graph::Graph &g) const override;
    BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const override;
    BidirectionalPath GetPath(const debruijn_graph::Graph &g) const override;

    EdgeId get() const;
};

class PathVertex : public InnerScaffoldVertex {
 private:
    BidirectionalPath *path_;

 public:
    explicit PathVertex(BidirectionalPath *path_);

    size_t GetId() const override;
    ScaffoldVertexT GetType() const override;
    size_t GetLengthFromGraph(const debruijn_graph::Graph &g) const override;
    double GetCoverageFromGraph(const debruijn_graph::Graph &g) const override;
    shared_ptr<InnerScaffoldVertex> GetConjugateFromGraph(const debruijn_graph::Graph &g) const override;
    VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const override;
    VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const override;
    optional<EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    optional<EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    std::string GetSequence(const debruijn_graph::Graph &g) const override;
    size_t GetSize() const override;

    EdgeId GetLastEdge() const override;
    EdgeId GetFirstEdge() const override;
    unordered_set<EdgeId> GetAllEdges() const override;

    string str(const debruijn_graph::Graph &g) const override;
    BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const override;
    BidirectionalPath GetPath(const debruijn_graph::Graph &g) const override;

    BidirectionalPath *get() const;
};

class ScaffoldVertex {
    shared_ptr<InnerScaffoldVertex> vertex_ptr_;

 public:
    explicit ScaffoldVertex(shared_ptr<InnerScaffoldVertex> vertex_ptr_);

    ScaffoldVertex(const ScaffoldVertex& other) = default;

    //make implicit for easy scaffold edge construction
    ScaffoldVertex(EdgeId edge);
    ScaffoldVertex(BidirectionalPath *path);

    ScaffoldVertex();

    //deviate from surrounding style to make compatible with generic graph algorithms
    size_t int_id() const;

    ScaffoldVertexT GetType() const;
    size_t GetLengthFromGraph(const debruijn_graph::Graph &g) const;
    double GetCoverageFromGraph(const debruijn_graph::Graph &g) const;
    ScaffoldVertex GetConjugateFromGraph(const debruijn_graph::Graph &g) const;
    VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const;
    VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const;
    boost::optional<EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const;
    boost::optional<EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const;
    std::string GetSequence(const debruijn_graph::Graph &g) const;
    size_t GetSize() const;

    EdgeId GetLastEdge() const;
    EdgeId GetFirstEdge() const;
    std::unordered_set<EdgeId> GetAllEdges() const;

    string str(const debruijn_graph::Graph &g) const;
    BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const;
    BidirectionalPath GetPath(const debruijn_graph::Graph &g) const;

    shared_ptr<InnerScaffoldVertex> GetInnerVertex() const;

    bool operator==(const ScaffoldVertex &rhs) const;
    bool operator!=(const ScaffoldVertex &rhs) const;
    bool operator<(const ScaffoldVertex &rhs) const;
    bool operator>(const ScaffoldVertex &rhs) const;
    bool operator<=(const ScaffoldVertex &rhs) const;
    bool operator>=(const ScaffoldVertex &rhs) const;
};

class EdgeGetter {
 public:
    EdgeId GetEdgeFromScaffoldVertex(const ScaffoldVertex& vertex) {
        VERIFY_DEV(vertex.GetType() == Edge);
        auto inner_vertex = std::static_pointer_cast<EdgeIdVertex>(vertex.GetInnerVertex());
        return inner_vertex->get();
    }
};

}
}

namespace std {
template<>
struct hash<path_extend::scaffold_graph::ScaffoldVertex> {
  size_t operator()(const path_extend::scaffold_graph::ScaffoldVertex& vertex) const {
      return vertex.int_id();
  }
};

template<>
struct less<path_extend::scaffold_graph::ScaffoldVertex> {
  bool operator()(const path_extend::scaffold_graph::ScaffoldVertex& lhs,
                  const path_extend::scaffold_graph::ScaffoldVertex& rhs) const {
      return lhs < rhs;
  }
};

}
