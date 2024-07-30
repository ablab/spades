//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/core/graph.hpp"

#include "boost/optional/optional.hpp"

#include <unordered_set>

namespace scaffold_graph {

enum ScaffoldVertexT { Edge = 0, Path = 1 };

class InnerScaffoldVertex {
 public:
    using EdgeId = debruijn_graph::EdgeId;

    virtual ~InnerScaffoldVertex() = default;

    virtual size_t GetId() const = 0;
    virtual ScaffoldVertexT GetType() const = 0;
    virtual size_t GetLengthFromGraph(const debruijn_graph::Graph &g) const = 0;
    virtual double GetCoverageFromGraph(const debruijn_graph::Graph &g) const = 0;
    virtual std::shared_ptr<InnerScaffoldVertex> GetConjugateFromGraph(const debruijn_graph::Graph &g) const = 0;
    virtual debruijn_graph::VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const = 0;
    virtual debruijn_graph::VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const = 0;
    virtual boost::optional<EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const = 0;
    virtual boost::optional<EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const = 0;
    virtual std::string GetSequence(const debruijn_graph::Graph &g) const = 0;
    virtual size_t GetSize() const = 0;

    virtual debruijn_graph::EdgeId GetLastEdge() const = 0;
    virtual EdgeId GetFirstEdge() const = 0;
    virtual std::unordered_set<EdgeId> GetAllEdges() const = 0;

    virtual path_extend::BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const = 0;
    virtual std::string str(const debruijn_graph::Graph& g) const = 0;
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
    std::shared_ptr<InnerScaffoldVertex> GetConjugateFromGraph(const debruijn_graph::Graph &g) const override;
    debruijn_graph::VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const override;
    debruijn_graph::VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const override;
    boost::optional<EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    boost::optional<EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    std::string GetSequence(const debruijn_graph::Graph &g) const override;
    size_t GetSize() const override;

    EdgeId GetLastEdge() const override;
    EdgeId GetFirstEdge() const override;
    std::unordered_set<EdgeId> GetAllEdges() const override;

    std::string str(const debruijn_graph::Graph &g) const override;
    path_extend::BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const override;

    EdgeId get() const;
};

class PathVertex : public InnerScaffoldVertex {
 private:
    path_extend::BidirectionalPath *path_;

 public:
    using VertexId = debruijn_graph::VertexId;
    using InnerScaffoldVertex::EdgeId;

    explicit PathVertex(path_extend::BidirectionalPath *path_);

    size_t GetId() const override;
    ScaffoldVertexT GetType() const override;
    size_t GetLengthFromGraph(const debruijn_graph::Graph &g) const override;
    double GetCoverageFromGraph(const debruijn_graph::Graph &g) const override;
    std::shared_ptr<InnerScaffoldVertex> GetConjugateFromGraph(const debruijn_graph::Graph &g) const override;
    VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const override;
    VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const override;
    boost::optional<EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    boost::optional<EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    std::string GetSequence(const debruijn_graph::Graph &g) const override;
    size_t GetSize() const override;

    EdgeId GetLastEdge() const override;
    EdgeId GetFirstEdge() const override;
    std::unordered_set<EdgeId> GetAllEdges() const override;

    std::string str(const debruijn_graph::Graph &g) const override;
    path_extend::BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const override;

    path_extend::BidirectionalPath *get() const;
};

class ScaffoldVertex {
    std::shared_ptr<InnerScaffoldVertex> vertex_ptr_;

 public:
    using EdgeId = debruijn_graph::EdgeId;
    using VertexId = debruijn_graph::VertexId;

    explicit ScaffoldVertex(std::shared_ptr<InnerScaffoldVertex> vertex_ptr_);

    ScaffoldVertex(const ScaffoldVertex& other) = default;

    //make implicit for easy scaffold edge construction
    ScaffoldVertex(EdgeId edge);
    ScaffoldVertex(path_extend::BidirectionalPath *path);

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

    std::string str(const debruijn_graph::Graph &g) const;
    path_extend::BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const;

    std::shared_ptr<InnerScaffoldVertex> GetInnerVertex() const;

    bool operator==(const ScaffoldVertex &rhs) const;
    bool operator!=(const ScaffoldVertex &rhs) const;
    bool operator<(const ScaffoldVertex &rhs) const;
    bool operator>(const ScaffoldVertex &rhs) const;
    bool operator<=(const ScaffoldVertex &rhs) const;
    bool operator>=(const ScaffoldVertex &rhs) const;
};

class EdgeGetter {
 public:
    debruijn_graph::EdgeId GetEdgeFromScaffoldVertex(const ScaffoldVertex& vertex) {
        VERIFY_DEV(vertex.GetType() == Edge);
        auto inner_vertex = std::static_pointer_cast<EdgeIdVertex>(vertex.GetInnerVertex());
        return inner_vertex->get();
    }
};

}

namespace std {
template<>
struct hash<scaffold_graph::ScaffoldVertex> {
  size_t operator()(const scaffold_graph::ScaffoldVertex& vertex) const {
      return vertex.int_id();
  }
};

template<>
struct less<scaffold_graph::ScaffoldVertex> {
  bool operator()(const scaffold_graph::ScaffoldVertex& lhs,
                  const scaffold_graph::ScaffoldVertex& rhs) const {
      return lhs < rhs;
  }
};

}
