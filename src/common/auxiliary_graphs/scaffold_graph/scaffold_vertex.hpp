//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
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
    explicit EdgeIdVertex(EdgeId edge_) : edge_(edge_) {}

    size_t GetId() const override { return edge_.int_id(); }
    ScaffoldVertexT GetType() const override { return Edge; }
    size_t GetLengthFromGraph(const debruijn_graph::Graph &g) const override { return g.length(edge_); }
    double GetCoverageFromGraph(const debruijn_graph::Graph &g) const override { return g.coverage(edge_); }
    std::shared_ptr<InnerScaffoldVertex> GetConjugateFromGraph(const debruijn_graph::Graph &g) const override;
    debruijn_graph::VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const override { return g.EdgeEnd(edge_); }
    debruijn_graph::VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const override { return g.EdgeStart(edge_); }
    boost::optional<EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    boost::optional<EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    std::string GetSequence(const debruijn_graph::Graph &g) const override;
    size_t GetSize() const override { return 1; }

    EdgeId GetLastEdge() const override { return edge_; }
    EdgeId GetFirstEdge() const override { return edge_; }
    std::unordered_set<EdgeId> GetAllEdges() const override;

    std::string str(const debruijn_graph::Graph &g) const override { return g.str(edge_); }
    path_extend::BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const override;

    EdgeId get() const { return edge_; }
};

class PathVertex : public InnerScaffoldVertex {
 private:
    path_extend::BidirectionalPath *path_;

 public:
    using VertexId = debruijn_graph::VertexId;
    using InnerScaffoldVertex::EdgeId;

    explicit PathVertex(path_extend::BidirectionalPath *path_) : path_(path_) {}

    size_t GetId() const override { return path_->GetId(); }
    ScaffoldVertexT GetType() const override { return Path; }
    size_t GetLengthFromGraph(const debruijn_graph::Graph &/*g*/) const override { return path_->Length(); }
    double GetCoverageFromGraph(const debruijn_graph::Graph &/*g*/) const override { return path_->Coverage(); }
    std::shared_ptr<InnerScaffoldVertex> GetConjugateFromGraph(const debruijn_graph::Graph &/*g*/) const override;
    VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const override;
    VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const override;
    boost::optional<EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    boost::optional<EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const override;
    std::string GetSequence(const debruijn_graph::Graph &g) const override;
    size_t GetSize() const override { return path_->Size(); }

    EdgeId GetLastEdge() const override;
    EdgeId GetFirstEdge() const override;
    std::unordered_set<EdgeId> GetAllEdges() const override;

    std::string str(const debruijn_graph::Graph &/*g*/) const override { return path_->str(); }
    path_extend::BidirectionalPath* ToPath(const debruijn_graph::Graph &/*g*/) const override { return path_; }

    path_extend::BidirectionalPath *get() const { return path_; }
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

    ScaffoldVertex() : vertex_ptr_(nullptr) {}

    //deviate from surrounding style to make compatible with generic graph algorithms
    size_t int_id() const { return vertex_ptr_->GetId(); }

    ScaffoldVertexT GetType() const { return vertex_ptr_->GetType(); }
    size_t GetLengthFromGraph(const debruijn_graph::Graph &g) const { return vertex_ptr_->GetLengthFromGraph(g); }
    double GetCoverageFromGraph(const debruijn_graph::Graph &g) const { return vertex_ptr_->GetCoverageFromGraph(g); }
    ScaffoldVertex GetConjugateFromGraph(const debruijn_graph::Graph &g) const;
    VertexId GetEndGraphVertex(const debruijn_graph::Graph &g) const { return vertex_ptr_->GetEndGraphVertex(g); }
    VertexId GetStartGraphVertex(const debruijn_graph::Graph &g) const { return vertex_ptr_->GetStartGraphVertex(g); }
    boost::optional<EdgeId> GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
        return vertex_ptr_->GetLastEdgeWithPredicate(pred);
    }
    boost::optional<EdgeId> GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
        return vertex_ptr_->GetFirstEdgeWithPredicate(pred);
    }
    std::string GetSequence(const debruijn_graph::Graph &g) const { return vertex_ptr_->GetSequence(g); }
    size_t GetSize() const { return vertex_ptr_->GetSize(); }

    EdgeId GetLastEdge() const { return vertex_ptr_->GetLastEdge(); }
    EdgeId GetFirstEdge() const { return vertex_ptr_->GetFirstEdge(); }
    std::unordered_set<EdgeId> GetAllEdges() const { return vertex_ptr_->GetAllEdges(); }

    std::string str(const debruijn_graph::Graph &g) const { return vertex_ptr_->str(g); }
    path_extend::BidirectionalPath* ToPath(const debruijn_graph::Graph &g) const { return vertex_ptr_->ToPath(g); }

    std::shared_ptr<InnerScaffoldVertex> GetInnerVertex() const { return vertex_ptr_; }

    bool operator==(const ScaffoldVertex &rhs) const { return GetType() == rhs.GetType() and int_id() == rhs.int_id(); }
    bool operator!=(const ScaffoldVertex &rhs) const { return !(rhs == *this); }
    bool operator<(const ScaffoldVertex &rhs) const {
        return GetType() < rhs.GetType() or (GetType() == rhs.GetType() and int_id() < rhs.int_id());
    }
    bool operator>(const ScaffoldVertex &rhs) const { return rhs < *this; }
    bool operator<=(const ScaffoldVertex &rhs) const { return !(rhs < *this); }
    bool operator>=(const ScaffoldVertex &rhs) const { return !(*this < rhs); }
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
