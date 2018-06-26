#pragma once

#include "common/assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/core/graph.hpp"

namespace path_extend {

namespace scaffold_graph {

enum ScaffoldVertexT { Edge = 0, Path = 1 };

class InnerScaffoldVertex {
 public:
    virtual ~InnerScaffoldVertex() {};

    virtual size_t getId() const = 0;
    virtual ScaffoldVertexT getType() const = 0;
    virtual size_t getLengthFromGraph(const debruijn_graph::Graph &g) const = 0;
    virtual double getCoverageFromGraph(const debruijn_graph::Graph &g) const = 0;
    virtual shared_ptr<InnerScaffoldVertex> getConjugateFromGraph(const debruijn_graph::Graph &g) const = 0;
    virtual debruijn_graph::VertexId getEndGraphVertex(const debruijn_graph::Graph &g) const = 0;
    virtual debruijn_graph::VertexId getStartGraphVertex(const debruijn_graph::Graph &g) const = 0;
    virtual boost::optional<debruijn_graph::EdgeId> getLastEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const = 0;
    virtual boost::optional<debruijn_graph::EdgeId> getFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const = 0;

    virtual debruijn_graph::EdgeId getLastEdge() const = 0;
    virtual EdgeId getFirstEdge() const = 0;

    virtual BidirectionalPath* toPath(const debruijn_graph::Graph &g) const = 0;
    virtual BidirectionalPath getPath(const debruijn_graph::Graph &g) const = 0;
    virtual string str(const debruijn_graph::Graph& g) const = 0;
};

class EdgeIdVertex : public InnerScaffoldVertex {
 private:
    debruijn_graph::EdgeId edge_;

 public:
    explicit EdgeIdVertex(EdgeId edge_);

    size_t getId() const override;
    ScaffoldVertexT getType() const override;
    size_t getLengthFromGraph(const debruijn_graph::Graph &g) const override;
    double getCoverageFromGraph(const debruijn_graph::Graph &g) const override;
    shared_ptr<InnerScaffoldVertex> getConjugateFromGraph(const debruijn_graph::Graph &g) const override;
    debruijn_graph::VertexId getEndGraphVertex(const debruijn_graph::Graph &g) const override;
    debruijn_graph::VertexId getStartGraphVertex(const debruijn_graph::Graph &g) const override;
    optional<EdgeId> getLastEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const override;
    optional<EdgeId> getFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const override;

    EdgeId getLastEdge() const override;
    EdgeId getFirstEdge() const override;

    string str(const debruijn_graph::Graph &g) const override;
    BidirectionalPath* toPath(const debruijn_graph::Graph &g) const override;
    BidirectionalPath getPath(const debruijn_graph::Graph &g) const override;

    EdgeId get() const;
};

class PathVertex : public InnerScaffoldVertex {
 private:
    BidirectionalPath *path_;

 public:
    explicit PathVertex(BidirectionalPath *path_);

    size_t getId() const override;
    ScaffoldVertexT getType() const override;
    size_t getLengthFromGraph(const debruijn_graph::Graph &g) const override;
    double getCoverageFromGraph(const debruijn_graph::Graph &g) const override;
    shared_ptr<InnerScaffoldVertex> getConjugateFromGraph(const debruijn_graph::Graph &g) const override;
    VertexId getEndGraphVertex(const debruijn_graph::Graph &g) const override;
    VertexId getStartGraphVertex(const debruijn_graph::Graph &g) const override;
    optional<EdgeId> getLastEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const override;
    optional<EdgeId> getFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const override;

    EdgeId getLastEdge() const override;
    EdgeId getFirstEdge() const override;

    string str(const debruijn_graph::Graph &g) const override;
    BidirectionalPath* toPath(const debruijn_graph::Graph &g) const override;
    BidirectionalPath getPath(const debruijn_graph::Graph &g) const override;

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

    ScaffoldVertexT getType() const;
    size_t getLengthFromGraph(const debruijn_graph::Graph &g) const;
    double getCoverageFromGraph(const debruijn_graph::Graph &g) const;
    ScaffoldVertex getConjugateFromGraph(const debruijn_graph::Graph &g) const;
    debruijn_graph::VertexId getEndGraphVertex(const debruijn_graph::Graph &g) const;
    debruijn_graph::VertexId getStartGraphVertex(const debruijn_graph::Graph &g) const;
    boost::optional<debruijn_graph::EdgeId> getLastEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const;
    boost::optional<debruijn_graph::EdgeId> getFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const;

    debruijn_graph::EdgeId getLastEdge() const;
    debruijn_graph::EdgeId getFirstEdge() const;

    string str(const debruijn_graph::Graph &g) const;
    BidirectionalPath* toPath(const debruijn_graph::Graph &g) const;
    BidirectionalPath getPath(const debruijn_graph::Graph &g) const;

    shared_ptr<InnerScaffoldVertex> getInnerVertex() const;

    bool operator==(const ScaffoldVertex &rhs) const;
    bool operator!=(const ScaffoldVertex &rhs) const;
    bool operator<(const ScaffoldVertex &rhs) const;
    bool operator>(const ScaffoldVertex &rhs) const;
    bool operator<=(const ScaffoldVertex &rhs) const;
    bool operator>=(const ScaffoldVertex &rhs) const;
};

//fixme tmp solution to EdgeId-based scaffolder algorithms, should be removed
class EdgeGetter {
 public:

    EdgeId GetEdgeFromScaffoldVertex(const ScaffoldVertex& vertex) {
        VERIFY(vertex.getType() == Edge);
        auto inner_vertex = std::static_pointer_cast<EdgeIdVertex>(vertex.getInnerVertex());
        return inner_vertex->get();
    }
};

class PathGetter {
 public:

    BidirectionalPath* GetPathFromScaffoldVertex(const ScaffoldVertex& vertex) {
        VERIFY(vertex.getType() == Path);
        auto inner_vertex = std::static_pointer_cast<PathVertex>(vertex.getInnerVertex());
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
