#include "scaffold_vertex.hpp"

namespace path_extend {

namespace scaffold_graph {

size_t EdgeIdVertex::getId() const {
    return edge_.int_id();
}
size_t EdgeIdVertex::getLengthFromGraph(const debruijn_graph::Graph &g) const {
    return g.length(edge_);
}
EdgeIdVertex::EdgeIdVertex(EdgeId edge_) : edge_(edge_) {}
shared_ptr<InnerScaffoldVertex> EdgeIdVertex::getConjugateFromGraph(const debruijn_graph::Graph &g) const {
    return make_shared<EdgeIdVertex>(g.conjugate(get()));
}
EdgeId EdgeIdVertex::get() const {
    return edge_;
}
ScaffoldVertexT EdgeIdVertex::getType() const {
    return Edge;
}
debruijn_graph::VertexId EdgeIdVertex::getStartGraphVertex(const debruijn_graph::Graph &g) const {
    return g.EdgeStart(edge_);
}
debruijn_graph::VertexId EdgeIdVertex::getEndGraphVertex(const debruijn_graph::Graph &g) const {
    return g.EdgeEnd(edge_);
}
string EdgeIdVertex::str(const debruijn_graph::Graph &g) const {
    return g.str(edge_);
}
double EdgeIdVertex::getCoverageFromGraph(const debruijn_graph::Graph &g) const {
    return g.coverage(edge_);
}
//fixme change to smart ptr along with paths in PathContainer
BidirectionalPath* EdgeIdVertex::toPath(const debruijn_graph::Graph &g) const {
    BidirectionalPath* result = new BidirectionalPath (g);
    Gap gap(0);
    result->PushBack(edge_, gap);
    return result;
}
EdgeId EdgeIdVertex::getLastEdge() const {
    return edge_;
}
EdgeId EdgeIdVertex::getFirstEdge() const {
    return edge_;
}
optional<EdgeId> EdgeIdVertex::getLastEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const {
    boost::optional<EdgeId> result;
    if (pred(edge_)) {
        result = edge_;
    }
    return result;
}
optional<EdgeId> EdgeIdVertex::getFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const {
    boost::optional<EdgeId> result;
    if (pred(edge_)) {
        result = edge_;
    }
    return result;
}
BidirectionalPath EdgeIdVertex::getPath(const debruijn_graph::Graph &g) const {
    BidirectionalPath result(g);
    Gap gap(0);
    result.PushBack(edge_, gap);
    return result;
}

size_t PathVertex::getId() const {
    return path_->GetId();
}
size_t PathVertex::getLengthFromGraph(const debruijn_graph::Graph &/*g*/) const {
    return path_->Length();
}
shared_ptr<InnerScaffoldVertex> PathVertex::getConjugateFromGraph(const debruijn_graph::Graph &/*g*/) const {
    return std::make_shared<PathVertex>(get()->GetConjPath());
}
PathVertex::PathVertex(BidirectionalPath *path_) : path_(path_) {}
BidirectionalPath *PathVertex::get() const {
    return path_;
}
ScaffoldVertexT PathVertex::getType() const {
    return Path;
}
VertexId PathVertex::getEndGraphVertex(const debruijn_graph::Graph &g) const {
    VERIFY(path_->Size() > 0);
    return g.EdgeEnd(path_->Back());
}
VertexId PathVertex::getStartGraphVertex(const debruijn_graph::Graph &g) const {
    VERIFY(path_->Size() > 0);
    return g.EdgeStart(path_->Front());
}
string PathVertex::str(const debruijn_graph::Graph &/*g*/) const {
    return path_->str();
}
double PathVertex::getCoverageFromGraph(const debruijn_graph::Graph &/*g*/) const {
    return path_->Coverage();
}
path_extend::BidirectionalPath* PathVertex::toPath(const debruijn_graph::Graph &/*g*/) const {
    return path_;
}
EdgeId PathVertex::getLastEdge() const {
    const size_t path_size = path_->Size();
    VERIFY(path_size > 0);
    return path_->Back();
}
EdgeId PathVertex::getFirstEdge() const {
    const size_t path_size = path_->Size();
    VERIFY(path_size > 0);
    return path_->Front();
}
optional<EdgeId> PathVertex::getLastEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const {
    boost::optional<EdgeId> result;
    for (int i = static_cast<int>(path_->Size()) - 1; i >= 0; --i) {
        EdgeId current = path_->At(i);
        if (pred(current)) {
            result = current;
            return result;
        }
    }
    return result;
}
optional<EdgeId> PathVertex::getFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId>& pred) const {
    boost::optional<EdgeId> result;
    for (size_t i = 0; i < path_->Size(); ++i) {
        EdgeId current = path_->At(i);
        if (pred(current)) {
            result = current;
            return result;
        }
    }
    return result;
}
BidirectionalPath PathVertex::getPath(const debruijn_graph::Graph &g) const {
    BidirectionalPath result(*path_);
    return result;
}

ScaffoldVertex::ScaffoldVertex(shared_ptr<InnerScaffoldVertex> vertex_ptr_) : vertex_ptr_(vertex_ptr_) {}
size_t ScaffoldVertex::int_id() const {
    return vertex_ptr_->getId();
}
size_t ScaffoldVertex::getLengthFromGraph(const Graph &g) const {
    return vertex_ptr_->getLengthFromGraph(g);
}
ScaffoldVertex ScaffoldVertex::getConjugateFromGraph(const debruijn_graph::Graph &g) const {
    auto inner_vertex = vertex_ptr_->getConjugateFromGraph(g);
    ScaffoldVertex result(inner_vertex);
    return result;
}
ScaffoldVertexT ScaffoldVertex::getType() const {
    return vertex_ptr_->getType();
}
debruijn_graph::VertexId ScaffoldVertex::getStartGraphVertex(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->getStartGraphVertex(g);
}
debruijn_graph::VertexId ScaffoldVertex::getEndGraphVertex(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->getEndGraphVertex(g);
}
ScaffoldVertex::ScaffoldVertex(EdgeId edge) : vertex_ptr_(std::make_shared<EdgeIdVertex>(edge)) {}
ScaffoldVertex::ScaffoldVertex(BidirectionalPath *path) : vertex_ptr_(std::make_shared<PathVertex>(path)) {}
bool ScaffoldVertex::operator==(const ScaffoldVertex &rhs) const {
    return getType() == rhs.getType() and int_id() == rhs.int_id();
}
bool ScaffoldVertex::operator!=(const ScaffoldVertex &rhs) const {
    return !(rhs == *this);
}
bool ScaffoldVertex::operator<(const ScaffoldVertex &rhs) const {
    return getType() < rhs.getType() or (getType() == rhs.getType() and int_id() < rhs.int_id());
}
bool ScaffoldVertex::operator>(const ScaffoldVertex &rhs) const {
    return rhs < *this;
}
bool ScaffoldVertex::operator<=(const ScaffoldVertex &rhs) const {
    return !(rhs < *this);
}
bool ScaffoldVertex::operator>=(const ScaffoldVertex &rhs) const {
    return !(*this < rhs);
}
string ScaffoldVertex::str(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->str(g);
}
double ScaffoldVertex::getCoverageFromGraph(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->getCoverageFromGraph(g);
}
shared_ptr<InnerScaffoldVertex> ScaffoldVertex::getInnerVertex() const {
    return vertex_ptr_;
}
BidirectionalPath* ScaffoldVertex::toPath(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->toPath(g);
}
debruijn_graph::EdgeId ScaffoldVertex::getFirstEdge() const {
    return vertex_ptr_->getFirstEdge();
}
debruijn_graph::EdgeId ScaffoldVertex::getLastEdge() const {
    return vertex_ptr_->getLastEdge();
}
ScaffoldVertex::ScaffoldVertex(): vertex_ptr_(nullptr) {}
boost::optional<debruijn_graph::EdgeId> ScaffoldVertex::getLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
    return vertex_ptr_->getLastEdgeWithPredicate(pred);
}
boost::optional<debruijn_graph::EdgeId> ScaffoldVertex::getFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
    return vertex_ptr_->getFirstEdgeWithPredicate(pred);
}
BidirectionalPath ScaffoldVertex::getPath(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->getPath(g);
}

}

}