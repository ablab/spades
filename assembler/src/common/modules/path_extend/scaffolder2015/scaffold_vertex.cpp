#include <common/assembly_graph/paths/bidirectional_path_io/io_support.hpp>
#include "scaffold_vertex.hpp"
#include "scaffold_graph_constructor.hpp"

namespace path_extend {

namespace scaffold_graph {

size_t EdgeIdVertex::GetId() const {
    return edge_.int_id();
}
size_t EdgeIdVertex::GetLengthFromGraph(const debruijn_graph::Graph &g) const {
    return g.length(edge_);
}
EdgeIdVertex::EdgeIdVertex(EdgeId edge_) : edge_(edge_) {}
shared_ptr<InnerScaffoldVertex> EdgeIdVertex::GetConjugateFromGraph(const debruijn_graph::Graph &g) const {
    return make_shared<EdgeIdVertex>(g.conjugate(get()));
}
EdgeId EdgeIdVertex::get() const {
    return edge_;
}
ScaffoldVertexT EdgeIdVertex::GetType() const {
    return Edge;
}
debruijn_graph::VertexId EdgeIdVertex::GetStartGraphVertex(const debruijn_graph::Graph &g) const {
    return g.EdgeStart(edge_);
}
debruijn_graph::VertexId EdgeIdVertex::GetEndGraphVertex(const debruijn_graph::Graph &g) const {
    return g.EdgeEnd(edge_);
}
string EdgeIdVertex::str(const debruijn_graph::Graph &g) const {
    return g.str(edge_);
}
double EdgeIdVertex::GetCoverageFromGraph(const debruijn_graph::Graph &g) const {
    return g.coverage(edge_);
}
//fixme change to smart ptr along with paths in PathContainer
BidirectionalPath* EdgeIdVertex::ToPath(const debruijn_graph::Graph &g) const {
    BidirectionalPath* result = new BidirectionalPath (g);
    Gap gap(0);
    result->PushBack(edge_, gap);
    return result;
}
EdgeId EdgeIdVertex::GetLastEdge() const {
    return edge_;
}
EdgeId EdgeIdVertex::GetFirstEdge() const {
    return edge_;
}
optional<EdgeId> EdgeIdVertex::GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
    boost::optional<EdgeId> result;
    if (pred(edge_)) {
        result = edge_;
    }
    return result;
}
optional<EdgeId> EdgeIdVertex::GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
    boost::optional<EdgeId> result;
    if (pred(edge_)) {
        result = edge_;
    }
    return result;
}
BidirectionalPath EdgeIdVertex::GetPath(const debruijn_graph::Graph &g) const {
    BidirectionalPath result(g);
    Gap gap(0);
    result.PushBack(edge_, gap);
    return result;
}
std::unordered_set<EdgeId> EdgeIdVertex::GetAllEdges() const {
    unordered_set<EdgeId> result;
    result.insert(edge_);
    return result;
}
std::string EdgeIdVertex::GetSequence(const debruijn_graph::Graph &g) const {
    return g.EdgeNucls(edge_).str();
}

size_t PathVertex::GetId() const {
    return path_->GetId();
}
size_t PathVertex::GetLengthFromGraph(const debruijn_graph::Graph &/*g*/) const {
    return path_->Length();
}
shared_ptr<InnerScaffoldVertex> PathVertex::GetConjugateFromGraph(const debruijn_graph::Graph &/*g*/) const {
    return std::make_shared<PathVertex>(get()->GetConjPath());
}
PathVertex::PathVertex(BidirectionalPath *path_) : path_(path_) {}
BidirectionalPath *PathVertex::get() const {
    return path_;
}
ScaffoldVertexT PathVertex::GetType() const {
    return Path;
}
VertexId PathVertex::GetEndGraphVertex(const debruijn_graph::Graph &g) const {
    VERIFY(path_->Size() > 0);
    return g.EdgeEnd(path_->Back());
}
VertexId PathVertex::GetStartGraphVertex(const debruijn_graph::Graph &g) const {
    VERIFY(path_->Size() > 0);
    return g.EdgeStart(path_->Front());
}
string PathVertex::str(const debruijn_graph::Graph &/*g*/) const {
    return path_->str();
}
double PathVertex::GetCoverageFromGraph(const debruijn_graph::Graph &/*g*/) const {
    return path_->Coverage();
}
path_extend::BidirectionalPath* PathVertex::ToPath(const debruijn_graph::Graph &/*g*/) const {
    return path_;
}
EdgeId PathVertex::GetLastEdge() const {
    const size_t path_size = path_->Size();
    VERIFY(path_size > 0);
    return path_->Back();
}
EdgeId PathVertex::GetFirstEdge() const {
    const size_t path_size = path_->Size();
    VERIFY(path_size > 0);
    return path_->Front();
}
optional<EdgeId> PathVertex::GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
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
optional<EdgeId> PathVertex::GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
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
BidirectionalPath PathVertex::GetPath(const debruijn_graph::Graph &/*g*/) const {
    BidirectionalPath result(*path_);
    return result;
}
unordered_set<EdgeId> PathVertex::GetAllEdges() const {
    unordered_set<EdgeId> result;
    for (const auto &edge: *path_) {
        result.insert(edge);
    }
    return result;
}
std::string PathVertex::GetSequence(const debruijn_graph::Graph &g) const {
    //fixme might be ineffective
    path_extend::ScaffoldSequenceMaker sequence_maker(g);
    return sequence_maker.MakeSequence(*path_);
}

ScaffoldVertex::ScaffoldVertex(shared_ptr<InnerScaffoldVertex> vertex_ptr_) : vertex_ptr_(vertex_ptr_) {}
size_t ScaffoldVertex::int_id() const {
    return vertex_ptr_->GetId();
}
size_t ScaffoldVertex::GetLengthFromGraph(const Graph &g) const {
    return vertex_ptr_->GetLengthFromGraph(g);
}
ScaffoldVertex ScaffoldVertex::GetConjugateFromGraph(const debruijn_graph::Graph &g) const {
    auto inner_vertex = vertex_ptr_->GetConjugateFromGraph(g);
    ScaffoldVertex result(inner_vertex);
    return result;
}
ScaffoldVertexT ScaffoldVertex::GetType() const {
    return vertex_ptr_->GetType();
}
debruijn_graph::VertexId ScaffoldVertex::GetStartGraphVertex(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->GetStartGraphVertex(g);
}
debruijn_graph::VertexId ScaffoldVertex::GetEndGraphVertex(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->GetEndGraphVertex(g);
}
ScaffoldVertex::ScaffoldVertex(EdgeId edge) : vertex_ptr_(std::make_shared<EdgeIdVertex>(edge)) {}
ScaffoldVertex::ScaffoldVertex(BidirectionalPath *path) : vertex_ptr_(std::make_shared<PathVertex>(path)) {}
bool ScaffoldVertex::operator==(const ScaffoldVertex &rhs) const {
    return GetType() == rhs.GetType() and int_id() == rhs.int_id();
}
bool ScaffoldVertex::operator!=(const ScaffoldVertex &rhs) const {
    return !(rhs == *this);
}
bool ScaffoldVertex::operator<(const ScaffoldVertex &rhs) const {
    return GetType() < rhs.GetType() or (GetType() == rhs.GetType() and int_id() < rhs.int_id());
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
double ScaffoldVertex::GetCoverageFromGraph(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->GetCoverageFromGraph(g);
}
shared_ptr<InnerScaffoldVertex> ScaffoldVertex::GetInnerVertex() const {
    return vertex_ptr_;
}
BidirectionalPath* ScaffoldVertex::ToPath(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->ToPath(g);
}
debruijn_graph::EdgeId ScaffoldVertex::GetFirstEdge() const {
    return vertex_ptr_->GetFirstEdge();
}
debruijn_graph::EdgeId ScaffoldVertex::GetLastEdge() const {
    return vertex_ptr_->GetLastEdge();
}
ScaffoldVertex::ScaffoldVertex(): vertex_ptr_(nullptr) {}
boost::optional<debruijn_graph::EdgeId> ScaffoldVertex::GetLastEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
    return vertex_ptr_->GetLastEdgeWithPredicate(pred);
}
boost::optional<debruijn_graph::EdgeId> ScaffoldVertex::GetFirstEdgeWithPredicate(const func::TypedPredicate<EdgeId> &pred) const {
    return vertex_ptr_->GetFirstEdgeWithPredicate(pred);
}
BidirectionalPath ScaffoldVertex::GetPath(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->GetPath(g);
}
std::unordered_set<EdgeId> ScaffoldVertex::GetAllEdges() const {
    return vertex_ptr_->GetAllEdges();
}
std::string ScaffoldVertex::GetSequence(const debruijn_graph::Graph &g) const {
    return vertex_ptr_->GetSequence(g);
}
}

}