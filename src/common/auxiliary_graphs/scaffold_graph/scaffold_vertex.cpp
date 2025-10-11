//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_vertex.hpp"

#include "assembly_graph/paths/bidirectional_path_io/io_support.hpp"

namespace scaffold_graph {

std::shared_ptr<InnerScaffoldVertex> EdgeIdVertex::GetConjugateFromGraph(const debruijn_graph::Graph &g) const {
    return std::make_shared<EdgeIdVertex>(g.conjugate(get()));
}

path_extend::BidirectionalPath* EdgeIdVertex::ToPath(const debruijn_graph::Graph &g) const {
    //fixme think of a better way
    VERIFY(false);
    path_extend::BidirectionalPath* result;
//    path_extend::BidirectionalPath* result = new path_extend::BidirectionalPath (g);
//    path_extend::Gap gap(0);
//    result->PushBack(edge_, gap);
    return result;
}

boost::optional<debruijn_graph::EdgeId> EdgeIdVertex::GetLastEdgeWithPredicate(
        const func::TypedPredicate<EdgeId> &pred) const {
    boost::optional<EdgeId> result;
    if (pred(edge_)) {
        result = edge_;
    }
    return result;
}
boost::optional<debruijn_graph::EdgeId> EdgeIdVertex::GetFirstEdgeWithPredicate(
        const func::TypedPredicate<EdgeId> &pred) const {
    boost::optional<EdgeId> result;
    if (pred(edge_)) {
        result = edge_;
    }
    return result;
}
std::unordered_set<debruijn_graph::EdgeId> EdgeIdVertex::GetAllEdges() const {
    std::unordered_set<EdgeId> result;
    result.insert(edge_);
    return result;
}
std::string EdgeIdVertex::GetSequence(const debruijn_graph::Graph &g) const {
    return g.EdgeNucls(edge_).str();
}

std::shared_ptr<InnerScaffoldVertex> PathVertex::GetConjugateFromGraph(const debruijn_graph::Graph &/*g*/) const {
    return std::make_shared<PathVertex>(get()->GetConjPath());
}

debruijn_graph::VertexId PathVertex::GetEndGraphVertex(const debruijn_graph::Graph &g) const {
    VERIFY(path_->Size() > 0);
    return g.EdgeEnd(path_->Back());
}
debruijn_graph::VertexId PathVertex::GetStartGraphVertex(const debruijn_graph::Graph &g) const {
    VERIFY(path_->Size() > 0);
    return g.EdgeStart(path_->Front());
}

debruijn_graph::EdgeId PathVertex::GetLastEdge() const {
    const size_t path_size = path_->Size();
    VERIFY(path_size > 0);
    return path_->Back();
}
debruijn_graph::EdgeId PathVertex::GetFirstEdge() const {
    const size_t path_size = path_->Size();
    VERIFY(path_size > 0);
    return path_->Front();
}
boost::optional<debruijn_graph::EdgeId> PathVertex::GetLastEdgeWithPredicate(
        const func::TypedPredicate<EdgeId> &pred) const {
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
boost::optional<debruijn_graph::EdgeId> PathVertex::GetFirstEdgeWithPredicate(
        const func::TypedPredicate<EdgeId> &pred) const {
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
std::unordered_set<debruijn_graph::EdgeId> PathVertex::GetAllEdges() const {
    std::unordered_set<EdgeId> result;
    for (const auto &edge: *path_) {
        result.insert(edge);
    }
    return result;
}
std::string PathVertex::GetSequence(const debruijn_graph::Graph &g) const {
    path_extend::ScaffoldSequenceMaker sequence_maker(g);
    return sequence_maker.MakeSequence(*path_);
}

ScaffoldVertex::ScaffoldVertex(std::shared_ptr<InnerScaffoldVertex> vertex_ptr_) : vertex_ptr_(vertex_ptr_) {}

ScaffoldVertex ScaffoldVertex::GetConjugateFromGraph(const debruijn_graph::Graph &g) const {
    auto inner_vertex = vertex_ptr_->GetConjugateFromGraph(g);
    ScaffoldVertex result(inner_vertex);
    return result;
}

ScaffoldVertex::ScaffoldVertex(EdgeId edge) : vertex_ptr_(std::make_shared<EdgeIdVertex>(edge)) {}
ScaffoldVertex::ScaffoldVertex(path_extend::BidirectionalPath *path) : vertex_ptr_(std::make_shared<PathVertex>(path)) {}
}
