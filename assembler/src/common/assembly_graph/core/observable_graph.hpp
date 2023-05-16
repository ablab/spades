//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/logger/logger.hpp"
#include "graph_core.hpp"
#include "graph_iterators.hpp"

#include <vector>
#include <set>
#include <cstring>

namespace omnigraph {

template<class DataMaster>
class ObservableGraph: public GraphCore<DataMaster> {
public:
    typedef GraphCore<DataMaster> base;
    typedef typename base::DataMasterT DataMasterT;
    typedef typename base::VertexData VertexData;
    typedef typename base::EdgeData EdgeData;
    typedef typename base::EdgeId EdgeId;
    typedef typename base::VertexId VertexId;
    typedef typename base::VertexIt VertexIt;
    typedef typename base::edge_const_iterator edge_const_iterator;

    typedef HandlerApplier<VertexId, EdgeId> Applier;
    typedef SmartVertexIterator<ObservableGraph> SmartVertexIt;
    typedef SmartEdgeIterator<ObservableGraph> SmartEdgeIt;
    typedef ConstEdgeIterator<ObservableGraph> ConstEdgeIt;
    typedef ActionHandler<VertexId, EdgeId> Handler;

private:
   //todo switch to smart iterators
   mutable std::vector<Handler*> action_handler_list_;
   std::unique_ptr<const HandlerApplier<VertexId, EdgeId>> applier_;

public:
//todo move to graph core
    typedef ConstructionHelper<DataMaster> HelperT;

    HelperT GetConstructionHelper() {
//      TODO: fix everything and restore this check
//      VERIFY(this->VerifyAllDetached());
        return HelperT(*this);
    }

    const Applier& GetHandlerApplier() const {
        return *applier_;
    }

    void AddActionHandler(Handler* action_handler) const;

    bool RemoveActionHandler(const Handler* action_handler) const;

    bool AllHandlersThreadSafe() const;

   // TODO: for debug. remove.
    void PrintHandlersNames() const;

    void FireGameOver() const;
    
   //todo make Fire* protected once again with helper friend class
    void FireAddVertex(VertexId v) const;

    void FireAddEdge(EdgeId e) const;

    void FireDeleteVertex(VertexId v) const;

    void FireDeleteEdge(EdgeId e) const;

    void FireMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) const;

    void FireGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) const;

    void FireSplit(EdgeId edge, EdgeId new_edge1, EdgeId new_edge2) const;

    bool VerifyAllDetached();

    //smart iterators
    template<typename Priority>
    SmartVertexIterator<ObservableGraph, Priority> SmartVertexBegin(
            const Priority& priority, bool canonical_only = false) const {
        return SmartVertexIterator<ObservableGraph, Priority>(*this, priority, canonical_only);
    }

    SmartVertexIterator<ObservableGraph> SmartVertexBegin(bool canonical_only = false) const {
        return SmartVertexIterator<ObservableGraph>(*this, adt::identity(), canonical_only);
    }

    template<typename Priority>
    SmartEdgeIterator<ObservableGraph, Priority> SmartEdgeBegin(
            const Priority& priority, bool canonical_only = false) const {
        return SmartEdgeIterator<ObservableGraph, Priority>(*this, priority, canonical_only);
    }

    SmartEdgeIterator<ObservableGraph> SmartEdgeBegin(bool canonical_only = false) const {
        return SmartEdgeIterator<ObservableGraph>(*this, adt::identity(), canonical_only);
    }

    ConstEdgeIterator<ObservableGraph> ConstEdgeBegin(bool canonical_only = false) const {
        return ConstEdgeIterator<ObservableGraph>(*this, canonical_only);
    }

    void FireDeletePath(const std::vector<EdgeId>& edges_to_delete, const std::vector<VertexId>& vertices_to_delete) const;

    ObservableGraph(const DataMaster& master) :
            base(master), applier_(new PairedHandlerApplier<ObservableGraph>(*this)) {
    }

    virtual ~ObservableGraph();

    /////////////////////////graph operations
    //adding/removing vertices and edges

    VertexId AddVertex(VertexData data, VertexId id1 = 0, VertexId id2 = 0);

    void DeleteVertex(VertexId v);

    void ForceDeleteVertex(VertexId v);
    
    using base::conjugate;

    EdgeId AddEdge(EdgeData data, EdgeId id1 = 0, EdgeId id2 = 0);
    EdgeId AddEdge(VertexId v1, VertexId v2, EdgeData data,
                   EdgeId id1 = 0, EdgeId id2 = 0);

    void DeleteEdge(EdgeId e);

    void DeleteAllOutgoing(VertexId v);

    void DeleteAllIncoming(VertexId v);

    void clear();

    void CompressVertex(VertexId v);

    EdgeId UnsafeCompressVertex(VertexId v);

    std::vector<EdgeId> EdgesToDelete(const std::vector<EdgeId>& path) const;

    std::vector<VertexId> VerticesToDelete(const std::vector<EdgeId>& path) const;

    std::vector<EdgeId> CorrectMergePath(const std::vector<EdgeId>& path) const;

    EdgeId MergePath(const std::vector<EdgeId> &path,
                     bool safe_merging = true,
                     std::vector<uint32_t> overlaps = std::vector<uint32_t>());

    std::pair<EdgeId, EdgeId> SplitEdge(EdgeId edge, size_t position);

    EdgeId GlueEdges(EdgeId edge1, EdgeId edge2);

private:
    DECL_LOGGER("ObservableGraph")
};

template<class DataMaster>
typename ObservableGraph<DataMaster>::VertexId
ObservableGraph<DataMaster>::AddVertex(VertexData data, VertexId id1, VertexId id2) {
    VertexId v = base::HiddenAddVertex(std::move(data), id1, id2);
    FireAddVertex(v);
    return v;
}

template<class DataMaster>
void ObservableGraph<DataMaster>::DeleteVertex(VertexId v) {
    VERIFY(base::IsDeadEnd(v) && base::IsDeadStart(v));
    VERIFY(v != VertexId());
    FireDeleteVertex(v);
    base::HiddenDeleteVertex(v);
}

template<class DataMaster>
void ObservableGraph<DataMaster>::ForceDeleteVertex(VertexId v) {
    DeleteAllOutgoing(v);
    DeleteAllIncoming(v);
    DeleteVertex(v);
}

template<class DataMaster>
typename ObservableGraph<DataMaster>::EdgeId
ObservableGraph<DataMaster>::AddEdge(VertexId v1, VertexId v2, EdgeData data,
                                     EdgeId id1, EdgeId id2) {
    EdgeId e = base::HiddenAddEdge(v1, v2, std::move(data), id1, id2);
    FireAddEdge(e);
    return e;
}

template<class DataMaster>
typename ObservableGraph<DataMaster>::EdgeId
ObservableGraph<DataMaster>::AddEdge(EdgeData data, EdgeId id1, EdgeId id2) {
    EdgeId e = base::HiddenAddEdge(std::move(data), id1, id2);
    FireAddEdge(e);
    return e;
}

template<class DataMaster>
void ObservableGraph<DataMaster>::DeleteEdge(EdgeId e) {
    FireDeleteEdge(e);
    base::HiddenDeleteEdge(e);
}

template<class DataMaster>
void ObservableGraph<DataMaster>::DeleteAllOutgoing(VertexId v) {
    while (base::OutgoingEdgeCount(v) > 0) {
        EdgeId edge = *base::out_begin(v);
        DeleteEdge(edge);
    }
}

template<class DataMaster>
void ObservableGraph<DataMaster>::DeleteAllIncoming(VertexId v) {
    while (base::IncomingEdgeCount(v) > 0) {
        EdgeId edge = *base::in_begin(v);
        DeleteEdge(edge);
    }
}

template<class DataMaster>
void ObservableGraph<DataMaster>::CompressVertex(VertexId v) {
    //VERIFY(CanCompressVertex(v));
    if (base::CanCompressVertex(v)) {
        UnsafeCompressVertex(v);
    } else {
        TRACE("Vertex " << base::str(v) << " can't be compressed");
    }
}

template<class DataMaster>
typename ObservableGraph<DataMaster>::EdgeId ObservableGraph<DataMaster>::UnsafeCompressVertex(VertexId v) {
    VERIFY(base::CanCompressVertex(v));
    std::vector<EdgeId> edges_to_merge;
    edges_to_merge.push_back(base::GetUniqueIncomingEdge(v));
    edges_to_merge.push_back(base::GetUniqueOutgoingEdge(v));
    return MergePath(edges_to_merge);
}

template<class DataMaster>
std::vector<typename ObservableGraph<DataMaster>::EdgeId>
        ObservableGraph<DataMaster>::EdgesToDelete(const std::vector<EdgeId> &path) const {
    std::set<EdgeId> edgesToDelete;
    edgesToDelete.insert(path[0]);
    for (size_t i = 0; i + 1 < path.size(); i++) {
        EdgeId e = path[i + 1];
        if (edgesToDelete.find(base::conjugate(e)) == edgesToDelete.end())
            edgesToDelete.insert(e);
    }
    return std::vector<EdgeId>(edgesToDelete.begin(), edgesToDelete.end());
}

template<class DataMaster>
std::vector<typename ObservableGraph<DataMaster>::VertexId>
        ObservableGraph<DataMaster>::VerticesToDelete(const std::vector<EdgeId> &path) const {
    std::set<VertexId> verticesToDelete;
    for (size_t i = 0; i + 1 < path.size(); i++) {
        EdgeId e = path[i + 1];
        VertexId v = base::EdgeStart(e);
        if (verticesToDelete.find(base::conjugate(v)) == verticesToDelete.end())
            verticesToDelete.insert(v);
    }
    return std::vector<VertexId>(verticesToDelete.begin(), verticesToDelete.end());
}

template<class DataMaster>
void ObservableGraph<DataMaster>::AddActionHandler(Handler* action_handler) const {
#pragma omp critical(action_handler_list_modification)
    {
        TRACE("Action handler " << action_handler->name() << " added");
        if (std::find(action_handler_list_.begin(), action_handler_list_.end(), action_handler) != action_handler_list_.end()) {
            VERIFY_MSG(false, "Action handler " << action_handler->name() << " has already been added");
        } else {
            action_handler_list_.push_back(action_handler);
        }
    }
}

template<class DataMaster>
bool ObservableGraph<DataMaster>::RemoveActionHandler(const Handler* action_handler) const {
    bool result = false;
#pragma omp critical(action_handler_list_modification)
    {
        auto it = std::find(action_handler_list_.begin(), action_handler_list_.end(), action_handler);
        if (it != action_handler_list_.end()) {
            action_handler_list_.erase(it);
            TRACE("Action handler " << action_handler->name() << " removed");
            result = true;
        } else {
            TRACE("Action handler " << action_handler->name() << " wasn't found among graph action handlers");
        }
    }
    return result;
}

template<class DataMaster>
bool ObservableGraph<DataMaster>::AllHandlersThreadSafe() const {
    for (Handler* handler : action_handler_list_) {
        if (handler->IsAttached() && !handler->IsThreadSafe()) {
            return false;
        }
    }
    return true;
}

template<class DataMaster>
void ObservableGraph<DataMaster>::PrintHandlersNames() const {
    for (Handler* handler : action_handler_list_) {
        std::cout << handler->name() << " attached=" << handler->IsAttached() << std::endl;
    }
}

template<class DataMaster>
void ObservableGraph<DataMaster>::FireGameOver() const {
    auto action_handler_list_copy = action_handler_list_;  // We should copy list since GameOver initiates removing from the list
    for (Handler* handler_ptr : action_handler_list_copy) {
        TRACE("FireGameOver to handler " << handler_ptr->name());
        applier_->ApplyGameOver(*handler_ptr);
    }
}

template<class DataMaster>
void ObservableGraph<DataMaster>::FireAddVertex(VertexId v) const {
    for (Handler* handler_ptr : action_handler_list_) {
        if (handler_ptr->IsAttached()) {
            TRACE("FireAddVertex to handler " << handler_ptr->name());
            applier_->ApplyAdd(*handler_ptr, v);
        }
    }
}

template<class DataMaster>
void ObservableGraph<DataMaster>::FireAddEdge(EdgeId e) const {
    for (Handler* handler_ptr : action_handler_list_) {
        if (handler_ptr->IsAttached()) {
            TRACE("FireAddEdge to handler " << handler_ptr->name());
            applier_->ApplyAdd(*handler_ptr, e);
        }
    }
}

template<class DataMaster>
void ObservableGraph<DataMaster>::FireDeleteVertex(VertexId v) const {
    for (auto it = action_handler_list_.rbegin(); it != action_handler_list_.rend(); ++it) {
        if ((*it)->IsAttached()) {
            applier_->ApplyDelete(**it, v);
        }
    }
}

template<class DataMaster>
void ObservableGraph<DataMaster>::FireDeleteEdge(EdgeId e) const {
    for (auto it = action_handler_list_.rbegin(); it != action_handler_list_.rend(); ++it) {
        if ((*it)->IsAttached()) {
            applier_->ApplyDelete(**it, e);
        }
    };
}

template<class DataMaster>
void ObservableGraph<DataMaster>::FireMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) const {
    for (Handler* handler_ptr : action_handler_list_) {
        if (handler_ptr->IsAttached()) {
            applier_->ApplyMerge(*handler_ptr, old_edges, new_edge);
        }
    }
}

template<class DataMaster>
void ObservableGraph<DataMaster>::FireGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) const {
    for (Handler* handler_ptr : action_handler_list_) {
        if (handler_ptr->IsAttached()) {
            applier_->ApplyGlue(*handler_ptr, new_edge, edge1, edge2);
        }
    };
}

template<class DataMaster>
void ObservableGraph<DataMaster>::FireSplit(EdgeId edge, EdgeId new_edge1, EdgeId new_edge2) const {
    for (Handler* handler_ptr : action_handler_list_) {
        if (handler_ptr->IsAttached()) {
            applier_->ApplySplit(*handler_ptr, edge, new_edge1, new_edge2);
        }
    }
}

template<class DataMaster>
bool ObservableGraph<DataMaster>::VerifyAllDetached() {
    for (Handler* handler_ptr : action_handler_list_) {
        if (handler_ptr->IsAttached()) {
            return false;
        }
    }
    return true;
}

template<class DataMaster>
void ObservableGraph<DataMaster>::FireDeletePath(const std::vector<EdgeId> &edgesToDelete,
                                                 const std::vector<VertexId> &verticesToDelete) const {
    for (EdgeId e : edgesToDelete)
        FireDeleteEdge(e);
    for (VertexId v : verticesToDelete)
        FireDeleteVertex(v);
}

template<class DataMaster>
void ObservableGraph<DataMaster>::clear() {
    for (VertexId v : base::vertices())
        ForceDeleteVertex(v);
}

template<class DataMaster>
ObservableGraph<DataMaster>::~ObservableGraph<DataMaster>() {
    FireGameOver();
    clear();
}

template<class DataMaster>
std::vector<typename ObservableGraph<DataMaster>::EdgeId>
ObservableGraph<DataMaster>::CorrectMergePath(const std::vector<EdgeId>& path) const {
    for (size_t i = 0; i < path.size(); i++) {
        if (path[i] == base::conjugate(path[i])) {
            std::vector<EdgeId> result;
            if (i < path.size() - 1 - i) {
                for (size_t j = 0; j < path.size(); j++)
                    result.push_back(base::conjugate(path[path.size() - 1 - j]));
                i = path.size() - 1 - i;
            } else {
                result = path;
            }
            size_t size = 2 * i + 1;
            for (size_t j = result.size(); j < size; j++) {
                result.push_back(base::conjugate(result[size - 1 - j]));
            }
            return result;
        }
    }
    return path;
}

template<class DataMaster>
typename ObservableGraph<DataMaster>::EdgeId
        ObservableGraph<DataMaster>::MergePath(const std::vector<EdgeId> &path,
                                               bool safe_merging,
                                               std::vector<uint32_t> overlaps) {
    VERIFY(!path.empty());
    for (size_t i = 0; i < path.size(); i++)
        for (size_t j = i + 1; j < path.size(); j++) {
            VERIFY(path[i] != path[j]);
        }
    if (path.size() == 1) {
        TRACE(
                "Path of single edge " << base::str(*(path.begin())) << ". Nothing to merge.");
    };
    //      cerr << "Merging " << PrintDetailedPath(pObservableGraph<DataMaster><VertexIdT, EdgeIdT, VertexIt>ath) << endl;
    //      cerr << "Conjugate " << PrintConjugatePath(path) << endl;
    auto corrected_path = CorrectMergePath(path);
    VertexId v1 = base::EdgeStart(corrected_path[0]);
    VertexId v2 = base::EdgeEnd(corrected_path[corrected_path.size() - 1]);
    std::vector<const EdgeData *> to_merge;
    std::vector<uint32_t> local_overlaps;
    if (overlaps.empty()) {
        for (auto it1 = corrected_path.begin(), it2 = std::next(it1); it2 != corrected_path.end(); ++it1, ++it2) {
            to_merge.push_back(&(base::data(*it1)));
            VertexId end = base::EdgeEnd(*it1);
            VERIFY(end == base::EdgeStart(*it2));
            uint32_t overlap = base::data(end).overlap();
            overlaps.push_back(overlap);
        }
    }
    to_merge.push_back(&(base::data(corrected_path.back())));
    EdgeId new_edge = base::HiddenAddEdge(v1, v2, base::master().MergeData(to_merge, overlaps, safe_merging));
    FireMerge(corrected_path, new_edge);
    auto edges_to_delete = EdgesToDelete(corrected_path);
    auto vertices_to_delete = VerticesToDelete(corrected_path);
    FireDeletePath(edges_to_delete, vertices_to_delete);
    FireAddEdge(new_edge);
    base::HiddenDeletePath(edges_to_delete, vertices_to_delete);
    return new_edge;
}

template<class DataMaster>
std::pair<typename ObservableGraph<DataMaster>::EdgeId, typename ObservableGraph<DataMaster>::EdgeId>
        ObservableGraph<DataMaster>::SplitEdge(EdgeId edge, size_t position) {
    bool sc_flag = (edge == conjugate(edge));
    VERIFY_MSG(position > 0 && position < (sc_flag ? base::length(edge) / 2 + 1 : base::length(edge)),
            "Edge length is " << base::length(edge) << " but split pos was " << position);
    auto [vdata, edata1, edata2]  = base::master().SplitData(base::data(edge), position, sc_flag);
    VertexId splitVertex = base::HiddenAddVertex(std::move(vdata));
    EdgeId new_edge1 = base::HiddenAddEdge(base::EdgeStart(edge), splitVertex,
                                           std::move(edata1));
    EdgeId new_edge2 = base::HiddenAddEdge(splitVertex, sc_flag ? conjugate(splitVertex) : base::EdgeEnd(edge),
                                           std::move(edata2));
    VERIFY(!sc_flag || new_edge2 == conjugate(new_edge2))
    FireSplit(edge, new_edge1, new_edge2);
    FireDeleteEdge(edge);
    FireAddVertex(splitVertex);
    FireAddEdge(new_edge1);
    FireAddEdge(new_edge2);
    base::HiddenDeleteEdge(edge);
    return {new_edge1, new_edge2};
}

template<class DataMaster>
typename ObservableGraph<DataMaster>::EdgeId ObservableGraph<DataMaster>::GlueEdges(EdgeId edge1, EdgeId edge2) {
    EdgeId new_edge = base::HiddenAddEdge(base::EdgeStart(edge2), base::EdgeEnd(edge2), base::master().GlueData(base::data(edge1), base::data(edge2)));
    FireGlue(new_edge, edge1, edge2);
    FireDeleteEdge(edge1);
    FireDeleteEdge(edge2);
    FireAddEdge(new_edge);
    VertexId start = base::EdgeStart(edge1);
    VertexId end = base::EdgeEnd(edge1);
    base::HiddenDeleteEdge(edge1);
    base::HiddenDeleteEdge(edge2);

    if (base::IsDeadStart(start) && base::IsDeadEnd(start))
        DeleteVertex(start);

    if (base::IsDeadStart(end) && base::IsDeadEnd(end))
        DeleteVertex(end);

    DEBUG("Delete vertex");
    return new_edge;
}

} // namespace omnigraph
