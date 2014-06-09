//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "logger/logger.hpp"
#include "omni_utils.hpp"
#include "abstract_graph.hpp"
#include "coverage.hpp"

#include <vector>
#include <set>
#include <cstring>

namespace omnigraph {

template<class DataMaster>
class AbstractConjugateGraph;

template<class DataMaster>
class PairedEdge;

template<class V>
class PairedVertexLock;

template<class Graph>
class ConstructionHelper {
private:
    typedef typename Graph::DataMaster::EdgeData EdgeData;
    typedef typename Graph::DataMaster::VertexData VertexData;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    Graph &graph_;

public:

    ConstructionHelper(Graph &graph) : graph_(graph) {
    }

    Graph &graph() {
        return graph_;
    }

    EdgeId AddEdge(const EdgeData &data) {
        return AddEdge(data, graph_.GetGraphIdDistributor());
    }

    EdgeId AddEdge(const EdgeData &data, IdDistributor &id_distributor) {
        return graph_.AddEdge(data, id_distributor);
    }

    void LinkIncomingEdge(VertexId v, EdgeId e) {
        graph_.LinkIncomingEdge(v, e);
    }

    void LinkOutgoingEdge(VertexId v, EdgeId e) {
        graph_.LinkOutgoingEdge(v, e);
    }

    void DeleteLink(VertexId v, EdgeId e) {
        v->RemoveOutgoingEdge(e);
    }

    void DeleteUnlinkedEdge(EdgeId e) {
        EdgeId rc = graph_.conjugate(e);
        if (e != rc) {
            delete rc.get();
        }
        delete e.get();
    }

    VertexId CreateVertex(const VertexData &data) {
        return CreateVertex(data, graph_.GetGraphIdDistributor());
    }

    VertexId CreateVertex(const VertexData &data, IdDistributor &id_distributor) {
        return graph_.CreateVertex(data, id_distributor);
    }

    template<class Iter>
    void AddVerticesToGraph(Iter begin, Iter end) {
        for(; begin != end; ++begin) {
            graph_.AddVertexToGraph(*begin);
        }
    }
};

template<class DataMaster>
class PairedVertex {
private:
    typedef restricted::pure_pointer<PairedVertex<DataMaster>> VertexId;
    typedef restricted::pure_pointer<PairedEdge<DataMaster>> EdgeId;
    typedef typename DataMaster::VertexData VertexData;
    typedef typename std::vector<EdgeId>::const_iterator edge_raw_iterator;

    class conjugate_iterator : public boost::iterator_facade<conjugate_iterator,
            EdgeId, boost::forward_traversal_tag, EdgeId> {
    public:
        explicit conjugate_iterator(edge_raw_iterator it,
                                    bool conjugate = false)
                : it_(it),
                  conjugate_(conjugate) {
        }

        //fixme seems to be ok and rather useful. Delete comment?
        // Should not exist. Temporary patch to write empty in_begin, out_end... methods for
        // ConcurrentGraphComponent which can not have such methods by definition
        conjugate_iterator()
                : conjugate_(false) {
    //        VERIFY_MSG(false, "There is no sense in using this. See comments.")
        }

    private:
        friend class boost::iterator_core_access;

        void increment() {
            it_++;
        }

        bool equal(const conjugate_iterator &other) const {
            return other.it_ == it_ && other.conjugate_ == conjugate_;
        }

        EdgeId dereference() const {
            return (conjugate_ ? (*it_)->conjugate() : *it_);
        }

        edge_raw_iterator it_;
        bool conjugate_;
    };

public:
    typedef conjugate_iterator edge_const_iterator;

private:
    friend class AbstractGraph<
            restricted::pure_pointer<PairedVertex<DataMaster>>,
            restricted::pure_pointer<PairedEdge<DataMaster>>, DataMaster> ;
    friend class AbstractConjugateGraph<DataMaster> ;
    friend class ConstructionHelper<AbstractConjugateGraph<DataMaster>>;
    friend class PairedEdge<DataMaster> ;
    friend class PairedVertexLock<restricted::pure_pointer<PairedVertex<DataMaster>>>;
    friend class conjugate_iterator;

    std::vector<EdgeId> outgoing_edges_;

    VertexId conjugate_;

    VertexData data_;

    bool IsMinimal() const {
        return conjugate_->conjugate_ <= conjugate_;
    }


    void set_conjugate(VertexId conjugate) {
        conjugate_ = conjugate;
    }

    size_t OutgoingEdgeCount() const {
        return outgoing_edges_.size();
    }

    edge_const_iterator out_begin() const {
        return edge_const_iterator(outgoing_edges_.cbegin(), false);
    }

    edge_const_iterator out_end() const {
        return edge_const_iterator(outgoing_edges_.cend(), false);
    }

    const vector<EdgeId> OutgoingEdgesTo(VertexId v) const {
        vector<EdgeId> result;
        for (auto it = outgoing_edges_.begin(); it != outgoing_edges_.end();
                ++it) {
            if ((*it)->end() == v) {
                result.push_back(*it);
            }
        }
        return result;
    }

    size_t IncomingEdgeCount() const {
        return conjugate_->OutgoingEdgeCount();
    }

    size_t IncomingEdgesCount() const {
        return (conjugate_->OutgoingEdgeCount());
    }

    edge_const_iterator in_begin() const {
        return edge_const_iterator(conjugate_->outgoing_edges_.cbegin(), true);
    }

    edge_const_iterator in_end() const {
        return edge_const_iterator(conjugate_->outgoing_edges_.cend(), true);
    }

    PairedVertex(VertexData data)
            : data_(data) {
    }

    VertexData &data() {
        return data_;
    }

    void set_data(VertexData data) {
        data_ = data;
    }

    void AddOutgoingEdge(EdgeId e) {
        VERIFY(this != 0);
        outgoing_edges_.insert(std::upper_bound(outgoing_edges_.begin(), outgoing_edges_.end(), e), e);
        //outgoing_edges_.push_back(e);
    }

    bool RemoveOutgoingEdge(const EdgeId e) {
        auto it = std::find(outgoing_edges_.begin(), outgoing_edges_.end(), e);
        if (it == outgoing_edges_.end())
            return false;

        outgoing_edges_.erase(it);
        return true;
    }

    VertexId conjugate() const {
        return conjugate_;
    }

    ~PairedVertex() {
        TRACE("PairedVertex destructor");
        VERIFY(outgoing_edges_.size() == 0);
        TRACE("PairedVertex destructor ok");
    }
};

template<class V>
class PairedVertexLock {
    restricted::PurePtrLock<V> inner_lock_;

    static bool IsMinimal(V v) {
        return !(v->conjugate_ < v);
    }

    static V MinimalFromPair(V v) {
        if (IsMinimal(v)) {
            return v;
        } else {
            return v->conjugate_;
        }
    }

    static V& GetLockableElement(V v) {
        return v->conjugate_;
    }

public:
    PairedVertexLock(V v) :
        inner_lock_(GetLockableElement(MinimalFromPair(v)))
    {
    }

};

template<class DataMaster>
class PairedEdge {
 private:
    typedef restricted::pure_pointer<PairedVertex<DataMaster>> VertexId;
    typedef restricted::pure_pointer<PairedEdge<DataMaster>> EdgeId;
    typedef typename DataMaster::EdgeData EdgeData;
    friend class AbstractGraph<
            restricted::pure_pointer<PairedVertex<DataMaster>>,
            restricted::pure_pointer<PairedEdge<DataMaster>>, DataMaster> ;
    friend class AbstractConjugateGraph<DataMaster> ;
    friend class ConstructionHelper<AbstractConjugateGraph<DataMaster>>;
    //todo unfriend
    friend class PairedVertex<DataMaster> ;
    VertexId end_;

    EdgeData data_;

    EdgeId conjugate_;

    PairedEdge(VertexId end, const EdgeData &data)
            : end_(end),
              data_(data) {
    }

    EdgeData &data() {
        return data_;
    }

    void set_data(const EdgeData &data) {
        data_ = data;
    }

    VertexId end() const {
        return end_;
    }

    VertexId start() const {
        return conjugate_->end()->conjugate();
    }

    void set_conjugate(EdgeId conjugate) {
        conjugate_ = conjugate;
    }

    void SetEndVertex(VertexId end) {
        end_ = end;
    }

    ~PairedEdge() {
    }

public:
    EdgeId conjugate() const {
        return conjugate_;
    }

    size_t length(size_t k) const {
        return data_.size() - k;
    }
};

template<class DataMaster>
class AbstractConjugateGraph : public AbstractGraph<
        restricted::pure_pointer<PairedVertex<DataMaster>>,
        restricted::pure_pointer<PairedEdge<DataMaster>>, DataMaster> {
private:
    typedef AbstractGraph<restricted::pure_pointer<PairedVertex<DataMaster>>,
            restricted::pure_pointer<PairedEdge<DataMaster>>, DataMaster> base;
    typedef restricted::IdDistributor IdDistributor;

    friend class ConstructionHelper<AbstractConjugateGraph<DataMaster>>;

public:

    typedef typename base::VertexId VertexId;
    typedef typename base::EdgeId EdgeId;
    typedef typename base::VertexData VertexData;
    typedef typename base::EdgeData EdgeData;
    typedef typename base::VertexIterator VertexIterator;
    typedef ConstructionHelper<AbstractConjugateGraph<DataMaster>> HelperT;
    using base::str;
    using base::CreateVertex;
    using base::HiddenAddEdge;

protected:

    VertexId CreateVertex(const VertexData &data1, const VertexData &data2, restricted::IdDistributor &id_distributor) {
        VertexId vertex1(new PairedVertex<DataMaster>(data1), id_distributor);
        VertexId vertex2(new PairedVertex<DataMaster>(data2), id_distributor);

        vertex1->set_conjugate(vertex2);
        vertex2->set_conjugate(vertex1);

        return vertex1;
    }

    virtual VertexId CreateVertex(const VertexData &data, restricted::IdDistributor &id_distributor) {
        return CreateVertex(data, this->master().conjugate(data), id_distributor);
    }

    /*virtual */
    void DestroyVertex(VertexId vertex) {
        VertexId conjugate = vertex->conjugate();
        delete vertex.get();
        delete conjugate.get();
    }

    virtual void AddVertexToGraph(VertexId vertex) {
        this->vertices_.insert(vertex);
        this->vertices_.insert(conjugate(vertex));
    }

    virtual void DeleteVertexFromGraph(VertexId vertex) {
        this->vertices_.erase(vertex);
        this->vertices_.erase(conjugate(vertex));
    }

    virtual void HiddenDeleteVertex(VertexId vertex) {
        DeleteVertexFromGraph(vertex);
        DestroyVertex(vertex);
    }

    virtual VertexId HiddenAddVertex(const VertexData &data, restricted::IdDistributor &id_distributor) {
        VertexId vertex = CreateVertex(data, id_distributor);
        AddVertexToGraph(vertex);
        return vertex;
    }

    virtual EdgeId HiddenAddEdge(const EdgeData &data, restricted::IdDistributor &id_distributor) {
        EdgeId result = AddSingleEdge(VertexId(0), VertexId(0), data, id_distributor);

        if (this->master().isSelfConjugate(data)) {
            result->set_conjugate(result);
            return result;
        }
        EdgeId rcEdge = AddSingleEdge(VertexId(0), VertexId(0),
                                      this->master().conjugate(data),
                                      id_distributor);
        result->set_conjugate(rcEdge);
        rcEdge->set_conjugate(result);
        TRACE("Unlinked edges added");
        return result;
    }

    virtual EdgeId HiddenAddEdge(
            VertexId v1, VertexId v2, const EdgeData &data,
            restricted::IdDistributor &id_distributor) {
        TRACE("Adding edge between vertices " << this->str(v1) << " and " << this->str(v2));

//      todo was suppressed for concurrent execution reasons (see concurrent_graph_component.hpp)
//		VERIFY(this->vertices_.find(v1) != this->vertices_.end() && this->vertices_.find(v2) != this->vertices_.end());

        EdgeId result = AddSingleEdge(v1, v2, data, id_distributor);

        if (this->master().isSelfConjugate(data) && (v1 == conjugate(v2))) {
//              todo why was it removed???
//			Because of some split issues: when self-conjugate edge is split armageddon happends
//          VERIFY(v1 == conjugate(v2));
//          VERIFY(v1 == conjugate(v2));
            result->set_conjugate(result);
            return result;
        }
        EdgeId rcEdge = AddSingleEdge(v2->conjugate(), v1->conjugate(),
                                      this->master().conjugate(data),
                                      id_distributor);
        result->set_conjugate(rcEdge);
        rcEdge->set_conjugate(result);
        TRACE("Edges added");
        return result;
    }

    virtual void HiddenDeleteEdge(EdgeId edge) {
        EdgeId rcEdge = conjugate(edge);
        VertexId rcStart = conjugate(edge->end());
        VertexId start = conjugate(rcEdge->end());
        start->RemoveOutgoingEdge(edge);
        rcStart->RemoveOutgoingEdge(rcEdge);
        if (edge != rcEdge) {
            delete rcEdge.get();
        }
        delete edge.get();
    }

    virtual vector<EdgeId> CorrectMergePath(const vector<EdgeId>& path) const {
        for (size_t i = 0; i < path.size(); i++) {
            if (path[i] == conjugate(path[i])) {
                vector<EdgeId> result;
                if (i < path.size() - 1 - i) {
                    for (size_t j = 0; j < path.size(); j++)
                        result.push_back(conjugate(path[path.size() - 1 - j]));
                    i = path.size() - 1 - i;
                } else {
                    result = path;
                }
                size_t size = 2 * i + 1;
                for (size_t j = result.size(); j < size; j++) {
                    result.push_back(conjugate(result[size - 1 - j]));
                }
                return result;
            }
        }
        return path;
    }

    virtual vector<EdgeId> EdgesToDelete(const vector<EdgeId> &path) const {
        set<EdgeId> edgesToDelete;
        edgesToDelete.insert(path[0]);
        for (size_t i = 0; i + 1 < path.size(); i++) {
            EdgeId e = path[i + 1];
            if (edgesToDelete.find(conjugate(e)) == edgesToDelete.end())
                edgesToDelete.insert(e);
        }
        return vector<EdgeId>(edgesToDelete.begin(), edgesToDelete.end());
    }

    virtual vector<VertexId> VerticesToDelete(
            const vector<EdgeId> &path) const {
        set<VertexId> verticesToDelete;
        for (size_t i = 0; i + 1 < path.size(); i++) {
            EdgeId e = path[i + 1];
            VertexId v = this->EdgeStart(e);
            if (verticesToDelete.find(conjugate(v)) == verticesToDelete.end())
                verticesToDelete.insert(v);
        }
        return vector<VertexId>(verticesToDelete.begin(),
                                verticesToDelete.end());
    }

    EdgeId AddSingleEdge(VertexId v1, VertexId v2, const EdgeData &data,
                         IdDistributor &idDistributor) {
        EdgeId newEdge(new PairedEdge<DataMaster>(v2, data), idDistributor);
        if (v1 != VertexId(0))
            v1->AddOutgoingEdge(newEdge);
        return newEdge;
    }

    virtual void LinkIncomingEdge(VertexId v, EdgeId e) {
        VERIFY(this->EdgeEnd(e) == VertexId(0));
        this->conjugate(v)->AddOutgoingEdge(this->conjugate(e));
        e->SetEndVertex(v);
    }

    virtual void LinkOutgoingEdge(VertexId v, EdgeId e) {
        VERIFY(this->EdgeEnd(this->conjugate(e)) == VertexId(0));
        v->AddOutgoingEdge(e);
        this->conjugate(e)->SetEndVertex(this->conjugate(v));
    }

public:
    AbstractConjugateGraph(const DataMaster& master)
            : base(new PairedHandlerApplier<AbstractConjugateGraph>(*this),
                   master) {
    }

    virtual ~AbstractConjugateGraph() {
        TRACE("~AbstractConjugateGraph")
        for (auto it = this->SmartVertexBegin(); !it.IsEnd(); ++it) {
            this->ForceDeleteVertex(*it);
        }
        TRACE("~AbstractConjugateGraph ok")
    }

    VertexId conjugate(VertexId v) const {
        return v->conjugate();
    }

    EdgeId conjugate(EdgeId edge) const {
        return edge->conjugate();
    }

    HelperT GetConstructionHelper() {
//      TODO: fix everything and restore this check
//      VERIFY(this->VerifyAllDetached());
        return HelperT(*this);
    }

protected:
    /*virtual*/
    bool AdditionalCompressCondition(VertexId v) const {
        return !(this->EdgeEnd(this->GetUniqueOutgoingEdge(v)) == conjugate(v)
                && this->EdgeStart(this->GetUniqueIncomingEdge(v))
                        == conjugate(v));
    }
public:
    /*virtual*/
    bool RelatedVertices(VertexId v1, VertexId v2) const {
        return v1 == v2 || v1 == conjugate(v2);
    }

private:
    DECL_LOGGER("AbstractConjugateGraph")
};

}
