//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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

template<class DataMaster>
class PairedVertex {
 private:
    typedef restricted::pure_pointer<PairedVertex<DataMaster>> VertexId;
    typedef restricted::pure_pointer<PairedEdge<DataMaster>> EdgeId;
    typedef typename DataMaster::VertexData VertexData;
    typedef typename std::vector<EdgeId>::const_iterator edge_raw_iterator;

    class conjugate_iterator : public std::iterator<std::forward_iterator_tag,
            typename std::iterator_traits<edge_raw_iterator>::value_type> {
        edge_raw_iterator it_;
        const bool conjugate_;

     public:
        typedef typename std::iterator_traits<edge_raw_iterator>::difference_type difference_type;
        typedef typename std::iterator_traits<edge_raw_iterator>::pointer pointer;
        typedef typename std::iterator_traits<edge_raw_iterator>::value_type reference;
        typedef typename std::iterator_traits<edge_raw_iterator>::value_type value_type;

        explicit conjugate_iterator(edge_raw_iterator it,
                                    bool conjugate = false)
                : it_(it),
                  conjugate_(conjugate) {
        }

        // Should not exist. Temporary patch to write empty in_begin, out_end... methods for
        // ConcurrentGraphComponent which can not have such methods by definition
        conjugate_iterator()
                : conjugate_(false) {
            VERIFY_MSG(false, "There is no sense in using this. See comments.")
        }

        reference operator*() const {
            if (conjugate_)
                return (*it_)->conjugate();
            else
                return *it_;
        }

        reference operator[](difference_type n) const {
            return *(*this + n);
        }

        conjugate_iterator& operator++() {
            ++it_;
            return *this;
        }

        conjugate_iterator& operator--() {
            --it_;
            return *this;
        }

        conjugate_iterator operator++(int) {
            conjugate_iterator res = *this;
            ++it_;
            return res;
        }

        conjugate_iterator operator--(int) {
            conjugate_iterator res = *this;
            --it_;
            return res;
        }

        conjugate_iterator& operator+(const difference_type &n) {
            return conjugate_iterator(it_ + n, conjugate_);
        }

        conjugate_iterator& operator-(const difference_type &n) {
            return conjugate_iterator(it_ - n, conjugate_);
        }

        conjugate_iterator& operator+=(const difference_type &n) {
            it_ += n;
            return *this;
        }

        conjugate_iterator& operator-=(const difference_type &n) {
            it_ -= n;
            return *this;
        }

        friend bool operator==(const conjugate_iterator &r1,
                               const conjugate_iterator &r2) {
            return r1.it_ == r2.it_;
        }

        friend bool operator!=(const conjugate_iterator &r1,
                               const conjugate_iterator &r2) {
            return r1.it_ != r2.it_;
        }

        friend bool operator<(const conjugate_iterator &r1,
                              const conjugate_iterator &r2) {
            return r1.it_ < r2.it_;
        }

        friend bool operator<=(const conjugate_iterator &r1,
                               const conjugate_iterator &r2) {
            return r1.it_ <= r2.it_;
        }

        friend bool operator>=(const conjugate_iterator &r1,
                               const conjugate_iterator &r2) {
            return r1.it_ >= r2.it_;
        }

        friend bool operator>(const conjugate_iterator &r1,
                              const conjugate_iterator &r2) {
            return r1.it_ > r2.it_;
        }

        friend conjugate_iterator operator+(difference_type n,
                                            const conjugate_iterator &i1) {
            return i1 + n;
        }

        friend difference_type operator-(const conjugate_iterator &i1,
                                         const conjugate_iterator &i2) {
            return i1.it_ - i2.it_;
        }
    };

 public:
    typedef conjugate_iterator edge_const_iterator;

 private:
    friend class AbstractGraph<
            restricted::pure_pointer<PairedVertex<DataMaster>>,
            restricted::pure_pointer<PairedEdge<DataMaster>>, DataMaster> ;
    friend class AbstractConjugateGraph<DataMaster> ;
    friend class PairedEdge<DataMaster> ;
    friend class conjugate_iterator;

    std::vector<EdgeId> outgoing_edges_;

    VertexId conjugate_;

    VertexData data_;

    void set_conjugate(VertexId conjugate) {
        conjugate_ = conjugate;
    }

    size_t OutgoingEdgeCount() const {
        return outgoing_edges_.size();
    }

    const vector<EdgeId>& OutgoingEdges() const {
        return outgoing_edges_;
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

    const vector<EdgeId> IncomingEdges() const {
        vector<EdgeId> result = conjugate_->OutgoingEdges();
        for (size_t i = 0; i < result.size(); i++) {
            result[i] = result[i]->conjugate();
        }
        return result;
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
        outgoing_edges_.push_back(e);
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

template<class DataMaster>
class PairedEdge : public CoveredEdge {
 private:
    typedef restricted::pure_pointer<PairedVertex<DataMaster>> VertexId;
    typedef restricted::pure_pointer<PairedEdge<DataMaster>> EdgeId;
    typedef typename DataMaster::EdgeData EdgeData;
    friend class AbstractGraph<
            restricted::pure_pointer<PairedVertex<DataMaster>>,
            restricted::pure_pointer<PairedEdge<DataMaster>>, DataMaster> ;
    friend class AbstractConjugateGraph<DataMaster> ;
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
    EdgeId conjugate() {
        return conjugate_;
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

 public:

    typedef typename base::VertexId VertexId;
    typedef typename base::EdgeId EdgeId;
    typedef typename base::VertexData VertexData;
    typedef typename base::EdgeData EdgeData;
    typedef typename base::VertexIterator VertexIterator;

 private:

    VertexId CreateVertex(const VertexData &data1, const VertexData &data2) {
        VertexId vertex1(new PairedVertex<DataMaster>(data1));
        VertexId vertex2(new PairedVertex<DataMaster>(data2));

        vertex1->set_conjugate(vertex2);
        vertex2->set_conjugate(vertex1);

        return vertex1;
    }

    virtual VertexId CreateVertex(const VertexData &data) {
        return CreateVertex(data, this->master().conjugate(data));
    }

    /*virtual */void DestroyVertex(VertexId vertex) {
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

    virtual VertexId HiddenAddVertex(const VertexData &data) {
        VertexId vertex = CreateVertex(data);
        AddVertexToGraph(vertex);
        return vertex;
    }

//    virtual EdgeId HiddenAddEdge(const EdgeData &data) {
//        return HiddenAddEdge(data, restricted::GlobalIdDistributor::GetInstance());
//    }

//    virtual EdgeId HiddenAddEdge(VertexId v1, VertexId v2,
//                                 const EdgeData &data) {
//        return HiddenAddEdge(v1, v2, data,
//                             restricted::GlobalIdDistributor::GetInstance());
//    }

    virtual EdgeId HiddenAddEdge(const EdgeData &data,
                                 restricted::IdDistributor * id_distributor = restricted::GlobalIdDistributor::GetInstance()) {
        EdgeId result = AddSingleEdge(VertexId(0), VertexId(0), data, id_distributor);

        if (this->master().isSelfConjugate(data)) {

            result->set_conjugate(result);
            return result;
        }
        EdgeId rcEdge = AddSingleEdge(VertexId(0), VertexId(0), this->master().conjugate(data), id_distributor);
        result->set_conjugate(rcEdge);
        rcEdge->set_conjugate(result);
        TRACE("Unlinked edges added");
        return result;
    }

    virtual EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData &data,
                                 restricted::IdDistributor * id_distributor = restricted::GlobalIdDistributor::GetInstance()) {
        TRACE("Adding edge between vertices " << this->str(v1) << " and " << this->str(v2));

//      todo was suppressed for concurrent execution reasons (see concurrent_graph_component.hpp)
//		VERIFY(this->vertices_.find(v1) != this->vertices_.end() && this->vertices_.find(v2) != this->vertices_.end());

        EdgeId result = AddSingleEdge(v1, v2, data, id_distributor);

        if (this->master().isSelfConjugate(data) && (v1 == conjugate(v2))) {
//      	todo why was it removed???
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
                         IdDistributor * idDistributor) {
        EdgeId newEdge(new PairedEdge<DataMaster>(v2, data), idDistributor);
        if(v1 != VertexId(0))
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

    virtual bool SplitCondition(VertexId vertex,
                                const vector<EdgeId> &splitting_dges) const {
        FOREACH (EdgeId e, splitting_dges) {
            if (this->EdgeStart(e) == conjugate(this->EdgeEnd(e)))
                return false;
        }
        return true;
    }

 private:
    DECL_LOGGER("AbstractConjugateGraph")
};

}
