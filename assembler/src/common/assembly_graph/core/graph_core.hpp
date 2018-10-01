//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "id_distributor.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"
#include "order_and_law.hpp"
#include "utils/stl_utils.hpp"

#include "adt/iterator_range.hpp"
#include "adt/small_pod_vector.hpp"

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <btree/safe_btree_set.h>

#include <vector>
#include <set>

namespace omnigraph {

template<class DataMaster>
class GraphCore;

template<class DataMaster>
class ConstructionHelper;

template<class T>
class PairedElementManipulationHelper;

template<class DataMaster>
class PairedVertex;

template<class DataMaster>
class PairedEdge;

namespace impl {
struct Id {
    uint64_t id_;

    Id(uint64_t id = 0)
            : id_(id) {}

    uint64_t int_id() const { return id_; }
    size_t hash() const { return id_; }

    bool operator==(Id other) const { return id_ == other.id_; }
    bool operator!=(Id other) const { return id_ != other.id_; }
    bool operator<(Id rhs) const { return id_ < rhs.id_; }
    bool operator<=(Id rhs) const { return *this < rhs || *this == rhs; }
};

// FIXME: Move to .cpp
inline std::ostream &operator<<(std::ostream &stream, Id id) {
    stream << id.int_id();
    return stream;
}

struct VertexId : public Id {
    using Id::Id;
};

struct EdgeId : public Id {
    using Id::Id;
};

}

template<class It, class Graph>
class conjugate_iterator : public boost::iterator_adaptor<conjugate_iterator<It, Graph>,
                                                          It,
                                                          typename std::iterator_traits<It>::value_type,
                                                          typename std::iterator_traits<It>::iterator_category,
                                                          typename std::iterator_traits<It>::value_type> {
  public:
    typedef typename std::iterator_traits<It>::value_type value_type;

    explicit conjugate_iterator(It it,
                                const Graph *graph,
                                bool conjugate = false)
            : conjugate_iterator::iterator_adaptor(it),
              graph_(graph), conjugate_(conjugate) {}

    conjugate_iterator()
            : graph_(nullptr), conjugate_(false) {}

  private:
    friend class boost::iterator_core_access;

    bool equal(const conjugate_iterator &other) const {
        return this->base() == other.base() && other.conjugate_ == conjugate_;
    }

    value_type dereference() const {
        value_type v = *this->base();
        return (conjugate_ ? graph_->conjugate(v) : v);
    }

    const Graph *graph_;
    bool conjugate_;
};

template<class DataMaster>
class PairedEdge {
private:
    typedef typename DataMaster::EdgeData EdgeData;
    typedef impl::EdgeId EdgeId;
    typedef impl::VertexId VertexId;
    friend class GraphCore<DataMaster>;
    friend class ConstructionHelper<DataMaster>;
    friend class PairedElementManipulationHelper<EdgeId>;
    //todo unfriend
    friend class PairedVertex<DataMaster>;
    VertexId end_;
    EdgeData data_;
    EdgeId conjugate_;

    PairedEdge(VertexId end, const EdgeData &data)
            : end_(end), data_(data) {}

    EdgeData &data() { return data_; }
    const EdgeData &data() const { return data_; }
    void set_data(const EdgeData &data) { data_ = data; }

    VertexId end() const { return end_; }
    void SetEndVertex(VertexId end) { end_ = end; }

    EdgeId conjugate() const { return conjugate_; }
    void set_conjugate(EdgeId conjugate) { conjugate_ = conjugate; }
};

template<class DataMaster>
class PairedVertex {
private:
    typedef typename DataMaster::VertexData VertexData;
    typedef impl::EdgeId EdgeId;
    typedef impl::VertexId VertexId;
    typedef typename adt::SmallPODVector<EdgeId>::const_iterator edge_raw_iterator;

public:
    typedef conjugate_iterator<edge_raw_iterator, GraphCore<DataMaster>> edge_const_iterator;

private:
    friend class GraphCore<DataMaster>;
    friend class ConstructionHelper<DataMaster>;
    friend class PairedEdge<DataMaster>;
    friend class PairedElementManipulationHelper<VertexId>;
    template<class It, class Graph>
    friend class conjugate_iterator;

    adt::SmallPODVector<EdgeId> outgoing_edges_;

    VertexId conjugate_;
    VertexData data_;

    VertexId conjugate() const { return conjugate_; }
    void set_conjugate(VertexId conjugate) { conjugate_ = conjugate; }

    size_t OutgoingEdgeCount() const { return outgoing_edges_.size(); }

    edge_const_iterator out_begin(const GraphCore<DataMaster> *graph, bool conjugate = false) const {
        return edge_const_iterator(outgoing_edges_.cbegin(), graph, conjugate);
    }

    edge_const_iterator out_end(const GraphCore<DataMaster> *graph, bool conjugate = false) const {
        return edge_const_iterator(outgoing_edges_.cend(), graph, conjugate);
    }

    PairedVertex(VertexData data)
            : data_(data) {}

    VertexData &data() { return data_; }
    const VertexData &data() const { return data_; }
    void set_data(VertexData data) { data_ = data; }

    void AddOutgoingEdge(EdgeId e) {
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

    ~PairedVertex() {
        VERIFY(outgoing_edges_.size() == 0);
    }
};

template<class DataMaster>
class GraphCore: private boost::noncopyable {
public:
    typedef DataMaster DataMasterT;
    typedef typename DataMasterT::VertexData VertexData;
    typedef typename DataMasterT::EdgeData EdgeData;
    typedef impl::EdgeId EdgeId;
    typedef impl::VertexId VertexId;
    typedef typename PairedVertex<DataMaster>::edge_const_iterator edge_const_iterator;

private:
    restricted::LocalIdDistributor id_distributor_;
    DataMaster master_;

    friend class ConstructionHelper<DataMaster>;

private:
    class VertexStorage {
      public:
        typedef PairedVertex<DataMaster> Vertex;
        typedef ReclaimingIdDistributor::id_iterator id_iterator;

        VertexStorage(size_t bias)
                : size_(0), id_distributor_(bias) {}

        id_iterator id_begin() const { return id_distributor_.begin(); }
        id_iterator id_end() const { return id_distributor_.end(); }

        void reserve(size_t sz) {
            if (id_distributor_.size() > 2*sz)
                return;

            id_distributor_.resize(2*sz);
            storage_.reserve(2*sz);
        }

        // FIXME: Count!
        size_t size() const { return size_; }

        std::pair<VertexId, VertexId> create(const VertexData& data1, const VertexData& data2) {
            auto v1 = new Vertex(data1), v2 = new Vertex(data2);
            uint64_t vid1 = id_distributor_.allocate(), vid2 = id_distributor_.allocate();

            while (storage_.size() < std::max(vid1, vid2) + 1)
                storage_.resize(storage_.size() * 2 + 1);

            VERIFY(storage_[vid1] == nullptr);
            VERIFY(storage_[vid2] == nullptr);

            storage_[vid1] = v1;
            storage_[vid2] = v2;

            size_ += 2;

            INFO("Create " << vid1 << ":" << vid2);
            return { vid1, vid2 };
        }

        void erase(VertexId id) {
            Vertex *v = storage_[id.int_id()];
            VertexId cid = v->conjugate();
            Vertex *cv = storage_[cid.int_id()];

            INFO("Remove " << id << ":" << cid);

            delete v;
            delete cv;

            id_distributor_.release(id.int_id());
            id_distributor_.release(cid.int_id());
            storage_[id.int_id()] = nullptr;
            storage_[cid.int_id()] = nullptr;
            size_ -= 2;
        }

        Vertex* at(VertexId id) const {
            return storage_.at(id.int_id());
        }
      private:
        size_t size_;
        std::vector<Vertex*> storage_;
        ReclaimingIdDistributor id_distributor_;
    };

    class EdgeStorage {
      public:
        typedef PairedEdge<DataMaster> Edge;

        EdgeStorage(size_t bias)
                : id_distributor_(bias) {}

        void reserve(size_t sz) {
            if (id_distributor_.size() > sz)
                return;

            id_distributor_.resize(sz);
            storage_.reserve(sz);
        }

        EdgeId create(VertexId end, const EdgeData &data) {
            auto e = new Edge(end, data);
            uint64_t id = id_distributor_.allocate();

            while (storage_.size() < id + 1)
                storage_.resize(storage_.size() * 2 + 1);

            VERIFY(storage_[id] == nullptr);
            storage_[id] = e;
            return id;
        }

        void erase(EdgeId id) {
            Edge *e = storage_[id.int_id()];

            delete e;
            storage_[id.int_id()] = nullptr;
            id_distributor_.release(id.int_id());
        }

        Edge* at(EdgeId id) const {
            return storage_.at(id.int_id());
        }
      private:
        std::vector<Edge*> storage_;
        ReclaimingIdDistributor id_distributor_;
    };

    VertexStorage vstorage_;
    EdgeStorage estorage_;

    PairedVertex<DataMaster>* vertex(VertexId id) const {
        return vstorage_.at(id);
    }
    PairedVertex<DataMaster>* cvertex(VertexId id) const {
        return vertex(conjugate(id));
    }

    PairedEdge<DataMaster>* edge(EdgeId id) const {
        return estorage_.at(id.int_id());
    }
    PairedEdge<DataMaster>* cedge(EdgeId id) const {
        return edge(conjugate(id));
    }

public:
    class VertexIt : public boost::iterator_adaptor<VertexIt,
                                                    typename VertexStorage::id_iterator,
                                                    VertexId,
                                                    typename std::iterator_traits<typename VertexStorage::id_iterator>::iterator_category,
                                                    VertexId> {
      public:
        VertexIt(typename VertexStorage::id_iterator it)
                : VertexIt::iterator_adaptor(it) {}

      private:
        friend class boost::iterator_core_access;

        VertexId dereference() const {
            return *this->base();
        }
    };

    VertexIt begin() const { return vstorage_.id_begin(); }
    VertexIt end() const { return vstorage_.id_end(); }
    adt::iterator_range<VertexIt> vertices() const { return { begin(), end()}; }

    size_t size() const { return vstorage_.size(); }

    edge_const_iterator out_begin(VertexId v) const { return vertex(v)->out_begin(this); }
    edge_const_iterator out_end(VertexId v) const { return vertex(v)->out_end(this); }

    edge_const_iterator in_begin(VertexId v) const { return cvertex(v)->out_begin(this, true); }
    edge_const_iterator in_end(VertexId v) const { return cvertex(v)->out_end(this, true); }

private:
    VertexId CreateVertex(const VertexData& data1, const VertexData& data2,
                          restricted::IdDistributor&) {
        VertexId vid1, vid2;
        std::tie(vid1, vid2) = vstorage_.create(data1, data2);

        vertex(vid1)->set_conjugate(vid2);
        vertex(vid2)->set_conjugate(vid1);

        return vid1;
    }

    void DestroyVertex(VertexId v) {
        vstorage_.erase(v);
    }

    EdgeId AddSingleEdge(VertexId v1, VertexId v2, const EdgeData &data,
                         restricted::IdDistributor &) {
        EdgeId eid = estorage_.create(v2, data);
        if (v1.int_id())
            vertex(v1)->AddOutgoingEdge(eid);
        return eid;
    }

    void DestroyEdge(EdgeId e, EdgeId rc) {
        if (e != rc)
            estorage_.erase(rc);
        estorage_.erase(e);
    }

    VertexId CreateVertex(const VertexData &data,
                          restricted::IdDistributor &id_distributor) {
        return CreateVertex(data, master_.conjugate(data), id_distributor);
    }

    VertexId CreateVertex(const VertexData &data) {
        return CreateVertex(data, id_distributor_);
    }

    void AddVertexToGraph(VertexId) {}

    bool AdditionalCompressCondition(VertexId v) const {
        return !(EdgeEnd(GetUniqueOutgoingEdge(v)) == conjugate(v) &&
                 EdgeStart(GetUniqueIncomingEdge(v)) == conjugate(v));
    }

protected:
    VertexId HiddenAddVertex(const VertexData& data, restricted::IdDistributor& id_distributor) {
        return CreateVertex(data, id_distributor);
    }

    VertexId HiddenAddVertex(const VertexData& data) {
        return HiddenAddVertex(data, id_distributor_);
    }

    void HiddenDeleteVertex(VertexId vertex) {
        DestroyVertex(vertex);
    }

    EdgeId HiddenAddEdge(const EdgeData& data,
                         restricted::IdDistributor& id_distributor) {
        EdgeId result = AddSingleEdge(VertexId(), VertexId(), data, id_distributor);
        if (this->master().isSelfConjugate(data)) {
            edge(result)->set_conjugate(result);
            return result;
        }
        EdgeId rcEdge = AddSingleEdge(VertexId(), VertexId(), this->master().conjugate(data), id_distributor);
        edge(result)->set_conjugate(rcEdge);
        edge(rcEdge)->set_conjugate(result);
        return result;
    }

    EdgeId HiddenAddEdge(const EdgeData &data) {
        return HiddenAddEdge(data, id_distributor_);
    }

    EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData& data,
                         restricted::IdDistributor& id_distributor) {
        //      todo was suppressed for concurrent execution reasons (see concurrent_graph_component.hpp)
        //      VERIFY(this->vertices_.find(v1) != this->vertices_.end() && this->vertices_.find(v2) != this->vertices_.end());
        EdgeId result = AddSingleEdge(v1, v2, data, id_distributor);
        if (this->master().isSelfConjugate(data) && (v1 == conjugate(v2))) {
            //              todo why was it removed???
            //          Because of some split issues: when self-conjugate edge is split armageddon happends
            //          VERIFY(v1 == conjugate(v2));
            //          VERIFY(v1 == conjugate(v2));
            edge(result)->set_conjugate(result);
            return result;
        }
        EdgeId rcEdge = AddSingleEdge(vertex(v2)->conjugate(), vertex(v1)->conjugate(), this->master().conjugate(data), id_distributor);
        edge(result)->set_conjugate(rcEdge);
        edge(rcEdge)->set_conjugate(result);
        return result;
    }

    EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData &data) {
        return HiddenAddEdge(v1, v2, data, id_distributor_);
    }

    void HiddenDeleteEdge(EdgeId e) {
        TRACE("Hidden delete edge " << e.int_id());
        EdgeId rcEdge = conjugate(e);
        VertexId rcStart = conjugate(edge(e)->end());
        VertexId start = conjugate(edge(rcEdge)->end());
        vertex(start)->RemoveOutgoingEdge(e);
        vertex(rcStart)->RemoveOutgoingEdge(rcEdge);
        DestroyEdge(e, rcEdge);
    }

    void HiddenDeletePath(const std::vector<EdgeId>& edgesToDelete,
                          const std::vector<VertexId>& verticesToDelete) {
        for (EdgeId e: edgesToDelete)
            HiddenDeleteEdge(e);
        for (VertexId v: verticesToDelete)
            HiddenDeleteVertex(v);
    }

public:
    GraphCore(const DataMaster& master)
            : master_(master),
              vstorage_(/* bias */3), estorage_(/* bias */3) {}

    virtual ~GraphCore() { VERIFY(size() == 0); }

    restricted::LocalIdDistributor &GetGraphIdDistributor() { return id_distributor_; }
    const restricted::LocalIdDistributor &GetGraphIdDistributor() const { return id_distributor_; }

    size_t int_id(EdgeId edge) const { return edge.int_id(); }
    size_t int_id(VertexId vertex) const { return vertex.int_id(); }

    const DataMaster& master() const { return master_; }
    const EdgeData& data(EdgeId e) const { return edge(e)->data(); }
    const VertexData& data(VertexId v) const { return vertex(v)->data(); }
    EdgeData& data(EdgeId e) { return edge(e)->data(); }
    VertexData& data(VertexId v) { return vertex(v)->data(); }

    size_t OutgoingEdgeCount(VertexId v) const { return vertex(v)->OutgoingEdgeCount(); }
    size_t IncomingEdgeCount(VertexId v) const { return cvertex(v)->OutgoingEdgeCount(); }

    adt::iterator_range<edge_const_iterator> OutgoingEdges(VertexId v) const { return { out_begin(v), out_end(v) }; }
    adt::iterator_range<edge_const_iterator> IncomingEdges(VertexId v) const { return { in_begin(v), in_end(v) }; }

    std::vector<EdgeId> GetEdgesBetween(VertexId v, VertexId u) const {
        std::vector<EdgeId> result;
        for (auto e : OutgoingEdges(v)) {
            if (edge(e)->end() != u)
                continue;

            result.push_back(e);
        }
        return result;
    }

    bool RelatedVertices(VertexId v1, VertexId v2) const {
        return v1 == v2 || v1 == conjugate(v2);
    }

    //////////////////////// Edge information
    VertexId EdgeStart(EdgeId edge) const { return conjugate(EdgeEnd(conjugate(edge))); }
    VertexId EdgeEnd(EdgeId e) const { return edge(e)->end(); }

    VertexId conjugate(VertexId v) const { return vertex(v)->conjugate(); }
    EdgeId conjugate(EdgeId e) const { return edge(e)->conjugate(); }

    size_t length(EdgeId edge) const { return master_.length(data(edge)); }
    size_t length(VertexId v) const { return master_.length(data(v)); }

    ////////////////////// shortcut methods
    std::vector<EdgeId> IncidentEdges(VertexId v) const {
        std::vector<EdgeId> answer;
        utils::push_back_all(answer, IncomingEdges(v));
        utils::push_back_all(answer, OutgoingEdges(v));
        return answer;
    }

    EdgeId GetUniqueOutgoingEdge(VertexId v) const {
        VERIFY(CheckUniqueOutgoingEdge(v));
        return *out_begin(v);
    }

    bool CheckUniqueIncomingEdge(VertexId v) const {
        return IncomingEdgeCount(v) == 1;
    }

    EdgeId GetUniqueIncomingEdge(VertexId v) const {
        VERIFY(CheckUniqueIncomingEdge(v));
        return *in_begin(v);
    }

    bool CheckUniqueOutgoingEdge(VertexId v) const {
        return OutgoingEdgeCount(v) == 1;
    }

    bool IsDeadEnd(VertexId v) const { return OutgoingEdgeCount(v) == 0; }
    bool IsDeadStart(VertexId v) const { return IncomingEdgeCount(v) == 0; }

    bool CanCompressVertex(VertexId v) const {
        //      TRACE("Compress vertex check: ");
        //      TRACE("Outgoing check: " << (OutgoingEdgeCount(v) == 1));
        //      TRACE("Outgoing check: " << (CheckUniqueOutgoingEdge(v)));
        //      TRACE("Incoming check: " << (IncomingEdgeCount(v) == 1));
        //      TRACE("Incoming check: " << (CheckUniqueIncomingEdge(v) == 1));
        //      if((OutgoingEdgeCount(v) == 1) && (IncomingEdgeCount(v) == 1)) {
        //          TRACE("Loop check: " << (GetUniqueOutgoingEdge(v) != GetUniqueIncomingEdge(v)));
        //          TRACE("Additional check: " << AdditionalCompressCondition(v));
        //      }
        return OutgoingEdgeCount(v) == 1 && IncomingEdgeCount(v) == 1 &&
            GetUniqueOutgoingEdge(v) != GetUniqueIncomingEdge(v) &&
            AdditionalCompressCondition(v);
    }

    //////////////////////printing
    std::string str(const EdgeId e) const {
//      return master_.str(data(edge));
        std::stringstream ss;
        ss << int_id(e) << " (" << length(e) << ")";
        return ss.str();
    }

    std::string str(const VertexId v) const {
//      return master_.str(data(v));
        return std::to_string(int_id(v));
    }

    std::string detailed_str(const VertexId v) const {
        std::stringstream ss;
        ss << str(v) << ";";
        ss << "Incoming edges " << str(IncomingEdges(v)) << "; ";
        ss << "Outgoing edges " << str(OutgoingEdges(v)) << ";";
        return ss.str();
    }

    std::string detailed_str(const std::vector<EdgeId>& path) const {
        std::stringstream ss;
        ss << "Path: ";
        ss << "Vertex " << detailed_str(EdgeStart(path[0])) << " | ";
        for (auto it = path.begin(); it != path.end(); ++it) {
            EdgeId e = *it;
            ss << "Edge " << str(e) << " | ";
            ss << "Vertex " << detailed_str(EdgeEnd(e)) << " | ";
        }
        return ss.str();
    }

    template<class Container>
    std::string str(const Container& container) const {
        return str(container.begin(), container.end());
    }

    template<class It>
    std::string str(It begin, It end) const {
        std::stringstream ss;
        std::string delim = "";
        for (auto it = begin; it != end; ++it) {
            ss << delim << str(*it);
            delim = ", ";
        }
        return ss.str();
    }

private:
    DECL_LOGGER("GraphCore");
};

}

namespace std {

template<>
struct hash<omnigraph::impl::VertexId> {
    size_t operator()(omnigraph::impl::VertexId v) const {
        return v.hash();
    }
};

template<>
struct hash<omnigraph::impl::EdgeId> {
    size_t operator()(omnigraph::impl::EdgeId e) const {
        return e.hash();
    }
};

}
