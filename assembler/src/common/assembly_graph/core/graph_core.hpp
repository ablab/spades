//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "id_distributor.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"
#include "utils/stl_utils.hpp"

#include "adt/iterator_range.hpp"
#include "adt/small_pod_vector.hpp"

#include <boost/core/noncopyable.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <btree/safe_btree_set.h>

#include <atomic>
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
    explicit operator bool() const { return id_; }

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
    DataMaster master_;

    friend class ConstructionHelper<DataMaster>;

private:
    static constexpr unsigned ID_BIAS = 3;

    template<class T>
    class IdStorage {
      public:
        typedef omnigraph::ReclaimingIdDistributor::id_iterator id_iterator;
        typedef T value_type;

        IdStorage(uint64_t bias = ID_BIAS)
                : size_(0), bias_(bias), id_distributor_(bias) {
            storage_.resize(id_distributor_.size() + bias_);
        }

        id_iterator id_begin() const { return id_distributor_.begin(); }
        id_iterator id_end() const { return id_distributor_.end(); }

        void reserve(size_t sz) {
            if (storage_.size() >= sz + bias_)
                return;

            id_distributor_.resize(sz);
            storage_.resize(sz + bias_);
        }

        // FIXME: Count!
        size_t size() const { return size_; }

        bool contains(uint64_t id) const {
            return (id < storage_.size()) && storage_[id];
        }

        template<typename... ArgTypes>
        uint64_t create(ArgTypes &&... args) {
            uint64_t id = id_distributor_.allocate();

            while (storage_.size() < id + 1)
                storage_.resize(storage_.size() * 2 + 1);

            VERIFY(storage_[id] == nullptr);
            storage_[id] = new T(std::forward<ArgTypes>(args)...);;
            size_ += 1;

            // INFO("Create " << vid1 << ":" << vid2);
            return id;
        }

        template<typename... ArgTypes>
        uint64_t emplace(uint64_t at, ArgTypes &&... args) {
            // One MUST call reserve before using emplace()
            VERIFY(storage_.at(at) == nullptr);
            VERIFY(!id_distributor_.occupied(at));

            id_distributor_.acquire(at);
            storage_[at] = new T(std::forward<ArgTypes>(args)...);;
            size_.fetch_add(1);

            // INFO("Create " << vid1 << ":" << vid2);
            return at;
        }

        void erase(uint64_t id) {
            auto *v = storage_[id];

            // INFO("Remove " << id << ":" << cid);
            delete v;

            id_distributor_.release(id);
            storage_[id] = nullptr;
            size_ -= 1;
        }

        T* at(uint64_t id) const {
            return storage_[id];
        }

        uint64_t reserved() const { return id_distributor_.size(); }

      private:
        std::atomic<size_t> size_;
        uint64_t bias_;
        std::vector<T*> storage_;
        omnigraph::ReclaimingIdDistributor id_distributor_;
    };

    using VertexStorage = IdStorage<PairedVertex<DataMaster>>;
    VertexStorage vstorage_;
    using EdgeStorage = IdStorage<PairedEdge<DataMaster>>;
    EdgeStorage estorage_;

    PairedVertex<DataMaster>* vertex(VertexId id) const {
        return vstorage_.at(id.int_id());
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
    size_t e_size() const { return estorage_.size(); }

    class EdgeIt : public boost::iterator_adaptor<EdgeIt,
                                                  typename EdgeStorage::id_iterator,
                                                  EdgeId,
                                                  typename std::iterator_traits<typename EdgeStorage::id_iterator>::iterator_category,
                                                  EdgeId> {
      public:
        EdgeIt(typename EdgeStorage::id_iterator it)
                : EdgeIt::iterator_adaptor(it) {}

      private:
        friend class boost::iterator_core_access;

        EdgeId dereference() const {
            return *this->base();
        }
    };

    EdgeIt e_begin() const { return estorage_.id_begin(); }
    EdgeIt e_end() const { return estorage_.id_end(); }
    adt::iterator_range<EdgeIt> edges() const { return { e_begin(), e_end()}; }

    edge_const_iterator out_begin(VertexId v) const { return vertex(v)->out_begin(this); }
    edge_const_iterator out_end(VertexId v) const { return vertex(v)->out_end(this); }

    edge_const_iterator in_begin(VertexId v) const { return cvertex(v)->out_begin(this, true); }
    edge_const_iterator in_end(VertexId v) const { return cvertex(v)->out_end(this, true); }

private:
    VertexId CreateVertex(const VertexData &data, VertexId id1 = 0, VertexId id2 = 0) {
        return CreateVertex(data, master_.conjugate(data), id1, id2);
    }

    VertexId CreateVertex(const VertexData& data1, const VertexData& data2,
                          VertexId id1 = 0, VertexId id2 = 0) {
        if (id1 && !id2)
            id2 = id1.int_id() + 1;

        VertexId vid1 = (id1 ? vstorage_.emplace(id1.int_id(), data1) : vstorage_.create(data1));
        VertexId vid2 = (id2 ? vstorage_.emplace(id2.int_id(), data2) : vstorage_.create(data2));

        vertex(vid1)->set_conjugate(vid2);
        vertex(vid2)->set_conjugate(vid1);

        return vid1;
    }

    void DestroyVertex(VertexId v) {
        VertexId cv = conjugate(v);
        vstorage_.erase(v.int_id());
        vstorage_.erase(cv.int_id());
    }

    EdgeId AddSingleEdge(VertexId v1, VertexId v2,
                         const EdgeData &data, EdgeId id = 0) {
        EdgeId eid = (id ?
                      estorage_.emplace(id.int_id(), v2, data) :
                      estorage_.create(v2, data));
        if (v1.int_id())
            vertex(v1)->AddOutgoingEdge(eid);
        return eid;
    }

    void DestroyEdge(EdgeId e, EdgeId rc) {
        if (e != rc)
            estorage_.erase(rc.int_id());
        estorage_.erase(e.int_id());
    }

    bool AdditionalCompressCondition(VertexId v) const {
        return !(EdgeEnd(GetUniqueOutgoingEdge(v)) == conjugate(v) &&
                 EdgeStart(GetUniqueIncomingEdge(v)) == conjugate(v));
    }

protected:
    VertexId HiddenAddVertex(const VertexData& data, VertexId at1 = 0, VertexId at2 = 0) {
        return CreateVertex(data, at1, at2);
    }

    void HiddenDeleteVertex(VertexId vertex) {
        DestroyVertex(vertex);
    }

    EdgeId HiddenAddEdge(const EdgeData& data,
                         EdgeId at1 = 0, EdgeId at2 = 0) {
        EdgeId result = AddSingleEdge(VertexId(), VertexId(), data, at1);
        if (this->master().isSelfConjugate(data)) {
            edge(result)->set_conjugate(result);
            return result;
        }

        if (at1 && !at2)
            at2 = at1.int_id() + 1;
        EdgeId rcEdge = AddSingleEdge(VertexId(), VertexId(), this->master().conjugate(data),
                                      at2);
        edge(result)->set_conjugate(rcEdge);
        edge(rcEdge)->set_conjugate(result);
        return result;
    }

    EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData& data,
                         EdgeId at1 = 0, EdgeId at2 = 0) {
        //      todo was suppressed for concurrent execution reasons (see concurrent_graph_component.hpp)
        //      VERIFY(this->vertices_.find(v1) != this->vertices_.end() && this->vertices_.find(v2) != this->vertices_.end());
        EdgeId result = AddSingleEdge(v1, v2, data, at1);
        if (this->master().isSelfConjugate(data) && (v1 == conjugate(v2))) {
            //              todo why was it removed???
            //          Because of some split issues: when self-conjugate edge is split armageddon happends
            //          VERIFY(v1 == conjugate(v2));
            //          VERIFY(v1 == conjugate(v2));
            edge(result)->set_conjugate(result);
            return result;
        }

        if (at1 && !at2)
            at2 = at1.int_id() + 1;
        EdgeId rcEdge = AddSingleEdge(vertex(v2)->conjugate(), vertex(v1)->conjugate(),
                                      this->master().conjugate(data), at2);
        edge(result)->set_conjugate(rcEdge);
        edge(rcEdge)->set_conjugate(result);
        return result;
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
              vstorage_(ID_BIAS), estorage_(ID_BIAS) {}

    virtual ~GraphCore() { VERIFY(size() == 0); }

    void vreserve(size_t sz) { vstorage_.reserve(sz); }
    void ereserve(size_t sz) { estorage_.reserve(sz); }
    void reserve(size_t vertices, size_t edges) {
        vreserve(vertices);
        ereserve(edges);
    }

    size_t vreserved() const { return vstorage_.reserved(); }
    size_t ereserved() const { return estorage_.reserved(); }

    uint64_t min_id() const { return ID_BIAS; }

    bool contains(VertexId vertex) const {
        return vstorage_.contains(vertex.int_id());
    }

    bool contains(EdgeId edge) const {
        return estorage_.contains(edge.int_id());
    }

    size_t int_id(EdgeId edge) const { return edge.int_id(); }
    size_t int_id(VertexId vertex) const { return vertex.int_id(); }

    const DataMaster& master() const { return master_; }
    const EdgeData& data(EdgeId e) const { return edge(e)->data(); }
    const VertexData& data(VertexId v) const { return vertex(v)->data(); }
    EdgeData& data(EdgeId e) { return edge(e)->data(); }
    VertexData& data(VertexId v) { return vertex(v)->data(); }

    size_t OutgoingEdgeCount(VertexId v) const { return vertex(v)->OutgoingEdgeCount(); }
    size_t IncomingEdgeCount(VertexId v) const { return cvertex(v)->OutgoingEdgeCount(); }

    adt::iterator_range<edge_const_iterator> OutgoingEdges(VertexId v) const {
        auto vertex = this->vertex(v);
        return { vertex->out_begin(this), vertex->out_end(this) };
    }

    adt::iterator_range<edge_const_iterator> IncomingEdges(VertexId v) const {
        auto vertex = this->cvertex(v);
        return { vertex->out_begin(this, true), vertex->out_end(this, true) };
    }

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
    bool CheckUniqueOutgoingEdge(VertexId v) const {
        return OutgoingEdgeCount(v) == 1;
    }

    bool CheckUniqueIncomingEdge(VertexId v) const {
        return IncomingEdgeCount(v) == 1;
    }
    EdgeId GetUniqueIncomingEdge(VertexId v) const {
        VERIFY(CheckUniqueIncomingEdge(v));
        return *in_begin(v);
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
    std::string str(EdgeId e) const {
//      return master_.str(data(edge));
        std::stringstream ss;
        ss << int_id(e) << " (" << length(e) << ")";
        return ss.str();
    }

    std::string str(VertexId v) const {
//      return master_.str(data(v));
        return std::to_string(int_id(v));
    }

    std::string detailed_str(VertexId v) const {
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
