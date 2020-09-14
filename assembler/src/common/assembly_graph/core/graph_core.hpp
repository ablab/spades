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

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/noncopyable.hpp>
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

    uint64_t int_id() const noexcept { return id_; }
    size_t hash() const noexcept { return id_; }
    explicit operator bool() const noexcept { return id_; }

    bool operator==(Id other) const noexcept { return id_ == other.id_; }
    bool operator!=(Id other) const noexcept { return id_ != other.id_; }
    bool operator<(Id rhs) const noexcept { return id_ < rhs.id_; }
    bool operator<=(Id rhs) const noexcept { return *this < rhs || *this == rhs; }

    template <typename Archive>
    void BinArchive(Archive &ar) {
        ar(id_);
    }
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
    EdgeId conjugate_;
    EdgeData data_;

    PairedEdge(VertexId end, const EdgeData &data)
            : end_(end), data_(data) {}

    PairedEdge(PairedEdge &&that) = default;

    EdgeData &data() noexcept { return data_; }
    const EdgeData &data() const noexcept { return data_; }
    void set_data(const EdgeData &data)  noexcept { data_ = data; }

    VertexId end() const noexcept { return end_; }
    void SetEndVertex(VertexId end) noexcept { end_ = end; }

    EdgeId conjugate() const noexcept { return conjugate_; }
    void set_conjugate(EdgeId conjugate) noexcept { conjugate_ = conjugate; }
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

    VertexId conjugate() const noexcept { return conjugate_; }
    void set_conjugate(VertexId conjugate) noexcept { conjugate_ = conjugate; }

    size_t OutgoingEdgeCount() const noexcept { return outgoing_edges_.size(); }

    edge_const_iterator out_begin(const GraphCore<DataMaster> *graph, bool conjugate = false) const {
        return edge_const_iterator(outgoing_edges_.cbegin(), graph, conjugate);
    }

    edge_const_iterator out_end(const GraphCore<DataMaster> *graph, bool conjugate = false) const {
        return edge_const_iterator(outgoing_edges_.cend(), graph, conjugate);
    }

    PairedVertex(const VertexData &data)
            : data_(data) {}

    PairedVertex(PairedVertex &&that) = default;

    VertexData &data() noexcept { return data_; }
    const VertexData &data() const noexcept { return data_; }
    void set_data(VertexData data) noexcept { data_ = data; }

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
      private:
        void resize(size_t N) {
            VERIFY(N > storage_size_);
            T *new_storage = (T*)malloc(N * sizeof(T));
            for (uint64_t id = bias_; id < storage_size_; ++id) {
                if (!id_distributor_.occupied(id))
                    continue;

                T *val = &storage_[id];
                ::new(new_storage + id) T(std::move(*val));
                val->~T();
            }

            if (storage_)
                free(storage_);
            storage_ = new_storage;
            storage_size_ = N;
        }

      public:
        typedef omnigraph::ReclaimingIdDistributor::id_iterator id_iterator;
        typedef T value_type;

        IdStorage(uint64_t bias = ID_BIAS)
                : size_(0), bias_(bias), storage_(nullptr), storage_size_(0), id_distributor_(bias) {
            resize(id_distributor_.size() + bias_);
        }

        ~IdStorage() {
            if (storage_)
                free(storage_);
        }

        id_iterator id_begin() const { return id_distributor_.begin(); }
        id_iterator id_end() const { return id_distributor_.end(); }
        uint64_t max_id() const { return id_distributor_.max_id(); }

        void reserve(size_t sz) {
            if (storage_size_ >= sz + bias_)
                return;

            id_distributor_.resize(sz);
            resize(sz + bias_);
        }

        // FIXME: Count!
        size_t size() const noexcept { return size_; }

        bool contains(uint64_t id) const {
            return id < storage_size_ && id_distributor_.occupied(id);
        }

        template<typename... ArgTypes>
        uint64_t create(ArgTypes &&... args) {
            uint64_t id = id_distributor_.allocate();

            while (storage_size_ < id + 1)
                resize(storage_size_ * 2 + 1);

            new(storage_ + id) T(std::forward<ArgTypes>(args)...);;
            size_ += 1;

            // INFO("Create " << id);
            return id;
        }

        template<typename... ArgTypes>
        uint64_t emplace(uint64_t at, ArgTypes &&... args) {
            // One MUST call reserve before using emplace()
            VERIFY(!id_distributor_.occupied(at));

            id_distributor_.acquire(at);
            new(storage_ + at) T(std::forward<ArgTypes>(args)...);;
            size_.fetch_add(1);

            // INFO("Emplace " << at);
            return at;
        }

        void erase(uint64_t id) {
            T *v = &storage_[id];

            // INFO("Remove " << id);
            v->~T();

            id_distributor_.release(id);
            size_ -= 1;
        }

        T& at(uint64_t id) const noexcept {
            return storage_[id];
        }

        uint64_t reserved() const { return id_distributor_.size(); }
        void clear_state() { id_distributor_.clear_state(); }

      private:
        std::atomic<size_t> size_;
        uint64_t bias_;
        T *storage_;
        size_t storage_size_;
        omnigraph::ReclaimingIdDistributor id_distributor_;
    };

    using VertexStorage = IdStorage<PairedVertex<DataMaster>>;
    VertexStorage vstorage_;
    using EdgeStorage = IdStorage<PairedEdge<DataMaster>>;
    EdgeStorage estorage_;

    PairedVertex<DataMaster>& vertex(VertexId id) const noexcept {
        return vstorage_.at(id.int_id());
    }
    PairedVertex<DataMaster>& cvertex(VertexId id) const noexcept {
        return vertex(conjugate(id));
    }

    PairedEdge<DataMaster>& edge(EdgeId id) const noexcept {
        return estorage_.at(id.int_id());
    }
    PairedEdge<DataMaster>& cedge(EdgeId id) const noexcept {
        return edge(conjugate(id));
    }

    struct EdgePredicate {
        typedef EdgeId operand_type;
        
        EdgePredicate(const GraphCore<DataMaster> &graph)
                : graph_(graph) {}

        bool operator()(EdgeId) const noexcept;

        std::reference_wrapper<const GraphCore<DataMaster>> graph_;
    };

    struct AllEdges : public EdgePredicate {
        using EdgePredicate::EdgePredicate;
        using typename EdgePredicate::operand_type;

        bool operator()(EdgeId) const noexcept {
            return true;
        }
    };

    struct CanonicalEdges : public EdgePredicate {
        using EdgePredicate::EdgePredicate;
        using typename EdgePredicate::operand_type;

        bool operator()(EdgeId e) const noexcept {
            const GraphCore<DataMaster> &g = this->graph_;
            return e <= g.conjugate(e);
        }
    };

    struct VertexPredicate {
        typedef VertexId operand_type;
        
        VertexPredicate(const GraphCore<DataMaster> &graph)
                : graph_(graph) {}

        bool operator()(VertexId) const;

        std::reference_wrapper<const GraphCore<DataMaster>> graph_;
    };

    struct AllVertices : public VertexPredicate {
        using VertexPredicate::VertexPredicate;
        using typename VertexPredicate::operand_type;

        bool operator()(VertexId) const {
            return true;
        }
    };

    struct CanonicalVertices : public VertexPredicate {
        using VertexPredicate::VertexPredicate;
        using typename VertexPredicate::operand_type;

        bool operator()(VertexId e) const {
            const GraphCore<DataMaster> &g = this->graph_;
            return e <= g.conjugate(e);
        }
    };
    
    template<class Predicate1, class Predicate2>
    struct And {
        And(Predicate1 p1, Predicate2 p2)
                : p1_(std::move(p1)), p2_(std::move(p2)) {}

        bool operator()(typename Predicate1::operand_type e) const noexcept {
            return p1_(e) && p2_(e);
        }

        Predicate1 p1_;
        Predicate2 p2_;
    };

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

    size_t size() const noexcept { return vstorage_.size(); }
    size_t e_size() const noexcept { return estorage_.size(); }
    size_t max_eid() const { return estorage_.max_id(); }
    size_t max_vid() const { return vstorage_.max_id(); }

    template<class Predicate, bool Canonical = false>
    auto begin(Predicate p) const {
        using BasePredicate = typename std::conditional<Canonical, CanonicalVertices, AllVertices>::type;
        return boost::make_filter_iterator(And<BasePredicate, Predicate>(BasePredicate(*this), std::move(p)),
                                           vstorage_.id_begin(), vstorage_.id_end());
    }
    template<class Predicate, bool Canonical = false>
    auto end(Predicate p) const {
        using BasePredicate = typename std::conditional<Canonical, CanonicalVertices, AllVertices>::type;
        return boost::make_filter_iterator(And<BasePredicate, Predicate>(BasePredicate(*this), std::move(p)),
                                           vstorage_.id_end(), vstorage_.id_end());
    }
    template<class Predicate, bool Canonical = false>
    auto vertices(Predicate p) const {
        return adt::make_range(begin<Predicate, Canonical>(std::move(p)),
                               end<Predicate, Canonical>(std::move(p)));
    }
    template<class Predicate>
    auto canonical_vertices(Predicate p) const {
        return vertices<Predicate, true>(std::move(p));
    }

    template<bool Canonical>
    auto begin(std::enable_if_t<Canonical, int> = 0) const {
        return boost::make_filter_iterator(CanonicalVertices(*this),
                                           vstorage_.id_begin(), vstorage_.id_end());
    }
    template<bool Canonical>
    auto end(std::enable_if_t<Canonical, int> = 0) const {
        return boost::make_filter_iterator(CanonicalVertices(*this),
                                           vstorage_.id_end(), vstorage_.id_end());
    }
    template<bool Canonical = false>
    VertexIt begin(std::enable_if_t<!Canonical, int> = 0) const {
        return vstorage_.id_begin();
    }
    template<bool Canonical = false>
    VertexIt end(std::enable_if_t<!Canonical, int> = 0) const {
        return vstorage_.id_end();
    }
    template<bool Canonical = false>
    auto vertices() const {
        return adt::make_range(begin<Canonical>(), end<Canonical>());
    }
    auto canonical_vertices() const {
        return vertices<true>();
    }
    
    template<class Predicate, bool Canonical = false>
    auto e_begin(Predicate p) const {
        using BasePredicate = typename std::conditional<Canonical, CanonicalEdges, AllEdges>::type;
        return boost::make_filter_iterator(And<BasePredicate, Predicate>(BasePredicate(*this), std::move(p)),
                                           estorage_.id_begin(), estorage_.id_end());
    }
    template<class Predicate, bool Canonical = false>
    auto e_end(Predicate p) const {
        using BasePredicate = typename std::conditional<Canonical, CanonicalEdges, AllEdges>::type;
        return boost::make_filter_iterator(And<BasePredicate, Predicate>(BasePredicate(*this), std::move(p)),
                                           estorage_.id_end(), estorage_.id_end());
    }
    template<class Predicate, bool Canonical = false>
    auto edges(Predicate p) const {
        return adt::make_range(e_begin<Predicate, Canonical>(std::move(p)),
                               e_end<Predicate, Canonical>(std::move(p)));
    }
    template<class Predicate>
    auto canonical_edges(Predicate p) const {
        return edges<Predicate, true>(std::move(p));
    }

    template<bool Canonical>
    auto e_begin(std::enable_if_t<Canonical, int> = 0) const {
        return boost::make_filter_iterator(CanonicalEdges(*this),
                                           estorage_.id_begin(), estorage_.id_end());
    }
    template<bool Canonical>
    auto e_end(std::enable_if_t<Canonical, int> = 0) const {
        return boost::make_filter_iterator(CanonicalEdges(*this),
                                           estorage_.id_end(), estorage_.id_end());
    }
    template<bool Canonical = false>
    auto e_begin(std::enable_if_t<!Canonical, int> = 0) const {
        return estorage_.id_begin();
    }
    template<bool Canonical = false>
    auto e_end(std::enable_if_t<!Canonical, int> = 0) const {
        return estorage_.id_end();
    }
    template<bool Canonical = false>
    auto edges() const {
        return adt::make_range(e_begin<Canonical>(), e_end<Canonical>());
    }
    auto canonical_edges() const {
        return edges<true>();
    }

    edge_const_iterator out_begin(VertexId v) const { return vertex(v).out_begin(this); }
    edge_const_iterator out_end(VertexId v) const { return vertex(v).out_end(this); }

    edge_const_iterator in_begin(VertexId v) const { return cvertex(v).out_begin(this, true); }
    edge_const_iterator in_end(VertexId v) const { return cvertex(v).out_end(this, true); }

    void clear_state() { estorage_.clear_state(); vstorage_.clear_state(); }

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

        vertex(vid1).set_conjugate(vid2);
        vertex(vid2).set_conjugate(vid1);

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
            vertex(v1).AddOutgoingEdge(eid);
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
            edge(result).set_conjugate(result);
            return result;
        }

        if (at1 && !at2)
            at2 = at1.int_id() + 1;
        EdgeId rcEdge = AddSingleEdge(VertexId(), VertexId(), this->master().conjugate(data),
                                      at2);
        edge(result).set_conjugate(rcEdge);
        edge(rcEdge).set_conjugate(result);
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
            edge(result).set_conjugate(result);
            return result;
        }

        if (at1 && !at2)
            at2 = at1.int_id() + 1;
        EdgeId rcEdge = AddSingleEdge(vertex(v2).conjugate(), vertex(v1).conjugate(),
                                      this->master().conjugate(data), at2);
        edge(result).set_conjugate(rcEdge);
        edge(rcEdge).set_conjugate(result);
        return result;
    }

    void HiddenDeleteEdge(EdgeId e) {
        TRACE("Hidden delete edge " << e.int_id());
        EdgeId rcEdge = conjugate(e);
        VertexId rcStart = conjugate(edge(e).end());
        VertexId start = conjugate(edge(rcEdge).end());
        vertex(start).RemoveOutgoingEdge(e);
        vertex(rcStart).RemoveOutgoingEdge(rcEdge);
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
              vstorage_(ID_BIAS), estorage_(ID_BIAS) {
        INFO("Graph created, vertex min_id: " << ID_BIAS << ", edge min_id: " << ID_BIAS);
        INFO("Vertex size: " << sizeof(PairedVertex<DataMaster>) << ", edge size: " << sizeof(PairedEdge<DataMaster>));
    }

    virtual ~GraphCore() { VERIFY(size() == 0); }

    void vreserve(size_t sz) { vstorage_.reserve(sz); }
    void ereserve(size_t sz) { estorage_.reserve(sz); }
    void reserve(size_t vertices, size_t edges) {
        vreserve(vertices);
        ereserve(edges);
    }

    size_t vreserved() const { return vstorage_.reserved(); }
    size_t ereserved() const { return estorage_.reserved(); }

    uint64_t min_id() const noexcept { return ID_BIAS; }

    bool contains(VertexId vertex) const {
        return vstorage_.contains(vertex.int_id());
    }

    bool contains(EdgeId edge) const {
        return estorage_.contains(edge.int_id());
    }

    size_t int_id(EdgeId edge) const noexcept { return edge.int_id(); }
    size_t int_id(VertexId vertex) const noexcept { return vertex.int_id(); }

    const DataMaster& master() const noexcept { return master_; }
    const EdgeData& data(EdgeId e) const noexcept { return edge(e).data(); }
    const VertexData& data(VertexId v) const noexcept { return vertex(v).data(); }
    EdgeData& data(EdgeId e) noexcept { return edge(e).data(); }
    VertexData& data(VertexId v) noexcept { return vertex(v).data(); }

    size_t OutgoingEdgeCount(VertexId v) const noexcept { return vertex(v).OutgoingEdgeCount(); }
    size_t IncomingEdgeCount(VertexId v) const noexcept { return cvertex(v).OutgoingEdgeCount(); }

    adt::iterator_range<edge_const_iterator> OutgoingEdges(VertexId v) const {
        auto &vertex = this->vertex(v);
        return { vertex.out_begin(this), vertex.out_end(this) };
    }

    adt::iterator_range<edge_const_iterator> IncomingEdges(VertexId v) const {
        auto &vertex = this->cvertex(v);
        return { vertex.out_begin(this, true), vertex.out_end(this, true) };
    }

    std::vector<EdgeId> GetEdgesBetween(VertexId v, VertexId u) const {
        std::vector<EdgeId> result;
        for (auto e : OutgoingEdges(v)) {
            if (edge(e).end() != u)
                continue;

            result.push_back(e);
        }
        return result;
    }

    bool RelatedVertices(VertexId v1, VertexId v2) const noexcept {
        return v1 == v2 || v1 == conjugate(v2);
    }

    //////////////////////// Edge information
    VertexId EdgeStart(EdgeId edge) const noexcept { return conjugate(EdgeEnd(conjugate(edge))); }
    VertexId EdgeEnd(EdgeId e) const noexcept { return edge(e).end(); }

    VertexId conjugate(VertexId v) const noexcept { return vertex(v).conjugate(); }
    EdgeId conjugate(EdgeId e) const noexcept { return edge(e).conjugate(); }

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
    bool CheckUniqueOutgoingEdge(VertexId v) const noexcept {
        return OutgoingEdgeCount(v) == 1;
    }

    bool CheckUniqueIncomingEdge(VertexId v) const noexcept {
        return IncomingEdgeCount(v) == 1;
    }
    EdgeId GetUniqueIncomingEdge(VertexId v) const {
        VERIFY(CheckUniqueIncomingEdge(v));
        return *in_begin(v);
    }

    bool IsDeadEnd(VertexId v) const noexcept { return OutgoingEdgeCount(v) == 0; }
    bool IsDeadStart(VertexId v) const noexcept { return IncomingEdgeCount(v) == 0; }

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
    size_t operator()(omnigraph::impl::VertexId v) const noexcept {
        return v.hash();
    }
};

template<>
struct hash<omnigraph::impl::EdgeId> {
    size_t operator()(omnigraph::impl::EdgeId e) const noexcept {
        return e.hash();
    }
};

}
