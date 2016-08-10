//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "histogram.hpp"
#include "histptr.hpp"

namespace omnigraph {

namespace de {

/**
 * @brief Index of paired reads information. For each pair of edges, we store so-called histogram which is a set
 *        of points with distance between those edges. Index is internally arranged as a map of map of histograms:
 *        edge1 -> (edge2 -> histogram)
 *        When we add a point (a,b)->p into the index, we automatically insert a conjugate point (b',a')->p',
 *        (self-conjugate edge pairs are the sole exception), so the index is always conjugate-symmetrical.
 *        Index provides access for a lot of different information:
 *        - if you need to have a histogram between two edges, use Get(edge1, edge2);
 *        - if you need to get a neighbourhood of some edge (second edges with corresponding histograms), use Get(edge1);
 *        - if you need to skip a symmetrical half of that neighbourhood, use GetHalf(edge1);
 *        Backward information (e.g., (b,a)->-p) is currently inaccessible.
 * @param G graph type
 * @param Traits Policy-like structure with associated types of inner and resulting points, and how to convert between them
 * @param C map-like container type (parameterized by key and value type)
 */

template<class Derived, class G, class Traits>
class PairedBufferBase {
  protected:
    typedef typename Traits::Gapped InnerPoint;

  public:
    typedef G Graph;
    typedef typename Graph::EdgeId EdgeId;
    typedef std::pair<EdgeId, EdgeId> EdgePair;
    typedef typename Traits::Expanded Point;

  public:
    PairedBufferBase(const Graph &g)
            : size_(0), graph_(g) {}

    //---------------- Data inserting methods ----------------
    /**
     * @brief Adds a point between two edges to the index,
     *        merging weights if there's already one with the same distance.
     */
    void Add(EdgeId e1, EdgeId e2, Point p) {
        InnerPoint sp = Traits::Shrink(p, CalcOffset(e1));
        InsertWithConj(e1, e2, sp);
    }

    /**
     * @brief Adds a whole set of points between two edges to the index.
     */
    template<typename TH>
    void AddMany(EdgeId e1, EdgeId e2, const TH& hist) {
        for (auto p : hist) {
            InnerPoint sp = Traits::Shrink(p, CalcOffset(e1));
            InsertWithConj(e1, e2, sp);
        }
    }
    //---------------- Miscellaneous ----------------

    /**
     * Returns the graph the index is based on. Needed for custom iterators.
     */
    const Graph &graph() const { return graph_; }

    /**
     * @brief Returns the physical index size (total count of all histograms).
     */
    size_t size() const { return size_; }

  public:
    /**
     * @brief Returns a conjugate pair for two edges.
     */
    EdgePair ConjugatePair(EdgeId e1, EdgeId e2) const {
        return std::make_pair(this->graph_.conjugate(e2), this->graph_.conjugate(e1));
    }
    /**
     * @brief Returns a conjugate pair for a pair of edges.
     */
    EdgePair ConjugatePair(EdgePair ep) const {
        return ConjugatePair(ep.first, ep.second);
    }

  private:
    void InsertWithConj(EdgeId e1, EdgeId e2, InnerPoint p) {
        EdgePair minep, maxep;
        std::tie(minep, maxep) = this->MinMaxConjugatePair({ e1, e2 });
        bool selfconj = this->IsSelfConj(e1, e2);

        auto res = static_cast<Derived*>(this)->InsertOne(minep.first, minep.second, p);
        size_t added = (selfconj ? res.second : 2 * res.second);
#       pragma omp atomic
        size_ += added;
        if (res.first && !selfconj)
            static_cast<Derived*>(this)->InsertHistView(maxep.first, maxep.second, res.first);
        else if (selfconj) // This would double the weigth of self-conjugate pairs
            static_cast<Derived*>(this)->InsertOne(minep.first, minep.second, p);
    }

  protected:
    template<class OtherHist>
    void Merge(EdgeId e1, EdgeId e2, const OtherHist &h) {
        EdgePair minep, maxep;
        std::tie(minep, maxep) = this->MinMaxConjugatePair({ e1, e2 });
        bool selfconj = this->IsSelfConj(e1, e2);

        auto res = static_cast<Derived*>(this)->InsertHist(minep.first, minep.second, h);
        size_t added = (selfconj ? res.second : 2 * res.second);
#       pragma omp atomic
        size_ += added;
        if (res.first && !selfconj)
            static_cast<Derived*>(this)->InsertHistView(maxep.first, maxep.second, res.first);
        else if (selfconj) // This would double the weigth of self-conjugate pairs
            static_cast<Derived*>(this)->InsertHist(minep.first, minep.second, h);
    }

    std::pair<EdgePair, EdgePair> MinMaxConjugatePair(EdgePair ep) const {
        EdgePair conj = ConjugatePair(ep);

        return (ep < conj ? std::make_pair(ep, conj) : std::make_pair(conj, ep));
    }

    bool IsSelfConj(EdgeId e1, EdgeId e2) const {
        return e1 == this->graph_.conjugate(e2);
    }

    size_t CalcOffset(EdgeId e) const {
        return this->graph().length(e);
    }

  protected:
    size_t size_;
    const Graph& graph_;
};


template<typename G, typename Traits, template<typename, typename> class Container>
class PairedBuffer : public PairedBufferBase<PairedBuffer<G, Traits, Container>,
                                             G, Traits> {
    typedef PairedBuffer<G, Traits, Container> self;
    typedef PairedBufferBase<self, G, Traits> base;

    friend class PairedBufferBase<self, G, Traits>;

  protected:
    using typename base::InnerPoint;
    typedef omnigraph::de::Histogram<InnerPoint> InnerHistogram;
    typedef omnigraph::de::StrongWeakPtr<InnerHistogram> InnerHistPtr;

  public:
    using typename base::Graph;
    using typename base::EdgeId;
    using typename base::EdgePair;
    using typename base::Point;

    typedef Container<EdgeId, InnerHistPtr> InnerMap;
    typedef Container<EdgeId, InnerMap> StorageMap;

  public:
    PairedBuffer(const Graph &g)
            : base(g) {
        clear();
    }

    //---------------- Miscellaneous ----------------

    /**
     * @brief Clears the whole index. Used in merging.
     */
    void clear() {
        storage_.clear();
        this->size_ = 0;
    }

    typename StorageMap::locked_table lock_table() {
        return storage_.lock_table();
    }

  private:
    std::pair<typename InnerHistPtr::pointer, size_t> InsertOne(EdgeId e1, EdgeId e2, InnerPoint p) {
        InnerMap& second = storage_[e1];
        typename InnerHistPtr::pointer inserted = nullptr;
        if (!second.count(e2)) {
            inserted = new InnerHistogram();
            second.insert(std::make_pair(e2, InnerHistPtr(inserted, /* owning */ true)));
        }

        size_t added = second[e2]->merge_point(p);

        return { inserted, added };
    }

    template<class OtherHist>
    std::pair<typename InnerHistPtr::pointer, size_t> InsertHist(EdgeId e1, EdgeId e2, const OtherHist &h) {
        InnerMap& second = storage_[e1];
        typename InnerHistPtr::pointer inserted = nullptr;
        if (!second.count(e2)) {
            inserted = new InnerHistogram();
            second.insert(std::make_pair(e2, InnerHistPtr(inserted, /* owning */ true)));
        }

        size_t added = second[e2]->merge(h);

        return { inserted, added };
    }

    void InsertHistView(EdgeId e1, EdgeId e2, typename InnerHistPtr::pointer p) {
        auto res = storage_[e1].insert(std::make_pair(e2, InnerHistPtr(p, /* owning */ false)));
        VERIFY_MSG(res.second, "Index insertion inconsistency");
    }

  protected:
    StorageMap storage_;
};

} // namespace de

} // namespace omnigraph
