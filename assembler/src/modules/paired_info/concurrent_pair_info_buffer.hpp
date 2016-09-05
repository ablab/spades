//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "histogram.hpp"
#include "histptr.hpp"

#include <btree/btree_map.h>
#include <cuckoo/cuckoohash_map.hh>

namespace omnigraph {

namespace de {

template<typename G, typename Traits, template<typename, typename> class Container>
class ConcurrentPairedBuffer : public PairedBufferBase<ConcurrentPairedBuffer<G, Traits, Container>,
                                                       G, Traits> {
    typedef ConcurrentPairedBuffer<G, Traits, Container> self;
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
    typedef cuckoohash_map<EdgeId, InnerMap> StorageMap;

  public:
    ConcurrentPairedBuffer(const Graph &g)
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
        if (!storage_.contains(e1))
            storage_.insert(e1, InnerMap()); // We can fail to insert here, it's ok

        size_t added = 0;
        typename InnerHistPtr::pointer inserted = nullptr;
        storage_.update_fn(e1,
                           [&](InnerMap &second) { // Now we will hold lock to the whole "subtree" starting from e1
                               if (!second.count(e2)) {
                                   inserted = new InnerHistogram();
                                   second.insert(std::make_pair(e2, InnerHistPtr(inserted, /* owning */ true)));
                               }
                               added = second[e2]->merge_point(p);
                           });

        return { inserted, added };
    }

    template<class OtherHist>
    std::pair<typename InnerHistPtr::pointer, size_t> InsertHist(EdgeId e1, EdgeId e2, const OtherHist &h) {
        if (!storage_.contains(e1))
            storage_.insert(e1, InnerMap()); // We can fail to insert here, it's ok

        size_t added = 0;
        typename InnerHistPtr::pointer inserted = nullptr;
        storage_.update_fn(e1,
                           [&](InnerMap &second) { // Now we will hold lock to the whole "subtree" starting from e1
                               if (!second.count(e2)) {
                                   inserted = new InnerHistogram();
                                   second.insert(std::make_pair(e2, InnerHistPtr(inserted, /* owning */ true)));
                               }
                               added = second[e2]->merge(h);
                           });

        return { inserted, added };
    }

    void InsertHistView(EdgeId e1, EdgeId e2, typename InnerHistPtr::pointer p) {
        if (!storage_.contains(e1))
            storage_.insert(e1, InnerMap()); // We can fail to insert here, it's ok

        storage_.update_fn(e1,
                           [&](InnerMap &second) { // Now we will hold lock to the whole "subtree" starting from e1
                               auto res = second.insert(std::make_pair(e2, InnerHistPtr(p, /* owning */ false)));
                               VERIFY_MSG(res.second, "Index insertion inconsistency");
                           });
    }

  protected:
    StorageMap storage_;
};

template<class Graph>
using ConcurrentPairedInfoBuffer = ConcurrentPairedBuffer<Graph, RawPointTraits, btree_map>;

} // namespace de

} // namespace omnigraph
