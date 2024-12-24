//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2016-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "histogram.hpp"
#include "histptr.hpp"
#include "paired_info.hpp"
#include "paired_info_buffer.hpp"

#include <parallel_hashmap/btree.h>
#include <cuckoo/cuckoohash_map.hh>

#define XXH_INLINE_ALL
#include "xxh/xxhash.h"

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

    struct EdgeIdHasher {
        size_t operator()(EdgeId e) const {
            uint64_t h1 = e.hash();
            return XXH3_64bits(&h1, sizeof(h1));
        }
    };

    typedef Container<EdgeId, InnerHistPtr> InnerMap;
    typedef libcuckoo::cuckoohash_map<EdgeId, InnerMap, EdgeIdHasher> StorageMap;

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

    void BinWrite(std::ostream &str) {
        using io::binary::BinWrite;
        BinWrite<size_t>(str, storage_.size());
        const auto table = lock_table();
        for (const auto &i : table) {
            BinWrite(str, i.first.int_id());
            for (const auto &j : i.second) {
                if (j.second.owning()) {
                    BinWrite(str, j.first.int_id());
                    io::binary::BinWrite(str, *(j.second));
                }
            }
            BinWrite(str, (size_t)0); //null-term
        }
    }

    void BinRead(std::istream &str) {
        clear();
        auto table = lock_table();
        using io::binary::BinRead;
        auto storage_size = BinRead<size_t>(str);
        while (storage_size--) {
            auto e1 = BinRead<uint64_t>(str);
            while (true) {
                auto e2 = BinRead<uint64_t>(str);
                if (!e2) //null-term
                    break;
                auto hist = new InnerHistogram();
                io::binary::BinRead(str, *hist);
                TRACE(e1 << "->" << e2 << ": " << hist->size() << "points");
                table[e1][e2] = InnerHistPtr(hist, /* owning */ true);
                bool selfconj = this->IsSelfConj(e1, e2);
                size_t added = hist->size() * (selfconj ? 1 : 2);
                this->size_ += added;
                if (!selfconj) {
                    auto conj = this->ConjugatePair(e1, e2);
                    table[conj.first][conj.second] = InnerHistPtr(hist, /* owning */ false);
                }
            }
        }
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

template<class Graph>
using ConcurrentClusteredPairedInfoBuffer = ConcurrentPairedBuffer<Graph, PointTraits, btree_map>;

} // namespace de

} // namespace omnigraph
