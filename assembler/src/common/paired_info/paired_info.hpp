//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/adt/iterator_range.hpp"
#include <boost/iterator/iterator_facade.hpp>
#include <btree/safe_btree_map.h>

#include "paired_info_buffer.hpp"

#include <type_traits>

namespace omnigraph {

namespace de {

template<typename G, typename Traits, template<typename, typename> class Container>
class PairedIndex : public PairedBuffer<G, Traits, Container> {
    typedef PairedIndex<G, Traits, Container> self;
    typedef PairedBuffer<G, Traits, Container> base;

    typedef typename base::InnerHistogram InnerHistogram;
    typedef typename base::InnerHistPtr InnerHistPtr;
    typedef typename base::InnerPoint InnerPoint;

    using typename base::EdgePair;

public:
    using typename base::Graph;
    using typename base::EdgeId;
    typedef typename base::InnerMap InnerMap;
    typedef typename base::StorageMap StorageMap;
    using typename base::Point;

    typedef omnigraph::de::Histogram<Point> Histogram;

    //--Data access types--

    typedef typename StorageMap::const_iterator ImplIterator;

    //---------------- Data accessing methods ----------------

    /**
     * @brief Underlying raw implementation data (for custom iterator helpers).
     */
    ImplIterator data_begin() const {
        return this->storage_.begin();
    }

    /**
     * @brief Underlying raw implementation data (for custom iterator helpers).
     */
    ImplIterator data_end() const {
        return this->storage_.end();
    }

    /**
     * @brief Smart proxy set representing a composite histogram of points between two edges.
     * @detail You can work with the proxy just like any constant set.
     *         The only major difference is that it returns all consisting points by value,
     *         because some of them don't exist in the underlying sets and are
     *         restored from the conjugate info on-the-fly.
     */
    class HistProxy {

    public:
        /**
         * @brief Iterator over a proxy set of points.
         */
        class Iterator: public boost::iterator_facade<Iterator, Point, boost::bidirectional_traversal_tag, Point> {

            typedef typename InnerHistogram::const_iterator InnerIterator;

        public:
            Iterator(InnerIterator iter, DEDistance offset, bool back = false)
                    : iter_(iter), offset_(offset), back_(back)
            {}

        private:
            friend class boost::iterator_core_access;

            Point dereference() const {
                auto i = iter_;
                if (back_) --i;
                Point result = Traits::Expand(*i, offset_);
                if (back_)
                    result.d = -result.d;
                return result;
            }

            void increment() {
                back_ ? --iter_ : ++iter_;
            }

            void decrement() {
                back_ ? ++iter_ : --iter_;
            }

            inline bool equal(const Iterator &other) const {
                return iter_ == other.iter_ && back_ == other.back_;
            }

            InnerIterator iter_; //current position
            DEDistance offset_; //edge length
            bool back_;
        };

        /**
         * @brief Returns a wrapper for a histogram.
         */
        HistProxy(const InnerHistogram& hist, DEDistance offset = 0, bool back = false)
            : hist_(hist), offset_(offset), back_(back)
        {}

        /**
         * @brief Returns an empty proxy (effectively a Null object pattern).
         */
        static const InnerHistogram& empty_hist() {
            static InnerHistogram res;
            return res;
        }

        Iterator begin() const {
            return Iterator(back_ ? hist_.end() : hist_.begin(), offset_, back_);
        }

        Iterator end() const {
            return Iterator(back_ ? hist_.begin() : hist_.end(), offset_, back_);
        }

        /**
         * @brief Finds the point with the minimal distance.
         */
        Point min() const {
            VERIFY(!empty());
            return *begin();
        }

        /**
         * @brief Finds the point with the maximal distance.
         */
        Point max() const {
            VERIFY(!empty());
            return *--end();
        }

        /**
         * @brief Returns the copy of all points in a simple flat histogram.
         */
        Histogram Unwrap() const {
            return Histogram(begin(), end());
        }

        size_t size() const {
            return hist_.size();
        }

        bool empty() const {
            return hist_.empty();
        }

    private:
        const InnerHistogram& hist_;
        DEDistance offset_;
        bool back_;
    };

    typedef typename HistProxy::Iterator HistIterator;

    //---- Traversing edge neighbours ----

    using EdgeHist = std::pair<EdgeId, HistProxy>;

    /**
     * @brief A proxy map representing neighbourhood of an edge,
     *        where `Key` is the graph edge ID and `Value` is the proxy histogram.
     * @detail You can work with the proxy just like with any constant map.
     *         The only major difference is that it returns all consisting pairs by value,
     *         because proxies are constructed on-the-fly.
     */
    class EdgeProxy {
    public:

        /**
         * @brief Iterator over a proxy map.
         * @detail For a full proxy, traverses both straight and conjugate pairs.
         *         For a half proxy, traverses only lesser pairs (i.e., (a,b) where (a,b)<=(b',a')) of edges.
         */
        class Iterator: public boost::iterator_facade<Iterator, EdgeHist, boost::forward_traversal_tag, EdgeHist> {

            typedef typename InnerMap::const_iterator InnerIterator;

            void Skip() { //For a half iterator, skip conjugate pairs
                while (half_ && iter_ != stop_ && index_.GreaterPair(edge_, iter_->first))
                    ++iter_;
            }

        public:
            Iterator(const PairedIndex &index, InnerIterator iter, InnerIterator stop, EdgeId edge, bool half)
                    : index_ (index)
                    , iter_(iter)
                    , stop_(stop)
                    , edge_(edge)
                    , half_(half)
            {
                Skip();
            }

            void increment() {
                ++iter_;
                Skip();
            }

            void operator=(const Iterator &other) {
                //TODO: is this risky without an assertion?
                //VERIFY(index_ == other.index_);
                //We shouldn't reassign iterators from one index onto another
                iter_ = other.iter_;
                stop_ = other.stop_;
                edge_ = other.edge_;
                half_ = other.half_;
            }

        private:
            friend class boost::iterator_core_access;

            bool equal(const Iterator &other) const {
                return iter_ == other.iter_;
            }

            EdgeHist dereference() const {
                const auto& hist = *iter_->second;
                return std::make_pair(iter_->first, HistProxy(hist, index_.CalcOffset(edge_)));
            }

        private:
            const PairedIndex &index_; //TODO: get rid of this somehow
            InnerIterator iter_, stop_;
            EdgeId edge_;
            bool half_;
        };

        EdgeProxy(const PairedIndex &index, const InnerMap& map, EdgeId edge, bool half = false)
            : index_(index), map_(map), edge_(edge), half_(half)
        {}

        Iterator begin() const {
            return Iterator(index_, map_.begin(), map_.end(), edge_, half_);
        }

        Iterator end() const {
            return Iterator(index_, map_.end(), map_.end(), edge_, half_);
        }

        HistProxy operator[](EdgeId e2) const {
            if (half_ && index_.GreaterPair(edge_, e2))
                return HistProxy::empty_hist();
            return index_.Get(edge_, e2);
        }

        bool empty() const {
            return map_.empty();
        }

    private:
        const PairedIndex& index_;
        const InnerMap& map_;
        EdgeId edge_;
        //When false, represents all neighbours (consisting both of directly added data and "restored" conjugates).
        //When true, proxifies only half of the added edges.
        bool half_;
    };

    typedef typename EdgeProxy::Iterator EdgeIterator;

    //---------------- Constructor ----------------

    PairedIndex(const Graph &graph)
        : base(graph)
    {}

private:
    bool GreaterPair(EdgeId e1, EdgeId e2) const {
        auto ep = std::make_pair(e1, e2);
        return ep > this->ConjugatePair(ep);
    }

public:
    /**
     * @brief Adds a lot of info from another index, using fast merging strategy.
     *        Should be used instead of point-by-point index merge.
     */
    template<class Buffer>
    void Merge(Buffer& index_to_add) {
        if (index_to_add.size() == 0)
            return;

        auto locked_table = index_to_add.lock_table();
        for (auto& kvpair : locked_table) {
            EdgeId e1_to_add = kvpair.first; auto& map_to_add = kvpair.second;

            for (auto& to_add : map_to_add) {
                EdgePair ep(e1_to_add, to_add.first), conj = this->ConjugatePair(e1_to_add, to_add.first);
                if (ep > conj)
                    continue;

                base::Merge(ep.first, ep.second, *to_add.second);
            }
        }
        VERIFY(this->size() >= index_to_add.size());
    }

    template<class Buffer>
    typename std::enable_if<std::is_convertible<typename Buffer::InnerMap, InnerMap>::value,
        void>::type MoveAssign(Buffer& from) {
        auto& base_index = this->storage_;
        base_index.clear();
        auto locked_table = from.lock_table();
        for (auto& kvpair : locked_table) {
            base_index[kvpair.first] = std::move(kvpair.second);
        }
        this->size_ = from.size();
    }

public:
    //---------------- Data deleting methods ----------------

    /**
     * @brief Removes the specific entry from the index, and its conjugate.
     * @warning Don't use it on unclustered index, because hashmaps require set_deleted_item
     * @return The number of deleted entries (0 if there wasn't such entry)
     */
    size_t Remove(EdgeId e1, EdgeId e2, Point p) {
        InnerPoint point = Traits::Shrink(p, this->graph_.length(e1));

        // We remove first "non-owning part"
        EdgePair minep, maxep;
        std::tie(minep, maxep) = this->MinMaxConjugatePair({ e1, e2 });

        size_t res = RemoveSingle(minep.first, minep.second, point);
        size_t removed = (this->IsSelfConj(e1, e2) ? res : 2 * res);
        this->size_ -= removed;

        Prune(maxep.first, maxep.second);
        Prune(minep.first, minep.second);

        return removed;
    }

    /**
     * @brief Removes the whole histogram from the index, and its conjugate.
     * @warning Don't use it on unclustered index, because hashmaps require set_deleted_item
     * @return The number of deleted entries
     */
    size_t Remove(EdgeId e1, EdgeId e2) {
        EdgePair minep, maxep;
        std::tie(minep, maxep) = this->MinMaxConjugatePair({ e1, e2 });

        size_t removed = RemoveAll(maxep.first, maxep.second);
        removed += RemoveAll(minep.first, minep.second);
        this->size_ -= removed;

        return removed;
    }

  private:
    void Prune(EdgeId e1, EdgeId e2) {
        auto i1 = this->storage_.find(e1);
        if (i1 == this->storage_.end())
            return;

        auto& map = i1->second;
        auto i2 = map.find(e2);
        if (i2 == map.end())
            return;

        if (!i2->second->empty())
            return;

        map.erase(e2);
        if (map.empty())
            this->storage_.erase(e1);
    }

    size_t RemoveSingle(EdgeId e1, EdgeId e2, InnerPoint point) {
        auto i1 = this->storage_.find(e1);
        if (i1 == this->storage_.end())
            return 0;

        auto& map = i1->second;
        auto i2 = map.find(e2);
        if (i2 == map.end())
            return 0;

        if (!i2->second->erase(point))
           return 0;

        return 1;
    }

    size_t RemoveAll(EdgeId e1, EdgeId e2) {
        auto i1 = this->storage_.find(e1);
        if (i1 == this->storage_.end())
            return 0;
        auto& map = i1->second;
        auto i2 = map.find(e2);
        if (i2 == map.end())
            return 0;

        size_t size_decrease = i2->second->size();
        map.erase(i2);
        if (map.empty()) //Prune empty maps
            this->storage_.erase(i1);
        return size_decrease;
    }

public:

    /**
     * @brief Removes all neighbourhood of an edge (all edges referring to it, and their histograms)
     * @warning To keep the symmetricity, it also deletes all conjugates, so the actual complexity is O(size).
     * @return The number of deleted entries
     */
    size_t Remove(EdgeId edge) {
        InnerMap &inner_map = this->storage_[edge];
        std::vector<EdgeId> to_remove;
        to_remove.reserve(inner_map.size());
        size_t old_size = this->size();
        for (const auto& ep : inner_map)
            to_remove.push_back(ep.first);
        for (auto e2 : to_remove)
            this->Remove(edge, e2);
        return old_size - this->size();
    }

private:
    //When there is no such edge, returns a fake empty map for safety
    const InnerMap& GetImpl(EdgeId e) const {
        auto i = this->storage_.find(e);
        if (i != this->storage_.end())
            return i->second;
        return empty_map_;
    }

    //When there is no such histogram, returns a fake empty histogram for safety
    const InnerHistogram& GetImpl(EdgeId e1, EdgeId e2) const {
        auto i = this->storage_.find(e1);
        if (i != this->storage_.end()) {
            auto j = i->second.find(e2);
            if (j != i->second.end())
                return *j->second;
        }
        return HistProxy::empty_hist();
    }

public:

    /**
     * @brief Returns a whole proxy map to the neighbourhood of some edge.
     * @param e ID of starting edge
     */
    EdgeProxy Get(EdgeId e) const {
        return EdgeProxy(*this, GetImpl(e), e);
    }

    /**
     * @brief Returns a half proxy map to the neighbourhood of some edge.
     * @param e ID of starting edge
     */
    EdgeProxy GetHalf(EdgeId e) const {
        return EdgeProxy(*this, GetImpl(e), e, true);
    }

    /**
     * @brief Operator alias of Get(id).
     */
    EdgeProxy operator[](EdgeId e) const {
        return Get(e);
    }

    /**
     * @brief Returns a histogram proxy for all points between two edges.
     */
    HistProxy Get(EdgeId e1, EdgeId e2) const {
        return HistProxy(GetImpl(e1, e2), this->CalcOffset(e1));
    }

    /**
     * @brief Operator alias of Get(e1, e2).
     */
    HistProxy operator[](EdgePair p) const {
        return Get(p.first, p.second);
    }

    /**
     * @brief Checks if an edge (or its conjugated twin) is consisted in the index.
     */
    bool contains(EdgeId edge) const {
        return this->storage_.count(edge) + this->storage_.count(this->graph_.conjugate(edge)) > 0;
    }

    /**
     * @brief Checks if there is a histogram for two points (or their conjugated pair).
     */
    bool contains(EdgeId e1, EdgeId e2) const {
        auto i1 = this->storage_.find(e1);
        if (i1 != this->storage_.end() && i1->second.count(e2))
            return true;
        return false;
    }

    /**
     * @brief Inits the index with graph data. For each edge, adds a loop with zero weight.
     * @warning Do not call this on non-empty indexes.
     */
    void Init() {
        //VERIFY(size() == 0);
        for (auto it = this->graph_.ConstEdgeBegin(); !it.IsEnd(); ++it)
            this->Add(*it, *it, Point());
    }

private:
    InnerMap empty_map_; //null object
};

template<class T>
class NoLockingAdapter : public T {
  public:
    class locked_table {
      public:
        using iterator = typename T::iterator;
        using const_iterator = typename T::const_iterator;

        locked_table(T& table)
                : table_(table) {}

        iterator begin() { return table_.begin();  }
        const_iterator begin() const { return table_.begin(); }
        const_iterator cbegin() const { return table_.begin(); }

        iterator end() { return table_.end(); }
        const_iterator end() const { return table_.end(); }
        const_iterator cend() const { return table_.end(); }

        size_t size() const { return table_.size(); }

      private:
        T& table_;
    };

    // Nothing to lock here
    locked_table lock_table() {
        return locked_table(*this);
    }
};

//Aliases for common graphs
template<typename K, typename V>
using safe_btree_map = NoLockingAdapter<btree::safe_btree_map<K, V>>; //Two-parameters wrapper
template<typename Graph>
using PairedInfoIndexT = PairedIndex<Graph, PointTraits, safe_btree_map>;

template<typename K, typename V>
using btree_map = NoLockingAdapter<btree::btree_map<K, V>>; //Two-parameters wrapper

template<typename Graph>
using UnclusteredPairedInfoIndexT = PairedIndex<Graph, RawPointTraits, btree_map>;

/**
 * @brief A collection of paired indexes which can be manipulated as one.
 *        Used as a convenient wrapper in parallel index processing.
 */
template<class Index>
class PairedIndices {
    typedef std::vector<Index> Storage;
    Storage data_;

public:
    PairedIndices() {}

    PairedIndices(const typename Index::Graph& graph, size_t lib_num) {
        for (size_t i = 0; i < lib_num; ++i)
            data_.emplace_back(graph);
    }

    /**
     * @brief Initializes all indexes with zero points.
     */
    void Init() { for (auto& it : data_) it.Init(); }

    /**
     * @brief Clears all indexes.
     */
    void Clear() { for (auto& it : data_) it.clear(); }

    Index& operator[](size_t i) { return data_[i]; }

    const Index& operator[](size_t i) const { return data_[i]; }

    size_t size() const { return data_.size(); }

    typename Storage::iterator begin() { return data_.begin(); }
    typename Storage::iterator end() { return data_.end(); }

    typename Storage::const_iterator begin() const { return data_.begin(); }
    typename Storage::const_iterator end() const { return data_.end(); }
};

template<class Graph>
using PairedInfoIndicesT = PairedIndices<PairedInfoIndexT<Graph>>;

template<class Graph>
using UnclusteredPairedInfoIndicesT = PairedIndices<UnclusteredPairedInfoIndexT<Graph>>;

template<typename K, typename V>
using unordered_map = NoLockingAdapter<std::unordered_map<K, V>>; //Two-parameters wrapper
template<class Graph>
using PairedInfoBuffer = PairedBuffer<Graph, RawPointTraits, unordered_map>;

template<class Graph>
using PairedInfoBuffersT = PairedIndices<PairedInfoBuffer<Graph>>;

}

}
