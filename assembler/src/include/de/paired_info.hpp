//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <adt/iterator_range.hpp>
#include <btree/safe_btree_map.h>
#include <sparsehash/sparse_hash_map>

#include <type_traits>

#include "index_point.hpp"

namespace omnigraph {

namespace de {

/**
 * @brief Paired reads info storage. Arranged as a map of map of info points.
 * @param G graph type
 * @param H map-like container type (parameterized by key and value type)
 */
template<typename G, typename H, template<typename, typename> class Container>
class PairedIndex {

public:
    typedef G Graph;
    typedef H Histogram;
    typedef typename Graph::EdgeId EdgeId;
    typedef std::pair<EdgeId, EdgeId> EdgePair;
    typedef typename Histogram::value_type Point;

    typedef Container<EdgeId, Histogram> InnerMap;
    typedef Container<EdgeId, InnerMap> StorageMap;

    //--Data access types--

    typedef typename StorageMap::const_iterator ImplIterator;

public:
    /**
     * @brief Smart proxy set representing a composite histogram of points between two edges.
     * @detail You can work with the proxy just like with any constant set.
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

            typedef typename Histogram::const_iterator InnerIterator;

        public:
            Iterator(InnerIterator iter, bool back)
                    : iter_(iter), back_(back)
            {}

        private:
            friend class boost::iterator_core_access;

            Point dereference() const {
                if (back_) {
                    Point result = *(iter_ - 1);
                    result.d = -result.d;
                    return result;
                } else
                    return *iter_;
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
            bool back_;
        };

        /**
         * @brief Returns a wrapper for a histogram.
         */
        HistProxy(const Histogram& hist, bool back = false)
            : hist_(hist), back_(back)
        {}

        /**
         * @brief Returns an empty proxy (effectively a Null object pattern).
         */
        static const Histogram& empty_hist() {
            static Histogram res;
            return res;
        }

        Iterator begin() const {
            return Iterator(back_ ? hist_.end() : hist_.begin(), back_);
        }

        Iterator end() const {
            return Iterator(back_ ? hist_.begin() : hist_.end(), back_);
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
         * @brief Returns the copy of all points in a simple histogram.
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
        const Histogram& hist_;
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
     *         becauses some of them don't exist in the underlying sets and are
     *         restored from the conjugate info on-the-fly.
     */
    class EdgeProxy {
    public:

        /**
         * @brief Iterator over a proxy map.
         * @param full When true, traverses both straight and conjugate pairs.
         *             When false, traverses only lesser pairs of edges.
         */
        class Iterator: public boost::iterator_facade<Iterator, EdgeHist, boost::forward_traversal_tag, EdgeHist> {

            typedef typename InnerMap::const_iterator InnerIterator;

            bool SkipPair(EdgeId e1, EdgeId e2) {
                auto ep = std::make_pair(e1, e2);
                return ep > index_.ConjugatePair(ep);
            }

            void Skip() { //For a half iterator, skip conjugate pairs
                while (half_ && iter_ != stop_ && SkipPair(edge_, iter_->first))
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
                const auto& hist = iter_->second;
                return std::make_pair(iter_->first, HistProxy(hist));
            }

        private:
            const PairedIndex &index_;
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
            return index_.Get(edge_, e2);
        }

        //Currently unused
        /*HistProxy<true> GetBack(EdgeId e2) const {
            return index_.GetBack(edge_, e2);
        }*/

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

    //--Constructor--

    PairedIndex(const Graph &graph)
        : size_(0), graph_(graph)
    {}

    //--Inserting--
public:
    /**
     * @brief Returns a conjugate pair for two edges.
     */
    EdgePair ConjugatePair(EdgeId e1, EdgeId e2) const {
        return std::make_pair(graph_.conjugate(e2), graph_.conjugate(e1));
    }
    /**
     * @brief Returns a conjugate pair for a pair of edges.
     */
    EdgePair ConjugatePair(EdgePair ep) const {
        return ConjugatePair(ep.first, ep.second);
    }

private:
    void SwapConj(EdgeId &e1, EdgeId &e2) const {
        auto tmp = e1;
        e1 = graph_.conjugate(e2);
        e2 = graph_.conjugate(tmp);
    }

    void SwapConj(EdgeId &e1, EdgeId &e2, Point &p) const {
        SwapConj(e1, e2);
        p.d += CalcOffset(e1, e2);
    }

    float CalcOffset(EdgeId e1, EdgeId e2) const {
        return float(graph_.length(e1)) - float(graph_.length(e2));
    }

public:
    /**
     * @brief Adds a point between two edges to the index,
     *        merging weights if there's already one with the same distance.
     */
    void Add(EdgeId e1, EdgeId e2, Point point) {
        InsertWithConj(e1, e2, point);
    }

    /**
     * @brief Adds a whole set of points between two edges to the index.
     */
    template<typename TH>
    void AddMany(EdgeId e1, EdgeId e2, const TH& hist) {
        for (auto point : hist) {
            InsertWithConj(e1, e2, point);
        }
    }

private:

    void InsertWithConj(EdgeId e1, EdgeId e2,
                       Point sp) {
        size_ += storage_[e1][e2].merge_point(sp);
        //TODO: deal with loops and self-conj
        SwapConj(e1, e2, sp);
        size_ += storage_[e1][e2].merge_point(sp);
    }

    static bool IsSymmetric(EdgeId e1, EdgeId e2, Point point) {
        return (e1 == e2) && math::eq(point.d, 0.f);
    }

    bool IsSelfConj(EdgeId e1, EdgeId e2) {
        return e1 == graph_.conjugate(e2);
    }

public:
    /**
     * @brief Adds a lot of info from another index, using fast merging strategy.
     *        Should be used instead of point-by-point index merge.
     */
    template<class Index>
    void Merge(const Index& index_to_add) {
        auto& base_index = storage_;
        for (auto AddI = index_to_add.data_begin(); AddI != index_to_add.data_end(); ++AddI) {
            EdgeId e1_to_add = AddI->first;
            const auto& map_to_add = AddI->second;
            InnerMap& map_already_exists = base_index[e1_to_add];
            MergeInnerMaps(map_to_add, map_already_exists);
        }
        VERIFY(size() >= index_to_add.size());
    }

private:
    template<class OtherMap>
    void MergeInnerMaps(const OtherMap& map_to_add,
                        InnerMap& map) {
        for (const auto& to_add : map_to_add) {
            Histogram &hist_exists = map[to_add.first];
            size_ += hist_exists.merge(to_add.second);
        }
    }

public:
    //--Data deleting methods--

    /**
     * @brief Removes the specific entry from the index.
     * @warning Don't use it on unclustered index, because hashmaps require set_deleted_item
     * @return The number of deleted entries (0 if there wasn't such entry)
     */
    size_t Remove(EdgeId e1, EdgeId e2, Point point) {
        auto res = RemoveSingle(e1, e2, point);
        //TODO: deal with loops and self-conj
        SwapConj(e1, e2, point);
        res += RemoveSingle(e1, e2, point);
        return res;
    }

    /**
     * @brief Removes the whole histogram from the index.
     * @warning Don't use it on unclustered index, because hashmaps require set_deleted_item
     * @return The number of deleted entries
     */
    size_t Remove(EdgeId e1, EdgeId e2) {
        auto res = RemoveAll(e1, e2);
        if (!IsSelfConj(e1, e2)) { //TODO: loops?
            SwapConj(e1, e2);
            res += RemoveAll(e1, e2);
        }
        return res;
    }

private:

    //TODO: remove duplicode
    size_t RemoveSingle(EdgeId e1, EdgeId e2, Point point) {
        auto i1 = storage_.find(e1);
        if (i1 == storage_.end())
            return 0;
        auto& map = i1->second;
        auto i2 = map.find(e2);
        if (i2 == map.end())
            return 0;
        Histogram& hist = i2->second;
        if (!hist.erase(point))
           return 0;
        --size_;
        if (hist.empty()) { //Prune empty maps
            map.erase(e2);
            if (map.empty())
                storage_.erase(e1);
        }
        return 1;
    }

    size_t RemoveAll(EdgeId e1, EdgeId e2) {
        auto i1 = storage_.find(e1);
        if (i1 == storage_.end())
            return 0;
        auto& map = i1->second;
        auto i2 = map.find(e2);
        if (i2 == map.end())
            return 0;
        Histogram& hist = i2->second;
        size_t size_decrease = hist.size();
        map.erase(i2);
        size_ -= size_decrease;
        if (map.empty()) //Prune empty maps
            storage_.erase(i1);
        return size_decrease;
    }

public:

    /**
     * @brief Removes all neighbourhood of an edge (all edges referring to it, and their histograms)
     * @warning Currently doesn't check the conjugate info (should it?), so it may actually
     *          skip some data.
     * @return The number of deleted entries
     */
    size_t Remove(EdgeId edge) {
        InnerMap &inner_map = storage_[edge];
        for (auto iter = inner_map.begin(); iter != inner_map.end(); ++iter) {
            EdgeId e2 = iter->first;
            if (edge != e2) {
                this->Remove(e2, edge);
            }
        }
        size_t size_of_removed = inner_map.size();
        storage_.erase(edge);
        size_ -= size_of_removed;
        return size_of_removed;
    }

    // --Accessing--

    /**
     * @brief Underlying raw implementation data (for custom iterator helpers).
     */
    ImplIterator data_begin() const {
        return storage_.begin();
    }

    /**
     * @brief Underlying raw implementation data (for custom iterator helpers).
     */
    ImplIterator data_end() const {
        return storage_.end();
    }

    adt::iterator_range<ImplIterator> data() const {
        return adt::make_range(data_begin(), data_end());
    }

private:
    //When there is no such edge, returns a fake empty map for safety
    const InnerMap& GetImpl(EdgeId e) const {
        auto i = storage_.find(e);
        if (i != storage_.end())
            return i->second;
        return empty_map_;
    }

    //When there is no such histogram, returns a fake empty histogram for safety
    const Histogram& GetImpl(EdgeId e1, EdgeId e2) const {
        auto i = storage_.find(e1);
        if (i != storage_.end()) {
            auto j = i->second.find(e2);
            if (j != i->second.end())
                return j->second;
        }
        return HistProxy::empty_hist();
    }

public:

    /**
     * @brief Returns a proxy map to the neighbourhood of some edge.
     * @param e ID of starting edge
     * @param half true if edge pairs (a,b) where a > b should be skipped,
     *        when you don't care for conjugate info, or don't want to process them twice.
     *        Default is false.
     */
    EdgeProxy Get(EdgeId e, bool half = false) const {
        return EdgeProxy(*this, GetImpl(e), e, half);
    }

    /**
     * @brief Operator alias of Get(id).
     */
    EdgeProxy operator[](EdgeId e) const {
        return Get(e);
    }

    /**
     * @brief Returns a full histogram proxy for all points between two edges.
     */
    HistProxy Get(EdgeId e1, EdgeId e2) const {
        return HistProxy(GetImpl(e1, e2));
    }

    /**
     * @brief Operator alias of Get(e1, e2).
     */
    HistProxy operator[](EdgePair p) const {
        return Get(p.first, p.second);
    }

    //Currently unused
    /**
     * @brief Returns a full backwards histogram proxy for all points between two edges.
     */
    /*HistProxy<true> GetBack(EdgeId e1, EdgeId e2) const {
        return HistProxy<true>(GetImpl(e2, e1));
    }*/
    
    /**
     * @brief Checks if an edge (or its conjugated twin) is consisted in the index.
     */
    bool contains(EdgeId edge) const {
        return storage_.count(edge) + storage_.count(graph_.conjugate(edge)) > 0;
    }

    /**
     * @brief Checks if there is a histogram for two points.
     */
    bool contains(EdgeId e1, EdgeId e2) const {
        auto i1 = storage_.find(e1);
        if (i1 != storage_.end() && i1->second.count(e2))
            return true;
        return false;
    }

    // --Miscellaneous--

    /**
     * Returns the graph the index is based on. Needed for custom iterators.
     */
    const Graph &graph() const { return graph_; }

    /**
     * @brief Inits the index with graph data. Used in clustered indexes.
     * @warning Do not call this on non-empty indexes.
     */
    void Init() {
        //VERIFY(size() == 0);
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it)
            Add(*it, *it, Point());
    }

    /**
     * @brief Clears the whole index. Used in merging.
     */
    void Clear() {
        storage_.clear();
        size_ = 0;
    }

    /**
     * @brief Returns the physical index size (total count of all edge pairs)
     */
    size_t size() const { return size_; }

private:
    size_t size_;
    const Graph& graph_;
    StorageMap storage_;
    InnerMap empty_map_; //null object
};

//Aliases for common graphs
template<typename K, typename V>
using safe_btree_map = btree::safe_btree_map<K, V>; //Two-parameters wrapper
template<typename Graph>
using PairedInfoIndexT = PairedIndex<Graph, HistogramWithWeight, safe_btree_map>;

template<typename K, typename V>
using sparse_hash_map = google::sparse_hash_map<K, V>; //Two-parameters wrapper
template<typename Graph>
using UnclusteredPairedInfoIndexT = PairedIndex<Graph, RawHistogram, sparse_hash_map>;

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
        data_.reserve(lib_num);
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
    void Clear() { for (auto& it : data_) it.Clear(); }

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
using unordered_map = std::unordered_map<K, V>; //Two-parameters wrapper
template<class Graph>
using PairedInfoBuffer = PairedIndex<Graph, RawHistogram, unordered_map>;

template<class Graph>
using PairedInfoBuffersT = PairedIndices<PairedInfoBuffer<Graph>>;

}

}
