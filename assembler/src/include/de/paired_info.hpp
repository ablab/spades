//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "conj_iterator.hpp"
#include "index_point.hpp"

#include <adt/iterator_range.hpp>

#include <btree/safe_btree_map.h>
#include <sparsehash/sparse_hash_map>


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

    /**
     * @brief Smart proxy set representing a composite histogram of points between two edges.
     * @param full When true, represents the whole histogram (consisting both of directly added points
     *             and "restored" conjugates).
     *             When false, proxifies only the added points.
     * @detail You can work with the proxy just like with any constant set.
     *         The only major difference is that it returns all consisting points by value,
     *         becauses some of them don't exist in the underlying sets and are
     *         restored from the conjugate info on-the-fly.
     */
    template<bool full = true>
    class HistProxy {

    public:
        /**
         * @brief Iterator over a proxy set of points.
         * @warning Generally, the proxy is unordered even if the set is ordered.
         *          If you require that, convert it into a flat histogram with Unwrap().
         * @param full When true, traverses both straight and conjugate points,
         *             and automatically recalculates the distance for latter.
         *             When false, traverses only the added points and skips the rest.
         */
        class Iterator: public boost::iterator_facade<Iterator, Point, boost::bidirectional_traversal_tag, Point> {

            typedef typename ConjProxy<Histogram>::Iterator InnerIterator;

        public:
            Iterator(InnerIterator iter, float offset)
                    : iter_(iter)
                    , offset_(offset)
            {}

        private:
            friend class boost::iterator_core_access;

            Point dereference() const {
                Point result = *iter_;
                if (iter_.Conj())
                    result.d += offset_;
                return result;
            }

            void increment() {
                ++iter_;
            }

            void decrement() {
                --iter_;
            }

            inline bool equal(const Iterator &other) const {
                return iter_ == other.iter_;
            }

            InnerIterator iter_; //current position
            float offset_;       //offset to be added for conjugate distance
        };

        HistProxy(const Histogram& hist, const Histogram& conj_hist, float offset = 0)
            : hist_(hist, conj_hist)
            , offset_(offset)
        {}

        /**
         * @brief Returns an empty proxy (effectively a Null object pattern).
         */
        static const Histogram& empty_hist() {
            static Histogram res;
            return res;
        }

        /**
         * @brief Returns a wrapper for an ordinary histogram (for implicit conversions)
         */
        HistProxy(const Histogram& hist, float offset = 0)
                : hist_(hist, HistProxy::empty_hist())
                , offset_(offset)
        {}

        Iterator begin() const {
            return Iterator(hist_.begin(), offset_);
        }

        Iterator end() const {
            //auto i = full ? hist_.end() : hist_.conj_begin();
            //return Iterator(i, offset_);
            return Iterator(hist_.end(), offset_);
        }

        /**
         * @brief Finds the point with the minimal distance.
         * @todo Simplify
         */
        Point min() const {
            //Our histograms are ordered, so the minimum is `begin` of either
            //straight or conjugate half, but we should beware of emptiness.
            VERIFY(!empty());
            auto i1 = begin();
            if (full) {
                auto i2 = Iterator(hist_.conj_begin(), offset_);
                if (i1 == i2 || i2 == end())
                    return *i1;
                return std::min(*i1, *i2);
            } else {
                return *i1;
            }
        }

        /**
         * @brief Finds the point with the maximal distance.
         * @todo Simplify
         */
        Point max() const {
            //Our histograms are ordered, so the maximum is `rbegin` of either
            //straight or conjugate half, but we should beware of emptiness.
            VERIFY(!empty());
            auto i1 = end();
            if (full) {
                auto i2 = Iterator(hist_.conj_begin(), offset_);
                if (i1 == i2 || i2 == begin())
                    return *--i1;
                return std::max(*--i1, *--i2);
            } else {
                return *--i1;
            }
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
        const ConjProxy<Histogram> hist_;
        float offset_;
    };

    /**
     * @brief Type synonym for full histogram proxies (with added and conjugated points)
     */
    typedef HistProxy<true> FullHistProxy;
    /**
     * @brief Type synonym for raw histogram proxies (only with directly added points)
     */
    typedef HistProxy<false> RawHistProxy;

    typedef typename HistProxy<true>::Iterator HistIterator;
    typedef typename HistProxy<false>::Iterator RawHistIterator;

    //---- Traversing edge neighbours ----

    template<bool full = true>
    using EdgeHist = std::pair<EdgeId, HistProxy<full>>;

    /**
     * @brief A proxy map representing neighbourhood of an edge,
     *        where `Key` is the graph edge ID and `Value` is the proxy histogram.
     * @param full When true, represents all neighbours (consisting both of directly added data
     *             and "restored" conjugates).
     *             When false, proxifies only the added edges.
     * @detail You can work with the proxy just like with any constant map.
     *         The only major difference is that it returns all consisting pairs by value,
     *         becauses some of them don't exist in the underlying sets and are
     *         restored from the conjugate info on-the-fly.
     */
    template<bool full = true>
    class EdgeProxy {
    public:

        /**
         * @brief Iterator over a proxy map.
         * @param full When true, traverses both straight and conjugate pairs,
         *             and automatically recalculates the distance for latter.
         *             When false, traverses only the added points and skips the rest.
         */
        class Iterator: public boost::iterator_facade<Iterator, EdgeHist<full>, boost::forward_traversal_tag, EdgeHist<full>> {

            typedef typename ConjProxy<InnerMap>::Iterator InnerIterator;

            void Skip() { //For a raw iterator, skip conjugate pairs
                while (!full && !iter_.Conj() && iter_->first < edge_)
                    ++iter_;
            }

        public:
            Iterator(const PairedIndex &index, InnerIterator iter, EdgeId edge)
                    : index_ (index)
                    , iter_(iter)
                    , edge_(edge)
            {
                Skip();
            }

            void increment() {
                ++iter_;
                Skip();
            }

            void operator=(const Iterator &other) {
                //TODO: is this risky without an assertion?
                //We shouldn't reassign iterators from one index onto another
                iter_ = other.iter_;
                edge_ = other.edge_;
            }

        private:
            friend class boost::iterator_core_access;

            bool equal(const Iterator &other) const {
                return iter_ == other.iter_;
            }

            EdgeHist<full> dereference() const {
                if (full) {
                    float offset = index_.CalcOffset(edge_, iter_->first);
                    EdgePair conj = index_.ConjugatePair(edge_, iter_->first);
                    if (iter_.Conj()) {
                        return std::make_pair(conj.first,
                            HistProxy<full>(index_.GetImpl(edge_, conj.first),
                                            index_.GetImpl(iter_->first, conj.second),
                                            offset));
                    } else {
                        return std::make_pair(iter_->first,
                            HistProxy<full>(iter_->second, index_.GetImpl(conj), offset));
                    }
                } else {
                    const auto& hist = iter_->second;
                    const auto& conj_hist = index_.GetImpl(index_.ConjugatePair(edge_, iter_->first));
                    return std::make_pair(iter_->first, HistProxy<full>(hist, conj_hist));
                }
            }

        private:
            const PairedIndex &index_;
            InnerIterator iter_;
            EdgeId edge_;
        };

        EdgeProxy(const PairedIndex &index, const InnerMap& map, const InnerMap& conj_map, EdgeId edge)
            : index_(index), map_(map, conj_map), edge_(edge)
        {}

        Iterator begin() const {
            return Iterator(index_, map_.begin(), edge_);
        }

        Iterator end() const {
            auto i = full ? map_.end() : map_.conj_begin();
            return Iterator(index_, i, edge_);
        }

        HistProxy<full> operator[](EdgeId e2) const {
            //TODO: optimize
            EdgeId e1 = edge_;
            auto offset = index_.CalcOffset(e1, e2);
            if (full) {
                const auto& hist = index_.GetImpl(edge_, e2);
                const auto& conj_hist = index_.GetImpl(index_.ConjugatePair(edge_, e2));
                return HistProxy<full>(hist, conj_hist, offset);
            } else {
                if (index_.SwapConj(e1, e2))
                    return HistProxy<full>(HistProxy<full>::empty_hist(), index_.GetImpl(e1, e2), offset);
                else
                    return HistProxy<full>(index_.GetImpl(e1, e2));
            }
        }

        inline bool empty() const {
            return map_.empty();
        }

    private:
        const PairedIndex& index_;
        const ConjProxy<InnerMap> map_;
        EdgeId edge_;
    };

    /*template<> HistProxy<true> EdgeProxy<true>::operator[](EdgeId e2) const {
        return index_.Get(edge_, e2);
    }

    template<> HistProxy<false> EdgeProxy<false>::operator[](EdgeId e2) const {
        return index_.RawGet(edge_, e2);
    }*/

    typedef typename EdgeProxy<true>::Iterator EdgeIterator;
    typedef typename EdgeProxy<false>::Iterator RawEdgeIterator;

    //--Constructor--

    PairedIndex(const Graph &graph)
        : size_(0), graph_(graph)
    {}

    //--Inserting--
public:
    /**
     * @brief Returns a conjugate pair for two edges.
     */
    inline EdgePair ConjugatePair(EdgeId e1, EdgeId e2) const {
        return std::make_pair(graph_.conjugate(e2), graph_.conjugate(e1));
    }
    /**
     * @brief Returns a conjugate pair for a pair of edges.
     */
    inline EdgePair ConjugatePair(EdgePair ep) const {
        return ConjugatePair(ep.first, ep.second);
    }

    bool SwapConj(EdgeId &e1, EdgeId &e2) const {
        EdgePair ep = {e1, e2}, ep_conj = ConjugatePair(ep);
        if (ep > ep_conj) {
            e1 = ep_conj.first;
            e2 = ep_conj.second;
            return true;
        }
        return false;
    }

private:
    bool SwapConj(EdgeId &e1, EdgeId &e2, Point &p) const {
        if (SwapConj(e1, e2)) {
            p.d += CalcOffset(e1, e2);
            return true;
        }
        return false;
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
        SwapConj(e1, e2, point);
        InsertOrMerge(e1, e2, point);
    }

    /**
     * @brief Adds a whole set of points between two edges to the index.
     */
    template<typename TH>
    void AddMany(EdgeId e1, EdgeId e2, const TH& hist) {
        float offset = SwapConj(e1, e2) ? CalcOffset(e1, e2) : 0.0f;
        for (auto point : hist) {
            point.d += offset;
            InsertOrMerge(e1, e2, point);
        }
    }

private:

    void InsertOrMerge(EdgeId e1, EdgeId e2,
                       const Point &sp) {
        auto& straight = storage_[e1][e2];
        auto si = straight.find(sp);
        auto rp = -sp;
        if (si != straight.end()) {
            MergeData(straight, si, sp);
            if (!IsSymmetric(e1, e2, sp)) {
                auto& reversed = storage_[e2][e1];
                auto ri = reversed.find(rp);
                MergeData(reversed, ri, rp);
            }
        } else {
            InsertPoint(straight, sp);
            if (!IsSymmetric(e1, e2, sp)) {
                auto &reversed = storage_[e2][e1];
                InsertPoint(reversed, rp);
            }
        }
    }

    //Would be faster, but unstable for hash_map due to the iterator invalidation
    /*void InsertOrMerge(Histogram& straight, Histogram& reversed,
                       const Point &sp) {
        auto si = straight.find(sp);
        auto rp = -sp;
        if (si != straight.end()) {
            MergeData(straight, si, sp);
            auto ri = reversed.find(rp);
            MergeData(reversed, ri, rp);
        }
        else {
            InsertPoint(reversed, rp);
            InsertPoint(straight, sp);
            //if (!IsSymmetric(e1, e2, point)) TODO

        }
    }*/

    static bool IsSymmetric(EdgeId e1, EdgeId e2, Point point) {
        return (e1 == e2) && math::eq(point.d, 0.f);
    }

    // modifying the histogram
    inline void InsertPoint(Histogram& histogram, Point point) {
        histogram.insert(point);
        ++size_;
    }

    void MergeData(Histogram& hist, typename Histogram::iterator to_update, const Point& to_merge) {
        //We can't just modify the existing point, because if variation is non-zero,
        //resulting distance will differ
        auto to_add = *to_update + to_merge;
        auto after_removed = hist.erase(to_update);
        hist.insert(after_removed, to_add);
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
        auto res = RemoveImpl(e1, e2, point);
        auto conj = ConjugatePair(e1, e2);
        point.d += CalcOffset(e2, e1);
        res += RemoveImpl(conj.first, conj.second, point);
        return res;
    }

    /**
     * @brief Removes the whole histogram from the index.
     * @warning Don't use it on unclustered index, because hashmaps require set_deleted_item
     * @return The number of deleted entries
     */
    size_t Remove(EdgeId e1, EdgeId e2) {
        SwapConj(e1, e2);
        auto res = RemoveAll(e1, e2);
        if (e1 != e2)
            res += RemoveAll(e2, e1);
        return res;
    }

private:

    size_t RemoveImpl(EdgeId e1, EdgeId e2, Point point) {
        auto res = RemoveSingle(e1, e2, point);
        if (!IsSymmetric(e1, e2, point))
            res += RemoveSingle(e2, e1, -point);
        return res;
    }

    //TODO: remove duplicode
    size_t RemoveSingle(EdgeId e1, EdgeId e2, Point point) {
        auto i1 = storage_.find(e1);
        if (i1 != storage_.end()) {
            auto& map = i1->second;
            auto i2 = map.find(e2);
            if (i2 != map.end()) {
                Histogram& hist = i2->second;
                if (hist.erase(point)) {
                    --size_;
                    if (hist.empty()) {
                        map.erase(e2);
                        if (map.empty())
                            storage_.erase(e1);
                    }
                    return 1;
                }
                return 0;
            }
        }
        return 0;
    }

    size_t RemoveAll(EdgeId e1, EdgeId e2) {
        auto i1 = storage_.find(e1);
        if (i1 != storage_.end()) {
            auto& map = i1->second;
            auto i2 = map.find(e2);
            if (i2 != map.end()) {
                Histogram& hist = i2->second;
                size_t size_decrease = hist.size();
                map.erase(i2);
                size_ -= size_decrease;
                if (map.empty())
                    storage_.erase(i1);
                return size_decrease;
            }
        }
        return 0;
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

    /**
     * @brief Returns a full proxy map to the neighbourhood of some edge.
     */
    EdgeProxy<> Get(EdgeId id) const {
        return EdgeProxy<>(*this, GetImpl(id), GetImpl(graph_.conjugate(id)), id);
    }

    /**
     * @brief Returns a raw proxy map to neighboring edges
     * @detail You should use it when you don't care for backward
     *         and conjugate info, or don't want to process them twice.
     */
    EdgeProxy<false> RawGet(EdgeId id) const {
        return EdgeProxy<false>(*this, GetImpl(id), empty_map_, id);
    }

    /**
     * @brief Operator alias of Get(id).
     */
    EdgeProxy<> operator[](EdgeId id) const {
        return Get(id);
    }

private:
    //When there is no such edge, returns a fake empty map for safety
    const InnerMap& GetImpl(EdgeId e1) const {
        auto i = storage_.find(e1);
        if (i == storage_.end())
            return empty_map_;
        return i->second;
    }

    //When there is no such histogram, returns a fake empty histogram for safety
    const Histogram& GetImpl(EdgeId e1, EdgeId e2) const {
        auto i = storage_.find(e1);
        if (i != storage_.end()) {
            auto j = i->second.find(e2);
            if (j != i->second.end())
                return j->second;
        }
        return HistProxy<true>::empty_hist();
    }

    inline const Histogram& GetImpl(EdgePair e) const {
        return GetImpl(e.first, e.second);
    }

public:

    /**
     * @brief Returns a full histogram proxy for all points between two edges.
     */
    HistProxy<> Get(EdgeId e1, EdgeId e2) const {
        auto offset = CalcOffset(e1, e2);
        return HistProxy<>(GetImpl(e1, e2), GetImpl(ConjugatePair(e1, e2)), offset);
    }

    /**
     * @brief Operator alias of Get(e1, e2).
     */
    inline HistProxy<> operator[](EdgePair p) const {
        return Get(p.first, p.second);
    }

    /**
     * @brief Returns a raw histogram proxy for only straight points between two edges.
     */
    HistProxy<false> RawGet(EdgeId e1, EdgeId e2) const {
        if (SwapConj(e1, e2))
            return HistProxy<false>(HistProxy<false>::empty_hist(), GetImpl(e1, e2), CalcOffset(e1, e2));
        else
            return HistProxy<false>(GetImpl(e1, e2), HistProxy<false>::empty_hist(), 0);
    }

    /**
     * @brief Checks if an edge (or its conjugated twin) is consisted in the index.
     */
    bool contains(EdgeId edge) const {
        return storage_.count(edge) + storage_.count(graph_.conjugate(edge)) > 0;
    }

    /**
     * @brief Checks if there is a histogram for two points (or their conjugate pair).
     */
    bool contains(EdgeId e1, EdgeId e2) const {
        auto conj = ConjugatePair(e1, e2);
        auto i1 = storage_.find(e1);
        if (i1 != storage_.end() && i1->second.count(e2))
            return true;
        auto i2 = storage_.find(conj.first);
        if (i2 != storage_.end() && i2->second.count(conj.second))
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
     */
    void Init() {
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
     * @warning (not really total, doesn't include the conjugate info)
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
     * @brief Inits all indexes.
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

/*
//Debug
template<typename T>
std::ostream& operator<<(std::ostream& str, const PairedInfoBuffer<T>& pi) {
    str << "--- PI of size " << pi.size() << "---\n";

    for (auto i = pi.data_begin(); i != pi.data_end(); ++i) {
        auto e1 = i->first;
        str << e1 << " has: \n";

        for (auto j = i->second.begin(); j != i->second.end(); ++j) {
            str << "- " << j->first << ": ";
            for (auto p : j->second)
                str << p << ", ";
            str << std::endl;
        }
    }

    str << "-------\n";
    return str;
}

//Debug
template<typename T>
std::ostream& operator<<(std::ostream& str, const PairedInfoIndexT<T>& pi) {
    str << "--- PI of size " << pi.size() << "---\n";

    for (auto i = pi.data_begin(); i != pi.data_end(); ++i) {
        auto e1 = i->first;
        str << e1 << " has: \n";

        for (auto j = i->second.begin(); j != i->second.end(); ++j) {
            str << "- " << j->first << ": ";
            for (auto p : j->second)
                str << p << ", ";
            str << std::endl;
        }
    }

    str << "-------\n";
    return str;
}
*/

}

}
