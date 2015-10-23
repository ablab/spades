//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/iterator/iterator_facade.hpp>
#include <btree/btree_set.h>
#include <btree/safe_btree_map.h>
#include <sparsehash/sparse_hash_map>

#include "index_point.hpp"

namespace omnigraph {

namespace de {

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

    template<bool full = true>
    class HistProxy {

    public:
        class Iterator: public boost::iterator_facade<Iterator, Point, boost::bidirectional_traversal_tag, Point> {
            typedef typename Histogram::const_iterator const_iterator;

        public:
            Iterator(const_iterator iter, const_iterator stop_iter, const_iterator end_iter, int offset, bool conj)
                    : iter_(iter)
                    , stop_iter_(stop_iter)
                    , jump_iter_(end_iter)
                    , conj_(conj)
                    , offset_(offset)
            {
                //Skip();
            }

            Point dereference() const {
                Point result = *iter_;
                if (conj_)
                    result.d += offset_;
                return result;
            }

            void increment() {
                ++iter_;
                if (full && !conj_ && iter_ == stop_iter_) {
                    conj_ = true;
                    iter_ = jump_iter_;
                }
            }

            void decrement() {
                if (full && conj_ && iter_ == jump_iter_) {
                    conj_ = false;
                    iter_ = stop_iter_;
                }
                --iter_;
            }

            bool equal(const Iterator &other) const {
                return iter_ == other.iter_;
            }
        private:
            const_iterator iter_, stop_iter_, jump_iter_;
            bool conj_;
            int offset_;
        };

        HistProxy(const Histogram& hist, const Histogram& conj_hist, int offset = 0)
            : hist_(hist)
            , conj_hist_(conj_hist)
            , offset_(offset)
        {}

        static const Histogram& empty_hist() {
            static Histogram res;
            return res;
        }

        HistProxy(const Histogram& hist, int offset = 0)
                : hist_(hist)
                , conj_hist_(HistProxy::empty_hist())
                , offset_(offset)
        {}

        Iterator begin() const {
            if (full) { //Only one branch will actually be in the binary
                auto conj = hist_.empty();
                auto start = conj ? conj_hist_.begin() : hist_.begin();
                return Iterator(start, hist_.end(), conj_hist_.begin(), offset_, conj);
            } else {
                auto stop = hist_.end();
                return Iterator(hist_.begin(), stop, stop, 0, false);
            }
        }
        Iterator end() const {
            if (full) { //Only one branch will actually be in the binary
                return Iterator(conj_hist_.end(), hist_.end(), conj_hist_.begin(), offset_, true);
            } else {
                auto stop = hist_.end();
                return Iterator(hist_.end(), stop, stop, 0, false);
            }
        }

        Point front() const {
            return begin().dereference();
        }

        Point back() const {
            return (--end()).dereference();
        }

        Histogram Unwrap() {
            return Histogram(begin(), end());
        }

        size_t size() const {
            return hist_.size() + conj_hist_.size();
        }

        bool empty() const {
            return hist_.empty() && conj_hist_.empty();
        }

    private:
        const Histogram& hist_, & conj_hist_;
        int offset_;
    };

    typedef HistProxy<true> FullHistProxy;
    typedef HistProxy<false> RawHistProxy;

    typedef typename HistProxy<true>::Iterator HistIterator;
    typedef typename HistProxy<false>::Iterator RawHistIterator;

    //---- Traversing edge neighbours ----

    template<bool full = true>
    using EdgeHist = std::pair<EdgeId, HistProxy<full>>;

    template<bool full = true>
    class EdgeProxy {
    public:

        class Iterator: public boost::iterator_facade<Iterator, EdgeHist<full>, boost::forward_traversal_tag, EdgeHist<full>> {
            typedef typename InnerMap::const_iterator const_iterator;

        public:
            Iterator(const PairedIndex &index, const_iterator iter,
                     const_iterator stop_iter, const_iterator jump_iter,
                     EdgeId edge, bool conj)
                    : index_ (index)
                    , iter_(iter)
                    , stop_iter_(stop_iter)
                    , jump_iter_(jump_iter)
                    , edge_(edge)
                    , conj_(conj)
            {}

            //friend class boost::iterator_core_access;

            void increment() {
                ++iter_;
                if (full) { //For a full iterator, jump from the straight submap onto a conjugate one
                    if (!conj_ && iter_ == stop_iter_) {
                        conj_ = true;
                        iter_ = jump_iter_;
                    }
                } else { //For a raw iterator, skip conjugate pairs
                    while (iter_ != stop_iter_ && iter_->first < edge_)
                        ++iter_;
                }
            }

            void operator=(const Iterator &other) {
                //VERIFY(index_ == other.index_); //TODO: is this risky?
                iter_ = other.iter_;
                stop_iter_ = other.stop_iter_;
                jump_iter_ = other.jump_iter_;
                edge_ = other.edge_;
                conj_ = other.conj_;
            }

            bool equal(const Iterator &other) const {
                return iter_ == other.iter_ && conj_ == other.conj_;
            }

            EdgeHist<full> dereference() const {
                if (full) {
                    int offset = index_.CalcOffset(edge_, iter_->first);
                    EdgePair conj = index_.ConjugatePair(edge_, iter_->first);
                    const auto& hist = conj_ ? index_.GetImpl(conj) : iter_->second;
                    const auto& conj_hist = conj_ ? iter_->second : index_.GetImpl(conj);
                    auto edge = conj_ ? conj.first : iter_->first;
                    return std::make_pair(edge, HistProxy<full>(hist, conj_hist, offset));
                } else {
                    const auto& hist = iter_->second;
                    const auto& conj_hist = index_.GetImpl(index_.ConjugatePair(edge_, iter_->first));
                    return std::make_pair(iter_->first, HistProxy<full>(hist, conj_hist));
                }
            }

        private:
            const PairedIndex &index_;
            const_iterator iter_, stop_iter_, jump_iter_;
            EdgeId edge_;
            bool conj_;
        };

        EdgeProxy(const PairedIndex &index, const InnerMap& map, const InnerMap& conj_map, EdgeId edge)
            : index_(index), map_(map), conj_map_(conj_map), edge_(edge)
        {}

        Iterator begin() const {
            if (full) { //Only one branch will actually be in the binary
                auto conj = map_.empty();
                auto start = conj ? conj_map_.begin() : map_.begin();
                return Iterator(index_, start, map_.end(), conj_map_.begin(), edge_, conj);
            } else {
                auto stop = map_.end();
                return Iterator(index_, map_.begin(), stop, stop, edge_, false);
            }
        }

        Iterator end() const {
            if (full) { //Only one branch will actually be in the binary
                auto stop = conj_map_.end();
                return Iterator(index_, stop, stop, stop, edge_, true);
            } else {
                auto stop = map_.end();
                return Iterator(index_, stop, stop, stop, edge_, false);
            }
        }

        HistProxy<full> operator[](EdgeId e2) const {
            //TODO: optimize
            const auto& hist = index_.GetImpl(edge_, e2);
            const auto& conj_hist = full ? index_.GetImpl(index_.ConjugatePair(edge_, e2)) : HistProxy<full>::empty_hist();
            return HistProxy<full>(hist, conj_hist, index_.CalcOffset(e2, edge_));
        }

        bool empty() const {
            return !index_.contains(edge_);
        }

        /*static const InnerMap& empty_map() {
            static InnerMap res;
            return res;
        }*/

    private:
        const PairedIndex& index_;
        const InnerMap& map_, conj_map_;
        EdgeId edge_;
    };

    typedef typename EdgeProxy<true>::Iterator EdgeIterator;
    typedef typename EdgeProxy<false>::Iterator RawEdgeIterator;

    /*typedef std::pair<EdgeId, EdgeProxy<true>> EdgeEdge;

    class Iterator : public boost::iterator_facade<Iterator, EdgeEdge, boost::forward_traversal_tag, EdgeEdge> {
    public:
        typedef typename InnerMap::const_iterator const_iterator;

        void increment() {

        }

        EdgeEdge dereference() const {


        }

    private:
        const PairedIndex& pi;
        const_iterator iter_, end_iter_;
        bool conj_;

    };*/

    //--Constructor--

    PairedIndex(const Graph &graph)
        : size_(0), graph_(graph)
    {}

    //--Inserting--
public:
    inline EdgePair ConjugatePair(EdgeId e1, EdgeId e2) const {
        return std::make_pair(graph_.conjugate(e2), graph_.conjugate(e1));
    }

    inline EdgePair ConjugatePair(EdgePair ep) const {
        return ConjugatePair(ep.first, ep.second);
    }
private:
    bool SwapConj(EdgeId &e1, EdgeId &e2) const {
        EdgePair ep = {e1, e2}, ep_conj = ConjugatePair(ep);
        if (ep > ep_conj) {
            e1 = ep_conj.first;
            e2 = ep_conj.second;
            return true;
        }
        return false;
    }

    bool SwapConj(EdgeId &e1, EdgeId &e2, Point &p) const {
        if (SwapConj(e1, e2)) {
            p.d += CalcOffset(e1, e2);
            return true;
        }
        return false;
    }

    int CalcOffset(EdgeId e1, EdgeId e2) const {
        return (int)graph_.length(e2) - (int)graph_.length(e1);
    }

public:
    //Adds a single pair info, merging histograms if there's already one
    void Add(EdgeId e1, EdgeId e2, Point point) {
        SwapConj(e1, e2, point);
        InsertOrMerge(e1, e2, point);
    }

    //Adds a whole histogram, merging histograms if there's already one
    void AddMany(EdgeId e1, EdgeId e2, const Histogram& hist) {
        bool swapped = SwapConj(e1, e2);
        for (auto point : hist) {
            if (swapped)
                point.d += CalcOffset(e1, e2);
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
        }
        else {
            InsertPoint(straight, sp);
            if (!IsSymmetric(e1, e2, sp)) {
                auto &reversed = storage_[e2][e1];
                InsertPoint(reversed, rp);
            }
        }
    }

    //Unstable for hash_map due to the iterator invalidation
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
        //UpdateSinglePoint(hist, to_update, *to_update + point_to_add);
        auto to_add = *to_update + to_merge;
        auto after_removed = hist.erase(to_update);
        hist.insert(after_removed, to_add);
    }

public:
    //Adds a lot of info, using fast merging strategy
    template<class Index>
    void Add(const Index& index_to_add) {
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
        typedef typename Histogram::iterator hist_iterator;
        for (auto &i : map_to_add) {
            Histogram &hist_exists = map[i.first];
            const auto& hist_to_add = i.second;

            for (auto new_point : hist_to_add) {
                const std::pair<hist_iterator, bool>& result = hist_exists.insert(new_point);
                if (!result.second) { // in this case we need to merge two points
                    MergeData(hist_exists, result.first, new_point);
                } else
                    ++size_;
            }
        }
    }

public:
    //--Deleting--
    //TODO: that's currently unsafe for unclustered index,
    //because hashmaps require set_deleted_item

    //Removes the specific entry
    // Returns the number of deleted entries
    size_t Remove(EdgeId e1, EdgeId e2, Point point) {
        SwapConj(e1, e2);
        auto res = RemoveSingle(e1, e2, point);
        if (!IsSymmetric(e1, e2, point))
            res += RemoveSingle(e2, e1, -point);
        return res;
    }
    // Removes the whole histogram
    // Returns the number of deleted entries
    size_t Remove(EdgeId e1, EdgeId e2) {
        SwapConj(e1, e2);
        auto res = RemoveSingle(e1, e2);
        if (e1 != e2)
            res += RemoveSingle(e2, e1);
        return res;
    }

private:

    //TODO: remove duplicode
    size_t RemoveSingle(EdgeId e1, EdgeId e2, Point point) {
        auto i1 = storage_.find(e1);
        if (i1 != storage_.end()) {
            auto& map = i1->second;
            auto i2 = map.find(e2);
            if (i2 != map.end()) {
                Histogram& hist = i2->second;
                if (hist.erase(point))
                    --size_;
                if (hist.empty())
                    map.erase(e2);
                if (map.empty())
                    storage_.erase(e1);
                return 1;
            }
        }
        return 0;
    }

    size_t RemoveSingle(EdgeId e1, EdgeId e2) {
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
                return 1;
            }
        }
        return 0;
    }

public:

    //Removes all points, which refer to this edge,
    //and all backward information
    //Returns the number of deleted entries
    size_t Remove(EdgeId edge) {
        InnerMap& inner_map = storage_[edge];
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

    //Underlying raw implementation data
    inline ImplIterator data_begin() const {
        return storage_.begin();
    }

    inline ImplIterator data_end() const {
        return storage_.end();
    }

    // Returns a proxy map to neighboring edges
    inline EdgeProxy<> Get(EdgeId id) const {
        return EdgeProxy<>(*this, GetImpl(id), GetImpl(graph_.conjugate(id)), id);
    }

    // Returns a proxy map to neighboring edges
    inline EdgeProxy<false> RawGet(EdgeId id) const {
        return EdgeProxy<false>(*this, GetImpl(id), empty_map_, id);
    }

    // Operator alias
    inline EdgeProxy<> operator[](EdgeId id) const {
        return Get(id);
    }

private:
    //Returns a fake map for safety
    const InnerMap& GetImpl(EdgeId e1) const {
        auto i = storage_.find(e1);
        if (i == storage_.end())
            return empty_map_;
        return i->second;
    }

    //Returns a fake histogram for safety
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

    // Returns a proxy set of points
    HistProxy<> Get(EdgeId e1, EdgeId e2) const {
        int offset = CalcOffset(e2, e1); //we need to calculate the offset with initial edges
        return HistProxy<>(GetImpl(e1, e2), GetImpl(ConjugatePair(e1, e2)), offset);
    }

    // Operator alias
    inline HistProxy<> operator[](EdgePair p) const {
        return Get(p.first, p.second);
    }


    HistProxy<false> RawGet(EdgeId e1, EdgeId e2) const {
        SwapConj(e1, e2);
        return HistProxy<false>(GetImpl(e1, e2), HistProxy<false>::empty_hist(), 0);
    }

    bool contains(EdgeId edge) const {
        return storage_.count(edge) + storage_.count(graph_.conjugate(edge)) > 0;
    }

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

    //Returns the graph the index is based on
    const Graph &graph() const { return graph_; }

    //Inits the index with graph data
    void Init() {
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it)
            Add(*it, *it, Point());
    }

    //Clears the whole index
    void Clear() {
        storage_.clear();
        size_ = 0;
    }

    // Returns the total index size
    size_t size() const { return size_; }

private:
    size_t size_;
    const Graph& graph_;
    StorageMap storage_;
    InnerMap empty_map_;   //null object
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

template<class Index>
class PairedIndices {
    typedef std::vector<Index> Storage;
    Storage data_;

public:

    PairedIndices(const typename Index::Graph& graph, size_t lib_num) {
        data_.reserve(lib_num);
        for (size_t i = 0; i < lib_num; ++i)
            data_.emplace_back(graph);
    }

    void Init() { for (auto& it : data_) it.Init(); }

    void Clear() { for (auto& it : data_) it.Clear(); }

    inline Index& operator[](size_t i) { return data_[i]; }

    inline const Index& operator[](size_t i) const { return data_[i]; }

    inline size_t size() const { return data_.size(); }

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
