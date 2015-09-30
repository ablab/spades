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

template<typename Graph, typename Hist, template<typename, typename> class Container>
class PairedIndex {

public:
    typedef Hist Histogram;
    typedef typename Graph::EdgeId EdgeId;
    typedef std::pair<EdgeId, EdgeId> EdgePair;
    typedef typename Hist::value_type Point;

    typedef Container<EdgeId, Hist> InnerMap;
    typedef Container<EdgeId, InnerMap> StorageMap;

    //--Data access types--

    typedef typename StorageMap::const_iterator ImplIterator;

    template<bool full = true>
    class HistProxy {

    public:
        class Iterator: public boost::iterator_facade<Iterator, Point, boost::forward_traversal_tag, Point> {
            typedef typename Histogram::const_iterator const_iterator;

        public:
            Iterator(const_iterator iter, const_iterator end_iter, int offset = 0)
                : iter_(iter)
                , end_iter_(end_iter)
                , offset_(offset)
            {}

            Point dereference() const {
                Point result = *iter_;
                result.d += offset_;
                return result;
            }

            void increment() {
                ++iter_;
                while (!full && (iter_ != end_iter_) && (math::le(dereference().d, 0.0f)))
                    ++iter_;
            }

            bool equal(const Iterator &other) const {
                return iter_ == other.iter_;
            }
        private:
            const_iterator iter_, end_iter_;
            int offset_;
        };

        HistProxy(const Histogram &hist, int offset = 0)
            : hist_(hist)
            , offset_(offset)
        {}

        Iterator begin() const { return Iterator(hist_.begin(), hist_.end(), offset_); }
        Iterator end() const { return Iterator(hist_.end(), hist_.end(), offset_); }

        Point front() const { return begin().dereference(); }
        Point back() const { return Iterator(--hist_.end(), hist_.end(), offset_).dereference(); }

        Histogram Unwrap() {
            return Histogram(begin(), end());
        }

        size_t size() const;

    private:
        const Histogram &hist_;
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
            Iterator(const PairedIndex &index, const_iterator iter, EdgeId edge, bool conj)
                    : index_ (index)
                    , iter_(iter)
                    , edge_(edge)
                    , conj_(conj)
            {}

            //friend class boost::iterator_core_access;

            void increment() {
                ++iter_;
                if (full && !conj_ && iter_ == index_.storage_.find(edge_)->second.end()) {
                    conj_ = true;
                    iter_ = index_.storage_.find(index_.graph_.conjugate(edge_))->second.begin();
                }
            }

            void operator=(const Iterator &other) {
                //VERIFY(index_ == other.index_); //TODO: is this risky?
                iter_ = other.iter_;
                edge_ = other.edge_;
                conj_ = other.conj_;
            }

            bool equal(const Iterator &other) const {
                return iter_ == other.iter_;
            }

            EdgeHist<full> dereference() const {
                if (conj_) {
                    int offset = (int)index_.graph_.length(edge_) - (int)index_.graph_.length(iter_->first);
                    return std::make_pair(index_.graph_.conjugate(iter_->first),
                                          HistProxy<full>(iter_->second, offset));
                }
                else
                    return std::make_pair(iter_->first, HistProxy<full>(iter_->second));
            }

        private:
            const PairedIndex &index_;
            const_iterator iter_;
            EdgeId edge_;
            bool conj_;
        };

        EdgeProxy(const PairedIndex &index, EdgeId edge)
            : index_(index), edge_(edge)
        {}

        Iterator begin() const {
            auto i = index_.storage_.find(edge_);
            return Iterator(index_, i->second.begin(), edge_, false);
        }
        Iterator end() const {
            auto e = full ? index_.graph_.conjugate(edge_) : edge_;
            auto i = index_.storage_.find(e);
            return Iterator(index_, i->second.end(), edge_, full);
        }

        HistProxy<full> operator[](EdgeId e2) const {
            auto i = index_.storage_.find(edge_);
            if (i != index_.storage_.end()) {
                auto j = i->second.find(e2);
                if (j != i->second.end()) {
                    return HistProxy<full>(j->second, 0);
                }
            }
            if (full) {
                i = index_.storage_.find(index_.graph_.conjugate(edge_));
                if (i != index_.storage_.end()) {
                    auto j = i->second.find(index_.graph_.conjugate(edge_));
                    if (j != i->second.end()) {
                        return HistProxy<full>(j->second, index_.CalcOffset(edge_, e2));
                    }
                }
            }
            VERIFY(false);
        }

        size_t size() const;
    private:
        const PairedIndex &index_;
        EdgeId edge_;
    };

    typedef typename EdgeProxy<true>::Iterator EdgeIterator;
    typedef typename EdgeProxy<false>::Iterator RawEdgeIterator;

    //--Constructor--

    PairedIndex(const Graph &graph)
        : size_(0), graph_(graph)
    {}

    //--Inserting--
public:
    EdgePair ConjugatePair(EdgePair ep) const {
        return std::make_pair(graph_.conjugate(ep.second), graph_.conjugate(ep.first));
    }
private:
    bool SwapConj(EdgeId &e1, EdgeId &e2) const {
        EdgePair ep(e1, e2), ep_conj = ConjugatePair(ep);
        if (ep > ep_conj) {
            e1 = ep.first;
            e2 = ep.second;
            return true;
        }
        return false;
    }

    bool SwapConj(EdgeId &e1, EdgeId &e2, Point &p) const {
        EdgePair ep(e1, e2), ep_conj = ConjugatePair(ep);
        if (ep > ep_conj) {
            p.d += CalcOffset(e1, e2);
            e1 = ep.first;
            e2 = ep.second;
            return true;
        }
        return false;
    }

    int CalcOffset(EdgeId e1, EdgeId e2) const {
        return (int)graph_.length(e1) - (int)graph_.length(e2);
    }

public:
    //Adds a single pair info, merging histograms if there's already one
    void Add(EdgeId e1, EdgeId e2, Point point) {
        SwapConj(e1, e2, point);
        Hist& histogram = storage_[e1][e2];
        auto iterator_to_point = histogram.find(point);

        if (iterator_to_point != histogram.end())
            MergeData(e1, e2, *iterator_to_point, point);
        else
            InsertPoint(e1, e2, histogram, point);
    }

    //Adds a whole histogram, merging histograms if there's already one
    void AddMany(EdgeId e1, EdgeId e2, const Histogram& hist) {
        bool swapped = SwapConj(e1, e2);
        Hist& histogram = storage_[e1][e2];
        for (auto point : hist) {
            if (swapped)
                point.d += CalcOffset(e1, e2);
            auto iterator_to_point = histogram.find(point);
            if (iterator_to_point != histogram.end())
                MergeData(e1, e2, *iterator_to_point, point);
            else
                InsertPoint(e1, e2, histogram, point);
        }
    }

private:

    static bool IsSymmetric(EdgeId e1, EdgeId e2, Point point) {
        return (e1 == e2) && math::eq(point.d, 0.f);
    }

    // modifying the histogram
    void InsertPoint(EdgeId e1, EdgeId e2, Hist& histogram, Point point) {
        // first backwards
        if (!IsSymmetric(e1, e2, point)) {
            storage_[e2][e1].insert(-point);
            ++size_;
        }

        histogram.insert(point);
        ++size_;
    }

    void UpdateSinglePoint(Hist& hist, typename Hist::iterator point_to_update, const Point& new_point) {
        typename Hist::iterator after_removed = hist.erase(point_to_update);
        hist.insert(after_removed, new_point);
    }

    void MergeData(EdgeId e1, EdgeId e2, const Point& point_to_update, const Point& point_to_add) {
        Hist& reversed = storage_[e2][e1];
        UpdateSinglePoint(reversed, reversed.find(-point_to_update), -(point_to_update + point_to_add));

        Hist& straight = storage_[e1][e2];
        UpdateSinglePoint(straight, straight.find(point_to_update), point_to_update + point_to_add);
    }

    void MergeData(Hist& hist, typename Hist::iterator to_update, const Point& point_to_add) {
        UpdateSinglePoint(hist, to_update, *to_update + point_to_add);
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
        typedef typename Hist::iterator hist_iterator;
        for (auto I = map_to_add.begin(), E = map_to_add.end(); I != E; ++I) {
            Hist &hist_exists = map[I->first];
            const auto& hist_to_add = I->second;

            for (auto p_it = hist_to_add.begin(), E = hist_to_add.end(); p_it != E; ++p_it) {
                Point new_point = *p_it;
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

    //Removes the specific entry
    // Returns the number of deleted entries
    size_t Remove(EdgeId e1, EdgeId e2, Point point) {
        return RemoveSingle(e1, e2, point) + RemoveSingle(e2, e1, -point);
    }
    // Removes the whole histogram
    // Returns the number of deleted entries
    size_t Remove(EdgeId e1, EdgeId e2) {
        return RemoveSingle(e1, e2) + RemoveSingle(e2, e1);
    }

private:
    //TODO: remove duplicode
    size_t RemoveSingle(EdgeId e1, EdgeId e2, Point point) {
        auto i1 = storage_.find(e1);
        if (i1 != storage_.end()) {
            auto &map = i1->second;
            auto i2 = map.find(e2);
            if (i2 != map.end()) {
                Hist& hist = i2->second;
                if (hist.erase(point))
                    --size_;
                if (hist.empty())
                    storage_.erase(i1);
                return 1;
            }
        }
        return 0;
    }

    size_t RemoveSingle(EdgeId e1, EdgeId e2) {
        auto i1 = storage_.find(e1);
        if (i1 != storage_.end()) {
            auto &map = i1->second;
            auto i2 = map.find(e2);
            if (i2 != map.end()) {
                Hist& hist = i2->second;
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
            if (edge != e2)
                this->Remove(e2, edge);
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

    /*inline const InnerMap& RawGet(EdgeId edge) const {
        VERIFY(storage_.count(edge) > 0);
        return storage_.find(edge)->second;
    }*/

    inline ImplIterator data_end() const {
        return storage_.end();
    }

    EdgeHist<> begin() const;

    EdgeHist<false> raw_begin() const;

    EdgeHist<> end() const;

    EdgeHist<false> raw_end() const;

    // Returns a proxy map to neighboring edges
    inline EdgeProxy<> Get(EdgeId id) const {
        return EdgeProxy<>(*this, id);
    }

    // Returns a proxy map to neighboring edges
    inline EdgeProxy<false> RawGet(EdgeId id) const {
        return EdgeProxy<false>(*this, id);
    }

    // Operator alias
    inline EdgeProxy<> operator[](EdgeId id) const {
        //TODO: what if doesn't contain?
        return Get(id);
    }

private:
    const Histogram& GetImpl(EdgeId e1, EdgeId e2) const {
        //TODO: make it stable
        return storage_.find(e1)->second.find(e2)->second;
    }
public:

    // Returns a proxy map to points
    HistProxy<> Get(EdgeId e1, EdgeId e2) const {
        int offset = CalcOffset(e1, e2);
        if (SwapConj(e1, e2)) {
            return HistProxy<>(GetImpl(e1, e2), 0);
        } else {
            return HistProxy<>(GetImpl(e1, e2), offset);
        }
    }

    // Operator alias
    inline HistProxy<> operator[](EdgePair p) const {
        //TODO: what if doesn't contain?
        return Get(p.first, p.second);
    }

    bool contains(EdgeId edge) const {
        return storage_.count(edge) + storage_.count(graph_.conjugate(edge)) > 0;
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

    // Prints the contents of index
    void PrintAll() const {
        size_t size = 0;
        for (const auto i: storage_) {
            EdgeId e1 = i.first.first; EdgeId e2 = i.first.second;
            const auto& histogram = i.second;
            size += histogram.size();
            INFO("Histogram for edges "
                 << graph_.int_id(e1) << " "
                 << graph_.int_id(e2));
            for (const auto& point : histogram) {
                INFO("    Entry " << point.str());
            }
        }
        VERIFY_MSG(this->size() == size, "Size " << size << " must have been equal to " << this->size());
    }

private:
    size_t size_;
    const Graph& graph_;
    StorageMap storage_;
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

template<class Graph>
struct PairedInfoIndicesT {
    typedef PairedInfoIndexT<Graph> IndexT;
    std::vector<IndexT> data_;

    PairedInfoIndicesT(const Graph& graph, size_t lib_num) {
        data_.reserve(lib_num);
        for (size_t i = 0; i < lib_num; ++i)
            data_.emplace_back(graph);
    }

    void Init() { for (auto& it : data_) it.Init(); }

    void Clear() { for (auto& it : data_) it.Clear(); }

    IndexT& operator[](size_t i) { return data_[i]; }

    const IndexT& operator[](size_t i) const { return data_[i]; }

    size_t size() const { return data_.size(); }
};

template<class Graph>
using UnclusteredPairedInfoIndicesT = std::vector<UnclusteredPairedInfoIndexT<Graph> >;

template<typename K, typename V>
using unordered_map = std::unordered_map<K, V>; //Two-parameters wrapper
template<class Graph>
using PairedInfoBuffer = PairedIndex<Graph, RawHistogram, unordered_map>;


}

}
