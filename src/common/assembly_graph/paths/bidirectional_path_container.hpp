//***************************************************************************
//* Copyright (c) 2017-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "bidirectional_path.hpp"
#include <memory>
#include <vector>
#include <set>
#include <map>

namespace path_extend {

class PathComparator {
public:
    bool operator()(const BidirectionalPath& p1, const BidirectionalPath& p2) const {
        return p1.GetId() < p2.GetId();
    }

    bool operator()(const BidirectionalPath* p1, const BidirectionalPath* p2) const {
        return p1->GetId() < p2->GetId();
    }
};

typedef std::set<BidirectionalPath*, PathComparator> BidirectionalPathSet;

template<class Value>
using BidirectionalPathMap = std::map<BidirectionalPath*, Value, PathComparator>;

typedef std::multiset<BidirectionalPath *, PathComparator> BidirectionalPathMultiset;

class PathContainer {
    typedef std::pair<std::unique_ptr<BidirectionalPath>,
                      std::unique_ptr<BidirectionalPath>> PathPair;
public:
    typedef std::vector<PathPair> PathContainerT;

    class Iterator : public PathContainerT::iterator {
    public:
        Iterator(const PathContainerT::iterator& iter)
            : PathContainerT::iterator(iter) {
        }
        BidirectionalPath& get() const { return *(this->operator *().first); }
        BidirectionalPath& getConjugate() const { return *(this->operator *().second); }
    };

    class ConstIterator : public PathContainerT::const_iterator {
    public:
        ConstIterator(const PathContainerT::const_iterator& iter)
            : PathContainerT::const_iterator(iter) { }

        ConstIterator(const PathContainer::Iterator& iter)
            : PathContainerT::const_iterator(PathContainerT::iterator(iter)) { }

        const BidirectionalPath& get() const { return *(this->operator *().first);  }
        const BidirectionalPath& getConjugate() const { return *(this->operator *().second); }
    };

    PathContainer() {}

    PathContainer(const PathContainer&) = delete;
    PathContainer& operator=(const PathContainer&) = delete;

    PathContainer(PathContainer&&) = default;
    PathContainer& operator=(PathContainer&&) = default;

    PathContainer(ConstIterator begin, ConstIterator end) {
        clear();
        for (ConstIterator it = begin; it != end; ++it) {
            AddPair(BidirectionalPath::clone(it.get()),
                    BidirectionalPath::clone(it.getConjugate()));
        }
    }

    BidirectionalPath& operator[](size_t index) {
        return *(data_[index].first);
    }
    const BidirectionalPath& operator[](size_t index) const {
        return *(data_[index].first);
    }

    BidirectionalPath& Get(size_t index) {
        return *data_[index].first;
    }
    const BidirectionalPath& Get(size_t index) const {
        return *data_[index].first;
    }

    BidirectionalPath& GetConjugate(size_t index) {
        return *data_[index].second;
    }
    const BidirectionalPath& GetConjugate(size_t index) const {
        return *data_[index].second;
    }

    void Swap(size_t index) {
        std::swap(data_[index].first, data_[index].second);
    }

    ~PathContainer() {}

    size_t size() const { return data_.size(); }
    void clear() { data_.clear(); }
    void reserve(size_t size) { data_.reserve(size); }

    std::pair<BidirectionalPath&, BidirectionalPath&>
    LinkPair(std::pair<BidirectionalPath&, BidirectionalPath&> ppair) const {
        ppair.first.SetConjPath(&ppair.second);
        ppair.second.SetConjPath(&ppair.first);
        ppair.first.Subscribe(ppair.second);
        ppair.second.Subscribe(ppair.first);

        return ppair;
    }
    
    // This guy acquires the ownership of paths
    std::pair<BidirectionalPath&, BidirectionalPath&>
    AddPair(std::unique_ptr<BidirectionalPath> p, std::unique_ptr<BidirectionalPath> cp) {
        data_.emplace_back(std::move(p), std::move(cp));
        auto &entry = data_.back();

        return LinkPair({ *entry.first, *entry.second });
    }

    std::pair<BidirectionalPath&, BidirectionalPath&>
    Add(std::unique_ptr<BidirectionalPath> p) {
        auto cp = (p->GetConjPath() ?
                   BidirectionalPath::clone(*p->GetConjPath()) : BidirectionalPath::clone_conjugate(p));
        return AddPair(std::move(p), std::move(cp));
    }

    template<typename... Args>
    std::pair<BidirectionalPath&, BidirectionalPath&>
    CreatePair(Args&&... args) {
        return Add(BidirectionalPath::create(std::forward<Args>(args)...));
    }

    template<typename... Args>
    BidirectionalPath&
    Create(Args&&... args) {
        return Add(BidirectionalPath::create(std::forward<Args>(args)...)).first;
    }

    void SortByLength(bool desc = true) {
        std::stable_sort(data_.begin(), data_.end(), [=](const PathPair& p1, const PathPair& p2) {
            if (p1.first->Empty() || p2.first->Empty() || p1.first->Length() != p2.first->Length()) {
                return desc ? p1.first->Length() > p2.first->Length()
                            : p1.first->Length() < p2.first->Length();
            }
            const debruijn_graph::Graph& g = p1.first->graph();
            return g.int_id(p1.first->Front()) < g.int_id(p2.first->Front());
        });
    }

    Iterator begin() { return Iterator(data_.begin()); }
    Iterator end() { return Iterator(data_.end()); }
    ConstIterator begin() const { return ConstIterator(data_.begin()); }
    ConstIterator end() const { return ConstIterator(data_.end()); }

    Iterator erase(ConstIterator iter) {
        return Iterator(data_.erase(iter));
    }

    void print() const {
        for (size_t i = 0; i < size(); ++i) {
            Get(i).PrintDEBUG();
            GetConjugate(i).PrintDEBUG();
        }
    }

    void FilterPaths(func::TypedPredicate<const BidirectionalPath&> pred) {
        DEBUG("Filtering paths based on predicate");
        for (auto &pp : data_) {
            if (pred(*pp.first)) {
                VERIFY(pred(*pp.second)); //do we need it?
                DeletePathPair(pp);
            }
        }

        const PathPair empty_pp(nullptr, nullptr);
        data_.erase(std::remove(data_.begin(), data_.end(), empty_pp), data_.end());
        DEBUG("Paths filtered");
    }

    void FilterEmptyPaths() {
        DEBUG("Removing empty paths");
        FilterPaths([](const BidirectionalPath &p) { return p.Empty(); });
    }

private:
    void DeletePathPair(PathPair &pp) {
        pp.first.reset();
        pp.second.reset();
    }

    PathContainerT data_;

protected:
    DECL_LOGGER("BidirectionalPath");

};

}
