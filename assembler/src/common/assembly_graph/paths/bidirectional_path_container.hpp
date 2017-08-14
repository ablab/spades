//
// Created by andrey on 14.08.17.
//

#pragma once

#include "bidirectional_path.hpp"
#include "modules/path_extend/path_filter.hpp"
#include <vector>
#include <set>
#include <map>

namespace path_extend {

using namespace std;

typedef std::pair<BidirectionalPath*, BidirectionalPath*> PathPair;

class PathComparator {
public:
    bool operator()(const BidirectionalPath& p1, const BidirectionalPath& p2) const {
        return p1.GetId() < p2.GetId();
    }

    bool operator()(const BidirectionalPath* p1, const BidirectionalPath* p2) const {
        return p1->GetId() < p2->GetId();
    }
};

typedef set<BidirectionalPath*, PathComparator> BidirectionalPathSet;

template<class Value>
using BidirectionalPathMap = map<BidirectionalPath*, Value, PathComparator>;

typedef multiset<BidirectionalPath *, PathComparator> BidirectionalPathMultiset;

class PathContainer {
public:

    typedef vector<PathPair> PathContainerT;

    class Iterator : public PathContainerT::iterator {
    public:
        Iterator(const PathContainerT::iterator& iter)
            : PathContainerT::iterator(iter) {
        }
        BidirectionalPath* get() const {
            return this->operator *().first;
        }
        BidirectionalPath* getConjugate() const {
            return this->operator *().second;
        }
    };

    class ConstIterator : public PathContainerT::const_iterator {
    public:
        ConstIterator(const PathContainerT::const_iterator& iter)
            : PathContainerT::const_iterator(iter) {
        }

        ConstIterator(const PathContainer::Iterator& iter)
            : PathContainerT::const_iterator(iter) {
        }

        BidirectionalPath* get() const {
            return this->operator *().first;
        }
        BidirectionalPath* getConjugate() const {
            return this->operator *().second;
        }
    };

    PathContainer() {
    }


    PathContainer(const PathContainer&) = delete;
    PathContainer& operator=(const PathContainer&) = delete;

    PathContainer(PathContainer&&) = default;
    PathContainer& operator=(PathContainer&&) = default;

    PathContainer(ConstIterator begin, ConstIterator end) {
        DeleteAllPaths();
        for (ConstIterator it = begin; it != end; ++it) {
            AddPair(new BidirectionalPath(*it.get()), new BidirectionalPath(*it.getConjugate()));
        }
    }

    BidirectionalPath& operator[](size_t index) const {
        return *(data_[index].first);
    }

    BidirectionalPath* Get(size_t index) const {
        return data_[index].first;
    }

    BidirectionalPath* GetConjugate(size_t index) const {
        return data_[index].second;
    }

    void Swap(size_t index) {
        std::swap(data_[index].first, data_[index].second);
    }

    void DeleteAllPaths() {
        for (size_t i = 0; i < data_.size(); ++i) {
            DeletePathPair(data_[i]);
        }
        clear();
    }

    ~PathContainer() {
        DeleteAllPaths();
    }

    size_t size() const {
        return data_.size();
    }

    void clear() {
        data_.clear();
    }

    void reserve(size_t size) {
        data_.reserve(size);
    }

    bool AddPair(BidirectionalPath* p, BidirectionalPath* cp) {
        p->SetConjPath(cp);
        cp->SetConjPath(p);
        p->Subscribe(cp);
        cp->Subscribe(p);
        data_.push_back(std::make_pair(p, cp));
        return true;
    }

    void SortByLength(bool desc = true) {
        std::stable_sort(data_.begin(), data_.end(), [=](const PathPair& p1, const PathPair& p2) {
            if (p1.first->Empty() || p2.first->Empty() || p1.first->Length() != p2.first->Length()) {
                return desc ? p1.first->Length() > p2.first->Length()
                            : p1.first->Length() < p2.first->Length();
            }
            const Graph& g = p1.first->graph();
            return g.int_id(p1.first->Front()) < g.int_id(p2.first->Front());
        });
    }

    Iterator begin() {
        return Iterator(data_.begin());
    }

    Iterator end() {
        return Iterator(data_.end());
    }

    ConstIterator begin() const {
        return ConstIterator(data_.begin());
    }

    ConstIterator end() const {
        return ConstIterator(data_.end());
    }

    Iterator erase(Iterator iter) {
        return Iterator(data_.erase(iter));
    }

    void print() const {
        for (size_t i = 0; i < size(); ++i) {
            Get(i)->PrintDEBUG();
            GetConjugate(i)->PrintDEBUG();
        }
    }

    void FilterPaths(const func::TypedPredicate<const BidirectionalPath&>& pred) {
        DEBUG("Removing empty paths");
        for (auto &pp : data_) {
            if (pred(*pp.first)) {
                VERIFY(pred(*pp.second)); //do we need it?
                DeletePathPair(pp);
            }
        }

        const PathPair empty_pp(nullptr, nullptr);
        data_.erase(std::remove(data_.begin(), data_.end(), empty_pp), data_.end());
        DEBUG("Empty paths removed");
    }

    void FilterEmptyPaths() {
        FilterPaths(EmptyPathCondition());
    }

private:

    void DeletePathPair(PathPair &pp) {
        delete pp.first;
        pp.first = nullptr;
        delete pp.second;
        pp.second = nullptr;
    }

    vector<PathPair> data_;

protected:
    DECL_LOGGER("BidirectionalPath");

};

}