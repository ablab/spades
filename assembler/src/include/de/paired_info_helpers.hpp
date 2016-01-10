//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "paired_info.hpp"
#include "boost/optional.hpp"

namespace omnigraph {

namespace de {

template<typename Index, bool full>
class EdgePairIterator :
        public boost::iterator_facade<EdgePairIterator<Index, full>,
                typename Index::FullHistProxy,
                boost::forward_traversal_tag,
                typename Index::FullHistProxy>
{
    typedef typename ConjProxy<typename Index::StorageMap>::Iterator OuterIterator;
    typedef boost::optional<typename Index::InnerMap::const_iterator> InnerIterator;

protected:
    //They're not intended to be constucted explicitly, only via begin/end.
    EdgePairIterator(const Index& index, OuterIterator i)
        : index_(index), i_(i)
    {
        StartOver();
    }

    bool FakePair() {
        auto ep = std::make_pair(i_->first, (*j_)->first);
        return ep > index_.ConjugatePair(ep);
    }

    inline void Skip() { //For a raw iterator, skip conjugate pairs
        while (!full && j_ && FakePair()) {
            IncImpl();
        }
    }

    void IncImpl() {
        ++(*j_);
        if (j_ == i_->second.end()) { //Traversed all neighbours, jump to the next edge
            ++i_;
            StartOver();
        }
    }

public:
    void increment() {
        IncImpl();
        Skip();
    }

private:
    void StartOver() {
        if (i_.Iter() == index_.data_end()) {
            j_.reset();
        } else {
            j_ = i_->second.begin();
            Skip();
        }
    }

public:

    typename Index::FullHistProxy dereference() const {
        return index_.Get(first(), second()); //TODO: optimize
    }

    bool equal(const EdgePairIterator &other) const {
        return (j_ == other.j_) && (i_.Conj() == other.i_.Conj());
    }

    typename Index::EdgeId first() const {
        if (i_.Conj())
            return index_.graph().conjugate((*j_)->first);
        return i_->first;
    }

    typename Index::EdgeId second() const {
        if (i_.Conj())
            return index_.graph().conjugate(i_->first);
        return (*j_)->first;
    }

    static EdgePairIterator begin(const Index& index) {
        auto i = OuterIterator(index.data_begin(), index.data_end(), index.data_begin(), !index.size());
        return EdgePairIterator(index, i);
    }

    static EdgePairIterator end(const Index& index) {
        auto stop = full ? index.data_end() : index.data_begin();
        auto i = OuterIterator(stop, index.data_end(), index.data_begin(), true);
        return EdgePairIterator(index, i);
    }

private:
    const Index &index_;
    OuterIterator i_;
    InnerIterator j_;
};

template<typename Storage>
inline EdgePairIterator<Storage, true> pair_begin(const Storage &s) {
    return EdgePairIterator<Storage, true>::begin(s);
}

template<typename Storage>
inline EdgePairIterator<Storage, true> pair_end(const Storage &s) {
    return EdgePairIterator<Storage, true>::end(s);
}

template<typename Storage>
inline EdgePairIterator<Storage, false> raw_pair_begin(const Storage &s) {
    return EdgePairIterator<Storage, false>::begin(s);
}

template<typename Storage>
inline EdgePairIterator<Storage, false> raw_pair_end(const Storage &s) {
    return EdgePairIterator<Storage, false>::end(s);
}

//Small wrapper for range-based loops
//Usage: for (auto i in PairsOf(index))
/*template <typename Storage>
class PairsOf {
public:
    EdgePairIterator<Storage> begin() const{
        return pair_begin(storage_);
    }

    EdgePairIterator<Storage> end() const{
        return pair_begin(storage_);
    }

    PairsOf(const Storage& storage)
            : storage_(storage) {}
private:
    const Storage& storage_;
};*/

}

}
