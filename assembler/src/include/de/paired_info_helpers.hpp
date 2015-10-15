//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/optional.hpp>
#include "paired_info.hpp"

namespace omnigraph {

namespace de {

template<typename Index>
class EdgePairIterator :
        public boost::iterator_facade<EdgePairIterator<Index>,
                typename Index::FullHistProxy,
                boost::forward_traversal_tag,
                typename Index::FullHistProxy>
{
    typedef typename Index::ImplIterator OuterIterator;
    typedef boost::optional<typename Index::EdgeIterator> InnerIterator;

protected:
    //They're not intended to be constucted explicitly, only via begin/end.
    EdgePairIterator(const Index& index, OuterIterator i, InnerIterator j, bool conj = false)
        : index_(index), i_(i), j_(j), conj_(conj)
    {}

public:
    void increment() {
        ++(*j_);
        if (j_ == index_.Get(i_->first).end()) {
            ++i_;
            if (!conj_ && i_ == index_.data_end()) {
                conj_ = true;
                i_ = index_.data_begin();
            }
            auto edge = conj_ ? index_.graph().conjugate(i_->first) : i_->first;
            j_ = (i_ == index_.data_end()) ? InnerIterator() : index_.Get(edge).begin();
        }
    }

    typename Index::FullHistProxy dereference() const {
        return (*j_)->second;
    }

    bool equal(const EdgePairIterator &other) const {
        return (j_ == other.j_) && (conj_ == other.conj_);
    }

    typename Index::EdgeId first() const {
        if (conj_)
            return (*j_)->first;
        return i_->first;
    }

    typename Index::EdgeId second() const {
        if (conj_)
            return index_.graph().conjugate(i_->first);
        return (*j_)->first;
    }

    static EdgePairIterator begin(const Index& index) {
        InnerIterator inner;
        if (index.size())
            inner = index.Get(index.data_begin()->first).begin();
        return EdgePairIterator(index, index.data_begin(), inner, !inner);
    }
    static EdgePairIterator end(const Index& index) {
        return EdgePairIterator(index, index.data_end(), InnerIterator(), true);
    }
private:
    const Index &index_;
    OuterIterator i_;
    InnerIterator j_;
    bool conj_;
};

template<typename Storage>
inline EdgePairIterator<Storage> pair_begin(const Storage &s) {
    return EdgePairIterator<Storage>::begin(s);
}

template<typename Storage>
inline EdgePairIterator<Storage> pair_end(const Storage &s) {
    return EdgePairIterator<Storage>::end(s);
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
