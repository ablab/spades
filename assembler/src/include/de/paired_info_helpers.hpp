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

template<typename Index>
class EdgePairIterator :
        public boost::iterator_facade<EdgePairIterator<Index>,
                typename Index::FullHistProxy,
                boost::forward_traversal_tag,
                typename Index::FullHistProxy>
{
    typedef typename Index::ImplIterator OuterIterator;
    typedef boost::optional<typename Index::InnerMap::const_iterator> InnerIterator;

protected:
    //They're not intended to be constucted explicitly, only via begin/end.
    EdgePairIterator(const Index& index, OuterIterator i, bool conj = false)
        : index_(index), i_(i), conj_(conj)
    {
        StartOver();
    }

public:
    void increment() {
        ++(*j_);
        if (j_ == i_->second.end()) { //Traversed all neighbours, jump to the next edge
            ++i_;
            if (!conj_ && i_ == index_.data_end()) {
                conj_ = true;
                i_ = index_.data_begin();
            }
            StartOver();
        }
    }

private:
    void StartOver() {
        if (i_ == index_.data_end()) {
            j_.reset();
        } else {
            j_ = i_->second.begin();
            /*
            InnerIterator k = j_;
            while (k != stop_j_)
                k->increment();
            INFO("CHECKED");
            */
        }
    }

public:

    typename Index::FullHistProxy dereference() const {
        return index_.Get(first(), second()); //TODO: optimize
    }

    bool equal(const EdgePairIterator &other) const {
        return (j_ == other.j_) && (conj_ == other.conj_);
    }

    typename Index::EdgeId first() const {
        if (conj_)
            return index_.graph().conjugate((*j_)->first);
        return i_->first;
    }

    typename Index::EdgeId second() const {
        if (conj_)
            return index_.graph().conjugate(i_->first);
        return (*j_)->first;
    }

    static EdgePairIterator begin(const Index& index) {
        return EdgePairIterator(index, index.data_begin(), !index.size());
    }
    static EdgePairIterator end(const Index& index) {
        return EdgePairIterator(index, index.data_end(), true);
    }
private:
    const Index &index_;
    OuterIterator i_;
    InnerIterator j_, stop_j_;
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
