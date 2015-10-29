//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/iterator/iterator_facade.hpp>

namespace omnigraph {

namespace de {

//Proxy for containers which are essentially splitted into two: the straight and the conjugate one.
template<typename C>
class ConjProxy {
public:
    typedef C Container;

    //Iterator for this splitted container.
    //It automatically switches onto the conjugate half when finished the straight.
    class Iterator :
            public boost::iterator_facade<Iterator, typename Container::const_reference, boost::bidirectional_traversal_tag> {
    public:
        typedef typename Container::const_iterator InnerIterator;

        Iterator(InnerIterator start_iter, InnerIterator stop_iter, InnerIterator jump_iter, bool conj)
                : iter_(start_iter), stop_iter_(stop_iter), jump_iter_(jump_iter), conj_(conj) { }

        void increment() {
            ++iter_;
            if (!conj_ && iter_ == stop_iter_) {
                conj_ = true;
                iter_ = jump_iter_;
            }
        }

        void decrement() {
            if (conj_ && iter_ == jump_iter_) {
                conj_ = false;
                iter_ = stop_iter_;
            }
            --iter_;
        }

        inline bool equal(const Iterator &other) const {
            return iter_ == other.iter_ && conj_ == other.conj_;
        }

        inline typename C::const_reference dereference() const {
            return *iter_;
        }

        inline InnerIterator Iter() const {
            return iter_;
        }

        inline bool Conj() const {
            return conj_;
        }

    private:
        InnerIterator iter_, stop_iter_, jump_iter_;
        bool conj_;
    };

    ConjProxy(const Container &cont, const Container &conj_cont) :
            cont_(cont),
            conj_cont_(conj_cont) { }

    Iterator begin() const {
        auto conj = cont_.empty();
        auto start = conj ? conj_cont_.begin() : cont_.begin();
        return Iterator(start, cont_.end(), conj_cont_.begin(), conj);
    }

    inline Iterator conj_begin() const {
        return Iterator(conj_cont_.begin(), cont_.end(), conj_cont_.begin(), true);
    }

    inline Iterator end() const {
        return Iterator(conj_cont_.end(), cont_.end(), conj_cont_.begin(), true);
    }

    inline size_t size() const {
        return cont_.size() + conj_cont_.size();
    }

    inline bool empty() const {
        return cont_.empty() && conj_cont_.empty();
    }

private:
    const Container &cont_, &conj_cont_;
};

}

}
