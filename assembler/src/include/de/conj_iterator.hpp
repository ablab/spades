//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/iterator/iterator_facade.hpp>

namespace omnigraph {

namespace de {

/**
 * @brief Proxy for containers which are essentially splitted into two: the straight and the conjugate one
 *        (any of which can be empty).
 * @param C the underlying container type
 */
template<typename C>
class ConjProxy {
public:
    typedef C Container;

    /**
     * @brief Iterator for this splitted container.
     *        It automatically switches onto the conjugate half when finished the straight.
     */
    class Iterator :
            public boost::iterator_facade<Iterator, typename Container::const_reference, boost::bidirectional_traversal_tag> {
    public:
        typedef typename Container::const_iterator InnerIterator;

        Iterator(InnerIterator start_iter, InnerIterator stop_iter, InnerIterator jump_iter, bool conj)
                : iter_(start_iter), stop_iter_(stop_iter), jump_iter_(jump_iter), conj_(conj) { }

        /**
         * @brief Increments the iterator.
         * @detail The underlying iterator is incremented; when it reaches the `stop` position,
         *         it jumps to the `jump` position.
         */
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

        /**
         * @brief Returns the container const_iterator to the current element.
         */
        inline InnerIterator Iter() const {
            return iter_;
        }

        /**
         * @brief Returns if the iterator is on the conjugate half.
         */
        inline bool Conj() const {
            return conj_;
        }

    private:
        InnerIterator iter_, //the current position
                stop_iter_,  //when to stop and jump
                jump_iter_;  //where to jump
        bool conj_;
    };

    ConjProxy(const Container &cont, const Container &conj_cont) :
            cont_(cont),
            conj_cont_(conj_cont) { }

    /**
     * @brief Iteration always starts from the beginning of the straight half.
     */
    Iterator begin() const {
        auto conj = cont_.empty();
        auto start = conj ? conj_cont_.begin() : cont_.begin();
        return Iterator(start, cont_.end(), conj_cont_.begin(), conj);
    }

    /**
     * @brief Raw iterator should end right after the jumping, i.e. on the beginning
     *        of the conjugate half.
     */
    inline Iterator conj_begin() const {
        return Iterator(conj_cont_.begin(), cont_.end(), conj_cont_.begin(), true);
    }

    /**
     * @brief Full iterator ends on the end of the conjugate half.
     */
    inline Iterator end() const {
        return Iterator(conj_cont_.end(), cont_.end(), conj_cont_.begin(), true);
    }

    /**
     * @brief Returns the total size of both halves.
     */
    inline size_t size() const {
        return cont_.size() + conj_cont_.size();
    }

    /**
     * @brief Returns if both halves are empty.
     */
    inline bool empty() const {
        return cont_.empty() && conj_cont_.empty();
    }

private:
    const Container &cont_, &conj_cont_;
};

}

}
