//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __ITERATOR_RANGE_H__
#define __ITERATOR_RANGE_H__

#include <utility>
#include <iterator>

namespace adt {

template<typename IteratorT>
class iterator_range {
    IteratorT begin_iterator, end_iterator;

public:
    template<typename Container>
    iterator_range(Container &&c)
    //TODO: Consider ADL/non-member begin/end calls.
            : begin_iterator(c.begin()), end_iterator(c.end()) { }

    iterator_range(IteratorT begin_iterator, IteratorT end_iterator)
            : begin_iterator(std::move(begin_iterator)),
              end_iterator(std::move(end_iterator)) { }

    IteratorT begin() const { return begin_iterator; }

    IteratorT end() const { return end_iterator; }
};

template<class T>
iterator_range<T> make_range(T x, T y) {
    return iterator_range<T>(std::move(x), std::move(y));
}

template<typename T>
iterator_range<T> make_range(std::pair<T, T> p) {
    return iterator_range<T>(std::move(p.first), std::move(p.second));
}

template<typename T>
iterator_range<decltype(begin(std::declval<T>()))> drop_begin(T &&t, int n) {
    return make_range(std::next(begin(t), n), end(t));
}
}

#endif
