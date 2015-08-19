//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * parallel_wrapper.hpp
 *
 *  Created on: Dec 14, 2013
 *      Author: anton
 */

#pragma once
#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

namespace parallel {
#ifdef USE_GLIBCXX_PARALLEL
    template <class RandomAccessIterator>
    void sort (RandomAccessIterator first, RandomAccessIterator last) {
        __gnu_parallel::sort(first, last);
    }
    template <class RandomAccessIterator, class Compare>
    void sort (RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
        __gnu_parallel::sort(first, last, comp);
    }
#else
    template <class RandomAccessIterator>
    void sort (RandomAccessIterator first, RandomAccessIterator last) {
        std::sort(first, last);
    }
    template <class RandomAccessIterator, class Compare>
    void sort (RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
        std::sort(first, last, comp);
    }
#endif
}

