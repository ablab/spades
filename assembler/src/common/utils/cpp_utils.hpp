//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * cpp_utils.hpp
 *
 *  Created on: Nov 14, 2011
 *      Author: valery
 */

#pragma once

namespace utils {

// arrays
template<class T, size_t N>
size_t array_size(T (&/*arr*/)[N]) {
    return N;
}

template<class T, size_t N>
T *array_end(T (&arr)[N]) {
    return &arr[N];
}

template<size_t EXPECTED_SIZE, class T, size_t N>
void check_array_size(T (&/*arr*/)[N]) {
    static_assert(EXPECTED_SIZE == N, "Unexpected array size");
}

template<class T>
T identity_function(const T &t) {
    return t;
}

} // namespace utils
