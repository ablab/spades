//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * seq_common.hpp
 *
 *  Created on: Jun 25, 2012
 *      Author: andrey
 */

#ifndef SEQ_COMMON_HPP_
#define SEQ_COMMON_HPP_

#include "k_range.hpp"

typedef u_int64_t seq_element_type;

constexpr size_t t_size(void) {
    return sizeof(seq_element_type);
}

constexpr size_t get_t_elements_number(size_t value) {
    return ((value - 1) / (t_size() << 2) + 1);
}

constexpr size_t get_k_by_ts(size_t value) {
    return (value * (t_size() << 2));
}

constexpr size_t get_upper_bound(size_t value) {
    return get_k_by_ts(get_t_elements_number(value));
}

const size_t UPPER_BOUND = get_upper_bound(runtime_k::MAX_K); //((MAX_K - 1) / (sizeof(seq_element_type) << 2) + 1) * (sizeof(seq_element_type) << 2);

const size_t MAX_TS = get_t_elements_number(runtime_k::MAX_K);

const size_t MIN_TS = get_t_elements_number(runtime_k::MIN_K);

#endif /* SEQ_COMMON_HPP_ */
