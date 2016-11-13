//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef RUNTIME_K_HPP_
#define RUNTIME_K_HPP_

#include "sequence.hpp"
#include "seq.hpp"
#include "simple_seq.hpp"
#include "rtseq.hpp"

#include "k_range.hpp"

namespace runtime_k {

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

const size_t UPPER_BOUND = get_upper_bound(MAX_K); //((MAX_K - 1) / (sizeof(seq_element_type) << 2) + 1) * (sizeof(seq_element_type) << 2);

const size_t MAX_TS = get_t_elements_number(MAX_K);

const size_t MIN_TS = get_t_elements_number(MIN_K);


typedef RuntimeSeq<UPPER_BOUND> RtSeq;

} /* namespace runtime_k */

#endif /* RUNTIME_K_HPP_ */
