//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef RUNTIME_K_HPP_
#define RUNTIME_K_HPP_

#include "data_structures/sequence/sequence.hpp"
#include "data_structures/sequence/seq.hpp"
#include "data_structures/sequence/simple_seq.hpp"
#include "data_structures/sequence/rtseq.hpp"

#include "k_range.hpp"

namespace runtime_k {

#define T_SIZE sizeof(seq_element_type)

#define GET_T_ELEMENTS_NUMBER(value) ((value - 1) / (T_SIZE << 2) + 1)

#define GET_K_BY_TS(value) (value * (T_SIZE << 2))

#define GET_UPPER_BOUND(value) GET_K_BY_TS(GET_T_ELEMENTS_NUMBER(value))


const size_t UPPER_BOUND = GET_UPPER_BOUND(MAX_K); //((MAX_K - 1) / (sizeof(seq_element_type) << 2) + 1) * (sizeof(seq_element_type) << 2);

const size_t MAX_TS = GET_T_ELEMENTS_NUMBER(MAX_K);

const size_t MIN_TS = GET_T_ELEMENTS_NUMBER(MIN_K);


typedef RuntimeSeq<UPPER_BOUND> RtSeq;

} /* namespace runtime_k */

#endif /* RUNTIME_K_HPP_ */
