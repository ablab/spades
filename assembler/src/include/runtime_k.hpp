//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * runtime_map.hpp
 *
 *  Created on: Jun 21, 2012
 *      Author: andrey
 */

#ifndef RUNTIME_K_HPP_
#define RUNTIME_K_HPP_

#include "sequence/sequence.hpp"
#include "sequence/seq.hpp"
#include "sequence/simple_seq.hpp"
#include "sequence/rtseq.hpp"


#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "k_range.hpp"

namespace runtime_k {

#define T_SIZE sizeof(seq_element_type)

#define GET_T_ELEMENTS_NUMBER(value) ((value - 1) / (T_SIZE << 2) + 1)

#define GET_K_BY_TS(value) (value * (T_SIZE << 2))

#define GET_UPPER_BOUND(value) GET_K_BY_TS( GET_T_ELEMENTS_NUMBER(value) )


const size_t UPPER_BOUND = GET_UPPER_BOUND(MAX_K); //((MAX_K - 1) / (sizeof(seq_element_type) << 2) + 1) * (sizeof(seq_element_type) << 2);

const size_t MAX_TS = GET_T_ELEMENTS_NUMBER(MAX_K);

const size_t MIN_TS = GET_T_ELEMENTS_NUMBER(MIN_K);


typedef RuntimeSeq<UPPER_BOUND> RtSeq;


//Basic types and sequence <---> kmer functions
template <size_t size_>
class TypeContainerImpl {
public:
    typedef SimpleSeq<size_> Kmer;

    typedef unordered_set<Kmer, typename Kmer::hash, typename Kmer::equal_to> set_type;

    typedef std::vector<Kmer> vector_type;

    static Kmer from_sequence(const RtSeq& seq) {
        return seq.get_sseq<size_>();
    }

    static RtSeq to_sequence(const Kmer& kmer, size_t k = size_) {
        return RtSeq(kmer, k);
    }
};


template <size_t size_, typename Value>
class TypeValueContainerImpl: public TypeContainerImpl<size_> {

public:
    typedef TypeContainerImpl<size_> base;

    typedef typename base::Kmer Kmer;

    typedef typename base::set_type set_type;

    typedef unordered_map<Kmer, Value, typename Kmer::hash, typename Kmer::equal_to> map_type;

};

} /* namespace runtime_k */

#endif /* RUNTIME_K_HPP_ */
