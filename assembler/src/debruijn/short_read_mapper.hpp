/*
 * short_read_mapper.hpp
 *
 *  Created on: Dec 4, 2013
 *      Author: andrey
 */

//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once


#include "sequence_mapper.hpp"
#include "pacbio/pac_index.hpp"

namespace debruijn_graph {

//  boost::optional<int> GetISFromLongEdge(const Sequence &left, const Sequence &right, size_t is, size_t edge_length_threshold) const {
//      auto pos_left = GetLastKmerPos(left);
//      auto pos_right = GetFirstKmerPos(right);
//      if (pos_left.second == -1u || pos_right.second == -1u || pos_left.first != pos_right.first || g_.length(pos_left.first) < edge_length_threshold) {
//          return boost::none;
//      }
//
//      return boost::optional<int>((int) (pos_right.second - pos_left.second - k_ - is + left.size() + right.size()));
//  }
//
//    template<class PairedRead>
//    pair<bool, int> GetISFromLongEdge(const PairedRead& read, size_t edge_length_threshold) const {
//        auto pos_left = GetLastKmerPos(read.first.sequence());
//        auto pos_right = GetFirstKmerPos(read.second.sequence());
//        if (pos_left.second == -1u || pos_right.second == -1u || pos_left.first != pos_right.first || g_.length(pos_left.first) < edge_length_threshold) {
//          return make_pair(false, 0);
//        }
//
//        return make_pair(true, (int) (pos_right.second - pos_left.second - small_k_ - read.insert_size() + left.size() + right.size()));
//    }
//
//  pair<EdgeId, size_t> GetFirstKmerPos(const Sequence &sequence) const {
//    if (sequence.size() < k_) {
//      return make_pair(EdgeId(0), -1u);
//    }
//
//    Kmer left = sequence.start<Kmer>(k_);
//    left = kmer_mapper_.Substitute(left);
//
//    return index_.get(left);
//  }
//
//  pair<EdgeId, size_t> GetLastKmerPos(const Sequence &sequence) const {
//    if (sequence.size() < k_) {
//      return make_pair(EdgeId(0), -1u);
//    }
//
//    Kmer right = sequence.end<Kmer>(k_);
//    right = kmer_mapper_.Substitute(right);
//
//    return index_.get(right);
//  }
//    pair<EdgeId, size_t> GetKmerPos(const Kmer& kmer) const {
//        VERIFY(kmer.size() == k_);
//        if (sequence.size() < small_k_) {
//          return make_pair(EdgeId(0), -1u);
//        }
//
//        runtime_k::RtSeq left = sequence.start<runtime_k::RtSeq>(small_k_);
//        return index_.GetUniqueKmerPos(left);
//    }
//
//    pair<EdgeId, size_t> GetLastKmerPos(const Sequence &sequence) const {
//        if (sequence.size() < small_k_) {
//          return make_pair(EdgeId(0), -1u);
//        }
//
//        runtime_k::RtSeq right = sequence.end<runtime_k::RtSeq>(small_k_);
//        return index_.GetUniqueKmerPos(right);
//    }
  
template<class Graph>
class SensitiveReadMapper: public SequenceMapper<Graph> {
    using SequenceMapper<Graph>::g_;
private:

    size_t small_k_;
    pacbio::PacBioMappingIndex<Graph> index_;

public:

    typedef typename SequenceMapper<Graph>::Kmer Kmer;

    SensitiveReadMapper(const Graph& g, size_t k, size_t graph_k) :
        SequenceMapper<Graph>(g), small_k_(k), index_(g, small_k_, graph_k, false) 
    {
    }

    MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
        return index_.GetShortReadAlignment(sequence);
    }

    pair<EdgeId, size_t> GetKmerPos(const Kmer& kmer) const {
        VERIFY(kmer.size() == small_k_);

        return index_.GetUniqueKmerPos(kmer);
    }

    size_t KmerSize() const {
        return small_k_;
    }

};


template<class graph_pack>
std::shared_ptr<SequenceMapper<typename graph_pack::graph_t>> ChooseProperMapper(const graph_pack& gp, size_t read_length) {
    if (read_length > gp.k_value) {
        INFO("Read length = " << read_length << ", selecting usual mapper");
        return MapperInstance(gp);
    }
    else {
        INFO("Read length = " << read_length << ", selecting short read mapper");
        return std::make_shared<SensitiveReadMapper<typename graph_pack::graph_t> >(gp.g, read_length/ 3, gp.k_value);
    }
}

}

