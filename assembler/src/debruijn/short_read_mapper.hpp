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

template<class Graph>
class ShortReadMapper: public SequenceMapper<Graph> {
    using SequenceMapper<Graph>::g_;

private:
    pacbio::PacBioMappingIndex<Graph> index_;

    size_t small_k_;

public:
    ShortReadMapper(const Graph& g, size_t k, size_t graph_k) :
        SequenceMapper<Graph>(g), index_(g, k, graph_k, false), small_k_(k)
    {
    }

    ~ShortReadMapper() {
    }

    MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
        return index_.GetShortReadAlignment(sequence);
    }

    MappingPath<EdgeId> MapRead(const io::SingleRead &read) const {
        VERIFY(read.IsValid());
        return MapSequence(read.sequence());
    }

    pair<EdgeId, size_t> GetFirstKmerPos(const Sequence &sequence) const {
        if (sequence.size() < small_k_) {
          return make_pair(EdgeId(0), -1u);
        }

        runtime_k::RtSeq left = sequence.start<runtime_k::RtSeq>(small_k_);
        return index_.GetUniqueKmerPos(left);
    }

    pair<EdgeId, size_t> GetLastKmerPos(const Sequence &sequence) const {
        if (sequence.size() < small_k_) {
          return make_pair(EdgeId(0), -1u);
        }

        runtime_k::RtSeq right = sequence.end<runtime_k::RtSeq>(small_k_);
        return index_.GetUniqueKmerPos(right);
    }

    pair<bool, int> GetISFromLongEdge(const Sequence &left, const Sequence &right, size_t is, size_t edge_length_threshold) const {
        auto pos_left = GetLastKmerPos(left);
        auto pos_right = GetFirstKmerPos(right);
        if (pos_left.second == -1u || pos_right.second == -1u || pos_left.first != pos_right.first || g_.length(pos_left.first) < edge_length_threshold) {
          return make_pair(false, 0);
        }

        return make_pair(true, (int) (pos_right.second - pos_left.second - small_k_ - is + left.size() + right.size()));
    }
};


template<class graph_pack>
class MapperFactory {

private:
    const graph_pack& gp_;

public:
    typedef SequenceMapper<typename graph_pack::graph_t> SequenceMapperT;

    MapperFactory(const graph_pack& gp):gp_(gp) {
    }

    std::shared_ptr<SequenceMapperT> GetSequenceMapper(size_t read_length) {
        if (read_length > gp_.k_value) {
            INFO("Read length = " << read_length << ", selecting usual mapper");
            return std::make_shared<NewExtendedSequenceMapper<typename graph_pack::graph_t, typename graph_pack::index_t> >(gp_.g, gp_.index, gp_.kmer_mapper);
        }
        else {
            INFO("Read length = " << read_length << ", selecting short read mapper");
            return std::make_shared<ShortReadMapper<typename graph_pack::graph_t> >(gp_.g, read_length/ 3, gp_.k_value);
        }
    }

};

}

