/*
 * short_read_mapper.hpp
 *
 *  Created on: Dec 4, 2013
 *      Author: andrey
 */

//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once


#include "sequence_mapper.hpp"
#include "pacbio/pac_index.hpp"

namespace debruijn_graph {
  
template<class Graph>
class SensitiveReadMapper: public SequenceMapper<Graph> {
    using SequenceMapper<Graph>::g_;
private:

    size_t small_k_;

    static map<size_t, pacbio::PacBioMappingIndex<Graph>* > indices_;
    static size_t active_mappers_;

    pacbio::PacBioMappingIndex<Graph>* index_;

public:

    typedef typename SequenceMapper<Graph>::Kmer Kmer;

    SensitiveReadMapper(const Graph& g, size_t k, size_t graph_k) :
        SequenceMapper<Graph>(g), small_k_(k)
    {
        if (indices_.find(small_k_) == indices_.end()) {
            indices_.insert(make_pair(small_k_,
                                      new pacbio::PacBioMappingIndex<Graph>(g, small_k_, graph_k, false)));
        }
        index_ = indices_[small_k_];
        ++active_mappers_;
    }

    MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
        return index_->GetShortReadAlignment(sequence);
    }

    pair<EdgeId, size_t> GetKmerPos(const Kmer& kmer) const {
        VERIFY(kmer.size() == small_k_);

        return index_->GetUniqueKmerPos(kmer);
    }

    size_t KmerSize() const {
        return small_k_;
    }

    ~SensitiveReadMapper() {
        --active_mappers_;
    }

    static void EraseIndices() {
        if (active_mappers_ > 0) {
            WARN("There are still active mappers");
        }
        for (auto iter = indices_.begin(); iter != indices_.end(); ++iter) {
            delete iter->second;
        }
        indices_.clear();
    }

};

template<class Graph>
map<size_t, pacbio::PacBioMappingIndex<Graph>* > SensitiveReadMapper<Graph>::indices_;

template<class Graph>
size_t SensitiveReadMapper<Graph>::active_mappers_ = 0;


template<class graph_pack, class SequencingLib>
std::shared_ptr<SequenceMapper<typename graph_pack::graph_t>> ChooseProperMapper(const graph_pack& gp, const SequencingLib& library) {
    if (library.type() == io::LibraryType::MatePairs) {
        INFO("Mapping mate-pair library, selecting sensitive read mapper with k=" << cfg::get().sensitive_map.k);
        return std::shared_ptr<SequenceMapper<typename graph_pack::graph_t>> (new SensitiveReadMapper<typename graph_pack::graph_t>(gp.g, cfg::get().sensitive_map.k, gp.k_value));
    }

    size_t read_length = library.data().read_length; 
    if (read_length < gp.k_value && library.type() == io::LibraryType::PairedEnd) {
        INFO("Read length = " << read_length << ", selecting short read mapper");
        return std::shared_ptr<SequenceMapper<typename graph_pack::graph_t>> (new SensitiveReadMapper<typename graph_pack::graph_t> (gp.g, read_length/ 3, gp.k_value));
    }

    INFO("Selecting usual mapper");
    return MapperInstance(gp);
}

}

