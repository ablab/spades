//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "sequence_mapper.hpp"
#include "bwa_index.hpp"
#include "assembly_graph/paths/mapping_path.hpp"

namespace alignment {
  
template<class Graph>
class BWAReadMapper: public debruijn_graph::AbstractSequenceMapper<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    using debruijn_graph::AbstractSequenceMapper<Graph>::g_;
public:
    explicit BWAReadMapper(const Graph& g,
                           BWAIndex::AlignmentMode mode = BWAIndex::AlignmentMode::Default,
                           size_t length_cutoff = 0)
            : debruijn_graph::AbstractSequenceMapper<Graph>(g),
            index_(g, mode, length_cutoff) {}

    omnigraph::MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
        return index_.AlignSequence(sequence);
    }

    BWAIndex index_;
};

}
