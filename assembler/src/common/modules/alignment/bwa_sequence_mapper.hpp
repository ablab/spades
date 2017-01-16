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
    BWAReadMapper(const Graph& g)
            : debruijn_graph::AbstractSequenceMapper<Graph>(g),
            index_(g) {}

    omnigraph::MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
        return index_.AlignSequence(sequence);
    }

    ~BWAReadMapper() {
    }

    BWAIndex index_;
};

}

