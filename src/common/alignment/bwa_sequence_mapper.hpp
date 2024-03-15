//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "bwa_index.hpp"
#include "sequence_mapper.hpp"

#include "assembly_graph/paths/mapping_path.hpp"

namespace alignment {
  
template<class Graph>
class BWAReadMapper: public debruijn_graph::AbstractSequenceMapper<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    using debruijn_graph::AbstractSequenceMapper<Graph>::g_;
public:
    explicit BWAReadMapper(const Graph& g,
                           BWAIndex::AlignmentMode mode = BWAIndex::AlignmentMode::Default,
                           BWAIndex::RetainAlignments retain = BWAIndex::RetainAlignments::Default)
            : debruijn_graph::AbstractSequenceMapper<Graph>(g),
              index_(g, mode, retain) {}

    omnigraph::MappingPath<EdgeId> MapSequence(const Sequence &sequence,
                                               bool only_simple = false) const override {
        return index_.AlignSequence(sequence, only_simple);
    }

    BWAIndex index_;
};

}
