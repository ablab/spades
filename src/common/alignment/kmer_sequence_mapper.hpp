//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "sequence_mapper.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include <memory>
#include <string>

namespace debruijn_graph {

template<class Graph>
class EdgeIndex;
}

namespace alignment {
  
class ShortKMerReadMapper: public debruijn_graph::AbstractSequenceMapper<debruijn_graph::Graph> {
    typedef typename debruijn_graph::EdgeId EdgeId;
    using EdgeIndex = debruijn_graph::EdgeIndex<debruijn_graph::Graph>;
public:
    explicit ShortKMerReadMapper(const debruijn_graph::Graph& g,
                                 const std::string &workdir,
                                 unsigned k = 31, unsigned min_occ = 2)
            : debruijn_graph::AbstractSequenceMapper<debruijn_graph::Graph>(g),
              k_(k), min_occ_(min_occ) {
        Init(workdir);
    }

    omnigraph::MappingPath<EdgeId> MapSequence(const Sequence &sequence,
                                               bool only_simple = false) const override;

private:
    void Init(const std::string &workdir);

    std::unique_ptr<EdgeIndex> index_;
    unsigned k_;
    unsigned min_occ_;
};

}
