//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/paths/mapping_path.hpp"

class Sequence;

namespace io {
class SingleRead;
}

namespace debruijn_graph {

template<class Graph>
class SequenceMapper {
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef RtSeq Kmer;

    virtual ~SequenceMapper() = default;

    virtual omnigraph::MappingPath<EdgeId> MapSequence(const Sequence &sequence,
                                                       bool only_simple = false) const = 0;

    virtual omnigraph::MappingPath<EdgeId> MapRead(const io::SingleRead &read,
                                                   bool only_simple = false) const = 0;
};

}
