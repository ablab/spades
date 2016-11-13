//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * sequencem_mapping_logger.h
 *
 *  Created on: Nov 27, 2012
 *      Author: alex
 */

#ifndef KMER_MAPPER_LOGGER_H_
#define KMER_MAPPER_LOGGER_H_

#include "sequence/sequence.hpp"
#include "common/assembly_graph/core/action_handlers.hpp"
#include "utils/standard_base.hpp"

namespace debruijn {

template<class Graph>
class KmerMapperLogger : public omnigraph::GraphActionHandler<Graph> {
public:
    typedef pair<Sequence, Sequence> MappedSeq;
    typedef typename Graph::EdgeId EdgeId;

    KmerMapperLogger(Graph& graph) : GraphActionHandler<Graph>(graph, "KmerMapperLogger") {}
    virtual ~KmerMapperLogger() {}

    virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
        log_.push_back(MappedSeq(this->g().EdgeNucls(edge1), this->g().EdgeNucls(edge2)));
    }

    const vector<MappedSeq>& log() const {
        return log_;
    }

    vector<MappedSeq> log_;
};

} /* namespace debruijn */
#endif /* KMER_MAPPER_LOGGER_H_ */
