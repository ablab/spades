//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "assembly_graph/debruijn_graph.hpp"
#include "graph_support/basic_edge_conditions.hpp"

namespace debruijn_graph {

class ForbiddenPatternCondition : public EdgeCondition<Graph> {
    typedef EdgeCondition<Graph> base;

    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    Sequence pattern_;
    size_t max_offset_;
public:
    ForbiddenPatternCondition(const Graph& g, Sequence pattern, size_t max_offset) : base(g), pattern_(pattern), max_offset_(max_offset) {
    }

    /*virtual*/ bool Check(EdgeId e) const {
        Sequence nucls = this->g().EdgeNucls(e);
        for(size_t i = 0; i < max_offset_ && i + pattern_.size() < nucls.size(); i++) {
            if(nucls.Subseq(i, i + pattern_.size()) == pattern_ || (!nucls).Subseq(i, i + pattern_.size()) == pattern_) {
                return false;
            }
        }
        return true;
    }

};

}
