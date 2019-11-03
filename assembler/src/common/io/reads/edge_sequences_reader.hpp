//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"

namespace io {

class EdgeSequencesStream {
 public:
    typedef SingleReadSeq ReadT;

    explicit EdgeSequencesStream(const debruijn_graph::Graph &g) : graph_(g), edge_iterator_(g.e_begin<true>()) {}

    bool is_open() const {
        return true;
    }

    bool eof() {
        return edge_iterator_ == graph_.e_end<true>();
    }

    EdgeSequencesStream &operator>>(SingleReadSeq &singleread) {
        singleread = SingleReadSeq(graph_.EdgeNucls(*edge_iterator_));
        ++edge_iterator_;
        return *this;
    }

    void close() const {}

    void reset() {
        edge_iterator_ = graph_.e_begin<true>();
    }

 private:
    const debruijn_graph::Graph &graph_;
    decltype(graph_.e_begin<true>()) edge_iterator_;
};
}
