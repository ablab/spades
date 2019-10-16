//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/**

* Reader<SingleRead> is the very base class that reads from one file
* through Parser object.
* Reader<PairedRead> is the class that reads data from two input
* files and gets PairedReads using this data and distance information.
*/

#pragma once

#include <common/assembly_graph/core/graph.hpp>

namespace io {

class EdgeSequencesStream {
 public:
    typedef SingleReadSeq ReadT;

    explicit EdgeSequencesStream(const debruijn_graph::Graph &g) : graph_(g), edge_iterator_(g.e_begin<true>()) {}

    bool is_open() {
        return true;
    }

    bool eof() {
        return edge_iterator_ == graph_.e_end<true>();
    }

    EdgeSequencesStream &operator>>(SingleReadSeq &singleread) {
        auto edge_nucl = graph_.EdgeNucls(*edge_iterator_);
        singleread = SingleReadSeq(edge_nucl);
        edge_iterator_++;
        return *this;
    }

    void close() {}

    void reset() {}

 private:
    const debruijn_graph::Graph &graph_;
    decltype(graph_.e_begin<true>()) edge_iterator_;
};
}
