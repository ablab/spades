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
#include "read_stream.hpp"
#include "single_read.hpp"
#include "header_naming.hpp"

namespace io {

class EdgeSequencesStream {
 public:
    typedef SingleRead ReadT;

    explicit EdgeSequencesStream(const debruijn_graph::Graph &g) : g(g), edge_iterator(g.e_begin<true>()) {}

    bool is_open() {
        return true;
    }

    bool eof() {
        return edge_iterator != g.e_end<true>();
    }

    EdgeSequencesStream &operator>>(SingleRead &singleread) {
        std::string edge_nucl = g.EdgeNucls(*edge_iterator).str();
        singleread = SingleRead(MakeContigId(id_++, edge_nucl.size(), g.coverage(*edge_iterator)), edge_nucl);
        edge_iterator++;
        return *this;
    }

    void close() {}

    void reset() {}

 private:
    const debruijn_graph::Graph &g;
    size_t id_;
    decltype(g.e_begin<true>()) edge_iterator;
};
}
