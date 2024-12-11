//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "single_read.hpp"

namespace io {

class EdgeSequencesStream {
 public:
    typedef SingleReadSeq ReadT;

    explicit EdgeSequencesStream(const debruijn_graph::Graph &g)
            : graph_(g), edge_iterator_(g.e_begin<true>()) {}

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

class ComponentEdgeSequencesStream {
 public:
    typedef SingleReadSeq ReadT;

    explicit ComponentEdgeSequencesStream(const omnigraph::GraphComponent<debruijn_graph::Graph> &c)
            : c_(c), edge_iterator_(c.e_begin()) {}

    bool is_open() const {
        return true;
    }

    bool eof() {
        return edge_iterator_ == c_.e_end();
    }

    ComponentEdgeSequencesStream &operator>>(SingleReadSeq &singleread) {
        singleread = SingleReadSeq(c_.g().EdgeNucls(*edge_iterator_));
        ++edge_iterator_;
        return *this;
    }

    void close() const {}

    void reset() {
        edge_iterator_ = c_.e_begin();
    }

 private:
    const omnigraph::GraphComponent<debruijn_graph::Graph> &c_;
    decltype(c_.e_begin()) edge_iterator_;
};

}
