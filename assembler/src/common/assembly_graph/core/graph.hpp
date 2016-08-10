//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "observable_graph.hpp"
#include "coverage.hpp"
#include "debruijn_data.hpp"

namespace debruijn_graph {

using omnigraph::CoverageIndex;
class DeBruijnGraph: public omnigraph::ObservableGraph<DeBruijnDataMaster> {
public:
    typedef omnigraph::ObservableGraph<DeBruijnDataMaster> base;
    typedef base::DataMasterT DataMasterT;
    typedef base::VertexData VertexData;
    typedef base::EdgeData EdgeData;
    typedef base::EdgeId EdgeId;
    typedef base::VertexId VertexId;
    typedef base::VertexIt VertexIt;
    typedef VertexIt VertexIterator;
    typedef VertexIterator iterator; // for for_each
    typedef const VertexIterator const_iterator; // for for_each
private:
    CoverageIndex<DeBruijnGraph> coverage_index_;

public:
    DeBruijnGraph(size_t k) :
            base(k), coverage_index_(*this) {
    }

    CoverageIndex<DeBruijnGraph>& coverage_index() {
        return coverage_index_;
    }

    const CoverageIndex<DeBruijnGraph>& coverage_index() const {
        return coverage_index_;
    }

    /**
     * Method returns average coverage of the edge
     */
    double coverage(EdgeId edge) const {
        return coverage_index_.coverage(edge);
    }

    using base::AddVertex;
    using base::AddEdge;

    VertexId AddVertex() {
        return AddVertex(VertexData());
    }

    EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
        VERIFY(nucls.size() > k());
        return AddEdge(from, to, EdgeData(nucls));
    }

    size_t k() const {
        return master().k();
    }

    /**
     * Method returns Sequence stored in the edge
     */
    const Sequence& EdgeNucls(EdgeId edge) const {
        return this->data(edge).nucls();
    }

    const Sequence VertexNucls(VertexId v) const {
        //todo add verify on vertex nucls consistency
        if (this->OutgoingEdgeCount(v) > 0) {
            return EdgeNucls(*(this->out_begin(v))).Subseq(0, k());
        } else if (this->IncomingEdgeCount(v) > 0) {
            EdgeId inc = *(this->in_begin(v));
            size_t length = EdgeNucls(inc).size();
            return EdgeNucls(inc).Subseq(length - k(), length);
        }
        VERIFY(false);
        return Sequence();
    }

    Sequence PathNucls(const vector<EdgeId> &path) const {
        if(path.empty())
            return Sequence("");
        SequenceBuilder result;
        result.append(Sequence(""));
        result.append(this->EdgeNucls(path[0]).Subseq(0, this->k()));
        for (size_t i = 0; i < path.size(); ++i) {
            result.append(this->EdgeNucls(path[i]).Subseq(this->k()));
        }

        return result.BuildSequence();
    }

private:
    DECL_LOGGER("DeBruijnGraph")
};

typedef DeBruijnGraph ConjugateDeBruijnGraph;

typedef ConjugateDeBruijnGraph Graph;
typedef Graph::EdgeId EdgeId;
typedef Graph::VertexId VertexId;
}
