//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gfa_writer.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/graph_iterators.hpp"

using namespace gfa;
using namespace debruijn_graph;

static void WriteSegment(const std::string& edge_id, const Sequence &seq, double cov,
                         std::ostream &os) {
    os << "S\t"
       << edge_id << '\t' << seq.str() << '\t'
       << "KC:i:" << size_t(math::round(cov)) << '\n';
}

static void WriteLink(EdgeId e1, EdgeId e2, size_t overlap_size,
                      std::ostream &os, io::CanonicalEdgeHelper<Graph> &namer) {
    os << "L\t"
       << namer.EdgeOrientationString(e1, "\t") << '\t'
       << namer.EdgeOrientationString(e2, "\t") << '\t'
       << overlap_size << "M\n";
}

void GFAWriter::WriteSegments() {
    for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        WriteSegment(edge_namer_.EdgeString(e), graph_.EdgeNucls(e),
                     graph_.coverage(e) * double(graph_.length(e)),
                     os_);
    }
}

void GFAWriter::WriteLinks() {
    //TODO switch to constant vertex iterator
    for (auto it = graph_.SmartVertexBegin(/*canonical only*/true); !it.IsEnd(); ++it) {
        VertexId v = *it;
        for (auto inc_edge : graph_.IncomingEdges(v)) {
            for (auto out_edge : graph_.OutgoingEdges(v)) {
                WriteLink(inc_edge, out_edge, graph_.k(),
                          os_, edge_namer_);
            }
        }
    }
}
