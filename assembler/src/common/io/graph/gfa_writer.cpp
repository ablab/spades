//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gfa_writer.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/graph_iterators.hpp"
#include "assembly_graph/components/graph_component.hpp"

using namespace gfa;
using namespace debruijn_graph;

template class omnigraph::GraphComponent<Graph>;

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


void GFAWriter::WriteSegments(const Component &gc) {
    for (EdgeId e : gc.edges()) {
        if (e <= graph_.conjugate(e)) {
            WriteSegment(edge_namer_.EdgeString(e), graph_.EdgeNucls(e),
                         graph_.coverage(e) * double(graph_.length(e)),
                         os_);
        }
    }
}

void GFAWriter::WriteLinks(const Component &gc) {
    for (VertexId v : gc.vertices()) {
        if (v <= graph_.conjugate(v) && !gc.IsBorder(v)) {
            for (auto inc_edge : graph_.IncomingEdges(v)) {
                for (auto out_edge : graph_.OutgoingEdges(v)) {
                    WriteLink(inc_edge, out_edge, graph_.k(),
                              os_, edge_namer_);

                }
            }
        }
    }
}

                    void GFAComponentWriter::WriteSegments() {
    for (auto e : component_.edges()) {
        if (e.int_id() > component_.g().conjugate(e).int_id())
            continue;
        WriteSegment(edge_namer_.EdgeString(e), component_.g().EdgeNucls(e),
                     component_.g().coverage(e) * double(component_.g().length(e)),
                     os_);
    }
}

void GFAComponentWriter::WriteLinks() {
    //TODO switch to constant vertex iterator
    for (auto v : component_.vertices()) {
        if (v.int_id() > component_.g().conjugate(v).int_id())
            continue;
        for (auto inc_edge : component_.g().IncomingEdges(v)) {
            if (component_.contains(inc_edge)) {
                for (auto out_edge : component_.g().OutgoingEdges(v)) {
                    if (component_.contains(out_edge)) {
                        WriteLink(inc_edge, out_edge, component_.g().k(),
                                  os_, edge_namer_);
                    }
                }
            }
        }
    }
}

void GFAWriter::WriteSegmentsAndLinks(const Component &gc) {
    //TODO remove and add optional check?
    auto rc_closure = Component::FromEdges(graph_, gc.edges().begin(), gc.edges().end(),
                                           /*add_conjugate*/true);
    WriteSegments(rc_closure);
    WriteLinks(rc_closure);
}

