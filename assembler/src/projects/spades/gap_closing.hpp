#pragma once

#include "assembly_graph/graph_core/graph.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"
#include "assembly_graph/graph_support/edge_removal.hpp"
#include "algorithms/simplification/compressor.hpp"

namespace debruijn_graph {

namespace gap_closing {
typedef omnigraph::GapDescription<Graph> GapDescription;

class GapJoiner {
    Graph& g_;
    omnigraph::EdgeRemover<Graph> edge_remover_;

    EdgeId ClipEnd(EdgeId e, size_t pos) {
        VERIFY(pos > 0);
        VERIFY(omnigraph::TerminalVertexCondition<Graph>(g_).Check(g_.EdgeEnd(e)));
        VERIFY(e != g_.conjugate(e));
        if (pos == g_.length(e)) {
            return e;
        } else {
            auto split_res = g_.SplitEdge(e, pos);
            edge_remover_.DeleteEdge(split_res.second);
            return split_res.first;
        }
    }

    EdgeId ClipStart(EdgeId e, size_t pos) {
        return g_.conjugate(ClipEnd(g_.conjugate(e), g_.length(e) - pos));
    }

public:
    GapJoiner(Graph& g) :
            g_(g),
            edge_remover_(g) {
    }

    EdgeId operator() (const GapDescription& gap, bool compress = true) {
        VERIFY(gap.start != gap.end && gap.start != g_.conjugate(gap.end));
        EdgeId start = ClipEnd(gap.start, gap.edge_gap_start_position);
        EdgeId end = ClipStart(gap.end, gap.edge_gap_end_position);
        VERIFY(g_.EdgeNucls(start).end<runtime_k::RtSeq>(g_.k()) == gap.gap_seq.start<runtime_k::RtSeq>(g_.k())
               && g_.EdgeNucls(end).start<runtime_k::RtSeq>(g_.k()) == gap.gap_seq.end<runtime_k::RtSeq>(g_.k()));
        EdgeId new_edge = g_.AddEdge(g_.EdgeEnd(start), g_.EdgeStart(end), gap.gap_seq);

        if (compress) {
            return omnigraph::Compressor<Graph>(g_).CompressVertexEdgeId(g_.EdgeStart(new_edge));
        } else {
            return new_edge;
        }
    }
};

}
}