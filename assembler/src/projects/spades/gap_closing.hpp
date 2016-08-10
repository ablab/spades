#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"
#include "assembly_graph/graph_support/edge_removal.hpp"
#include "modules/simplification/compressor.hpp"

namespace debruijn_graph {

namespace gap_closing {
typedef omnigraph::GapDescription<Graph> GapDescription;

class GapJoiner {
    Graph& g_;
    omnigraph::EdgeRemover<Graph> edge_remover_;
    bool add_flanks_;

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

    EdgeId AddEdge(VertexId v1, VertexId v2, const Sequence& gap_seq) {
        if (!add_flanks_) {
            VERIFY_MSG(g_.VertexNucls(v1) == gap_seq.Subseq(0, g_.k()), 
                       g_.VertexNucls(v1) << " not equal " << gap_seq.Subseq(0, g_.k()));
            VERIFY_MSG(g_.VertexNucls(v2) == gap_seq.Subseq(gap_seq.size() - g_.k()),
                       g_.VertexNucls(v2) << " not equal " << gap_seq.Subseq(gap_seq.size() - g_.k()));
            return g_.AddEdge(v1, v2, gap_seq);
        } else {
            DEBUG("Adding gap seq " << gap_seq);
            DEBUG("Between vertices " << g_.VertexNucls(v1) << " and " << g_.VertexNucls(v2));
            return g_.AddEdge(v1, v2, g_.VertexNucls(v1) + gap_seq + g_.VertexNucls(v2));
        }
    }

public:
    GapJoiner(Graph& g, bool add_flanks = false) :
            g_(g),
            edge_remover_(g),
            add_flanks_(add_flanks) {
    }

    EdgeId operator() (const GapDescription& gap, bool compress = true) {
        VERIFY(gap.start != gap.end && gap.start != g_.conjugate(gap.end));
        DEBUG("Processing gap " << gap.str(g_));
        EdgeId start = ClipEnd(gap.start, gap.edge_gap_start_position);
        EdgeId end = ClipStart(gap.end, gap.edge_gap_end_position);
        EdgeId new_edge = AddEdge(g_.EdgeEnd(start), g_.EdgeStart(end), gap.gap_seq);

        if (compress) {
            return omnigraph::Compressor<Graph>(g_).CompressVertexEdgeId(g_.EdgeStart(new_edge));
        } else {
            return new_edge;
        }
    }
private:
    DECL_LOGGER("GapJoiner");
};

}
}
