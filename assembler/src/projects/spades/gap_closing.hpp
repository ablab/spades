#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"
#include "assembly_graph/graph_support/edge_removal.hpp"
#include "modules/simplification/compressor.hpp"
#include "modules/alignment/gap_info.hpp"

namespace debruijn_graph {

namespace gap_closing {
typedef omnigraph::GapDescription<Graph> GapDescription;

class GapJoiner {
    Graph& g_;
    omnigraph::EdgeRemover<Graph> edge_remover_;

    EdgeId ClipEnd(EdgeId e, size_t to_trim) {
        VERIFY_MSG(to_trim < g_.length(e), "Asked to trim " << to_trim << " edge " << g_.str(e));
        VERIFY(omnigraph::TerminalVertexCondition<Graph>(g_).Check(g_.EdgeEnd(e)));
        VERIFY(e != g_.conjugate(e));
        if (to_trim == 0) {
            return e;
        } else {
            auto split_res = g_.SplitEdge(e, g_.length(e) - to_trim);
            edge_remover_.DeleteEdge(split_res.second);
            return split_res.first;
        }
    }

    EdgeId ClipStart(EdgeId e, size_t to_trim) {
        return g_.conjugate(ClipEnd(g_.conjugate(e), to_trim));
    }

    EdgeId AddEdge(VertexId v1, VertexId v2, const Sequence &gap_seq) {
        DEBUG("Adding gap seq " << gap_seq);
        DEBUG("Between vertices " << g_.VertexNucls(v1) << " and " << g_.VertexNucls(v2));
        return g_.AddEdge(v1, v2, g_.VertexNucls(v1) + gap_seq + g_.VertexNucls(v2));
    }

public:
    GapJoiner(Graph& g) :
            g_(g),
            edge_remover_(g) {
    }

    EdgeId operator() (const GapDescription& gap, bool compress = true) {
        VERIFY(gap.left() != gap.right() && gap.left() != g_.conjugate(gap.right()));
        DEBUG("Processing gap " << gap.str(g_));
        EdgeId start = ClipEnd(gap.left(), gap.left_trim());
        EdgeId end = ClipStart(gap.right(), gap.right_trim());
        EdgeId new_edge = AddEdge(g_.EdgeEnd(start), g_.EdgeStart(end), gap.filling_seq());

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
