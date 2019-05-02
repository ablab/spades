//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"

#include "assembly_graph/core/graph.hpp"
#include "common/sequence/sequence.hpp"

namespace io {

namespace binary {

template<typename Graph>
class GraphIO : public IOSingle<Graph> {
public:
    GraphIO()
            : IOSingle<Graph>("debruijn graph", ".grseq") {
    }

private:
    void Write(BinOStream &str, const Graph &graph) override {
        str << graph.vreserved() << graph.ereserved();

        for (auto v1 : graph) {
            str << v1.int_id() << graph.conjugate(v1).int_id();
            for (auto e1 : graph.OutgoingEdges(v1)) {
                auto e2 = graph.conjugate(e1);
                if (e2 < e1)
                    continue;
                str << e1.int_id() << e2.int_id()
                    << graph.EdgeEnd(e1).int_id() << graph.EdgeStart(e2).int_id()
                    << graph.EdgeNucls(e1);
            }
            str << (uint64_t)0; //null-term
        }
    }

    void Read(BinIStream &str, Graph &graph) override {
        graph.clear();

        uint64_t max_vid, max_eid;
        str >> max_vid >> max_eid;
        graph.reserve(max_vid, max_eid);

        auto TryAddVertex = [&](uint64_t ids[2]) {
            if (graph.contains(typename Graph::VertexId(ids[0])))
                return;
            TRACE("Vertex " << ids[0] << " ~ " << ids[1] << " .");
            auto new_id = graph.AddVertex(typename Graph::VertexData(), ids[0], ids[1]);
            VERIFY(new_id == ids[0]);
            VERIFY(graph.conjugate(new_id) == ids[1]);
        };

        while (str.has_data()) { //Read until the end
            // FIXME use two separate ids instead of C-array! C-array are error-prone and could be easily mixed up with pointers!
            uint64_t start_ids[2];
            str >> start_ids;
            TryAddVertex(start_ids);
            while (true) {
                uint64_t edge_ids[2];
                str >> edge_ids[0];
                if (!edge_ids[0]) //null-term
                    break;
                str >> edge_ids[1];
                uint64_t end_ids[2];
                Sequence seq;
                str >> end_ids >> seq;
                TRACE("Edge " << edge_ids[0] << " : " << start_ids[0] << " -> "
                              << end_ids[0] << " l = " << seq.size() << " ~ " << edge_ids[1]);
                TryAddVertex(end_ids);

                auto new_id = graph.AddEdge(start_ids[0], end_ids[0],
                        typename Graph::EdgeData(seq), edge_ids[0], edge_ids[1]);
                VERIFY(new_id == edge_ids[0]);
                VERIFY(graph.conjugate(new_id) == edge_ids[1]);
            }
        }
    }

    DECL_LOGGER("GraphIO");
};

template<>
struct IOTraits<debruijn_graph::Graph> {
    typedef GraphIO<debruijn_graph::Graph> Type;
};

} // namespace binary

} // namespace io
