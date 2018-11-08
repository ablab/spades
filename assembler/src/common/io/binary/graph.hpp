//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/id_mapper.hpp"
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

    const IdMapper<typename Graph::EdgeId> &GetEdgeMapper() {
        return edge_mapper_;
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
            str << (size_t)0; //null-term
        }
    }

    void Read(BinIStream &str, Graph &graph) override {
        graph.clear();

        uint64_t max_vid, max_eid;
        str >> max_vid >> max_eid;
        graph.reserve(max_vid, max_eid);

        auto TryAddVertex = [&](uint64_t ids[2]) {
            //FIXME: use fast check
            if (vertex_mapper_.count(ids[0]))
                return;
            TRACE("Vertex " << ids[0] << " ~ " << ids[1] << " .");
            auto new_id = graph.AddVertex(typename Graph::VertexData(), std::min(ids[0], ids[1]));
            VERIFY(new_id == std::min(ids[0], ids[1]));
            VERIFY(graph.conjugate(new_id) == std::max(ids[0], ids[1]));
            vertex_mapper_[ids[0]] = ids[0];
            vertex_mapper_[ids[1]] = ids[1];
        };

        while (str.has_data()) { //Read until the end
            // FIXME use two separate ids instead of C-array! C-array are error-prone and could be easily mixed up with pointers!
            size_t start_ids[2];
            str >> start_ids;
            TryAddVertex(start_ids);
            while (true) {
                size_t edge_ids[2];
                str >> edge_ids[0];
                if (!edge_ids[0]) //null-term
                    break;
                str >> edge_ids[1];
                size_t end_ids[2];
                Sequence seq;
                str >> end_ids >> seq;
                TRACE("Edge " << edge_ids[0] << " : " << start_ids[0] << " -> "
                              << end_ids[0] << " l = " << seq.size() << " ~ " << edge_ids[1]);
                TryAddVertex(end_ids);

                auto new_id = graph.AddEdge(start_ids[0], end_ids[0], seq, edge_ids[0], edge_ids[1]);
                VERIFY(new_id == edge_ids[0]);
                VERIFY(graph.conjugate(new_id) == edge_ids[1]);
            }
        }
    }

private:
    IdMapper<typename Graph::VertexId> vertex_mapper_;
    IdMapper<typename Graph::EdgeId> edge_mapper_;

    DECL_LOGGER("GraphIO");
};

template<>
struct IOTraits<debruijn_graph::Graph> {
    typedef GraphIO<debruijn_graph::Graph> Type;
};

} // namespace binary

} // namespace io
