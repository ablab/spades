//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/id_mapper.hpp"
#include "io_base.hpp"

#include "common/assembly_graph/components/graph_component.hpp"
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
    void SaveImpl(BinSaveFile &file, const Graph &graph) override {
        file << graph.GetGraphIdDistributor().GetMax();

        for (auto v1 : graph) {
            file << v1.int_id() << graph.conjugate(v1).int_id();
            std::unordered_set<EdgeId> to_write; //TODO: reserve bytes, then rewrite in the end?
            for (auto e : graph.OutgoingEdges(v1)) {
                if (e <= graph.conjugate(e))
                    to_write.insert(e);
            }
            file << (size_t)to_write.size();
            for (auto e1 : to_write) {
                auto e2 = graph.conjugate(e1);
                file << e1.int_id() << e2.int_id()
                     << graph.EdgeEnd(e1).int_id() << graph.EdgeStart(e2).int_id()
                     << graph.EdgeNucls(e1);
            }
        }
    }

    void LoadImpl(BinLoadFile &file, Graph &graph) override {
        size_t max_id;
        file >> max_id;
        auto id_storage = graph.GetGraphIdDistributor().Reserve(max_id, /*force_zero_shift*/true);

        auto TryAddVertex = [&](size_t ids[2]) {
            if (vertex_mapper_.count(ids[0]))
                return;
            TRACE("Vertex " << ids[0] << " ~ " << ids[1] << " .");
            auto id_distributor = id_storage.GetSegmentIdDistributor(ids, ids + 2);
            auto new_id = graph.AddVertex(typename Graph::VertexData(), id_distributor);
            vertex_mapper_[ids[0]] = new_id;
            vertex_mapper_[ids[1]] = graph.conjugate(new_id);
        };

        size_t start_ids[2];
        while (file >> start_ids) { //Read until the end
            TryAddVertex(start_ids);
            auto count = file.Read<size_t>();
            while (count--) {
                size_t edge_ids[2], end_ids[2];
                Sequence seq;
                file >> edge_ids >> end_ids >> seq;
                TRACE("Edge " << edge_ids[0] << " : " << start_ids[0] << " -> "
                              << end_ids[0] << " l = " << seq.size() << " ~ " << edge_ids[1]);
                TryAddVertex(end_ids);

                VERIFY_MSG(!edge_mapper_.count(edge_ids[0]), edge_ids[0] << " is not unique");
                auto id_distributor = id_storage.GetSegmentIdDistributor(edge_ids, edge_ids + 2);
                auto new_id = graph.AddEdge(vertex_mapper_[start_ids[0]], vertex_mapper_[end_ids[0]], seq, id_distributor);
                edge_mapper_[edge_ids[0]] = new_id;
                edge_mapper_[edge_ids[1]] = graph.conjugate(new_id);
            }
        }
    }

    IdMapper<typename Graph::VertexId> vertex_mapper_;
    IdMapper<typename Graph::EdgeId> edge_mapper_;

    DECL_LOGGER("GraphIO");
};

} // namespace binary

} // namespace io
