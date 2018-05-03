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

template<typename Graph>
class GraphIO : public IOSingle<Graph> {
public:
    GraphIO()
            : IOSingle<Graph>("debruijn graph", ".grp") {
    }

    const IdMapper<typename Graph::EdgeId> &GetEdgeMapper() {
        return edge_mapper_;
    }

private:
    void SaveImpl(SaveFile &file, const Graph &graph) override {
        file << graph.GetGraphIdDistributor().GetMax();

        auto component = omnigraph::GraphComponent<Graph>::WholeGraph(graph);
        file << component.v_size() << component.e_size();

        for (auto iter = component.v_begin(); iter != component.v_end(); ++iter) {
            auto v = *iter;
            file << v.int_id() << graph.conjugate(v).int_id();
        }

        for (auto iter = component.e_begin(); iter != component.e_end(); ++iter) {
            auto e = *iter;
            file << e.int_id() << graph.conjugate(e).int_id()
                 << graph.EdgeStart(e).int_id() << graph.EdgeEnd(e).int_id()
                 << graph.EdgeNucls(e);
        }
    }

    void LoadImpl(LoadFile &file, Graph &graph) override {
        size_t max_id;
        file >> max_id;
        auto id_storage = graph.GetGraphIdDistributor().ReserveUpTo(max_id);

        size_t vertex_count, edge_count;
        file >> vertex_count >> edge_count;

        while (vertex_count--) {
            size_t ids[2];
            file >> ids;
            TRACE("Vertex " << ids[0] << " ~ " << ids[1] << " .");

            if (!vertex_mapper_.count(ids[0])) {
                auto id_distributor = id_storage.GetSegmentIdDistributor(ids, ids + 2);
                auto new_id = graph.AddVertex(typename Graph::VertexData(), id_distributor);
                vertex_mapper_[ids[0]] = new_id;
                vertex_mapper_[ids[1]] = graph.conjugate(new_id);
            }
        }

        while (edge_count--) {
            size_t ids[2], vertex_ids[2];
            Sequence seq;
            file >> ids >> vertex_ids >> seq;
            TRACE("Edge " << ids[0] << " : " << vertex_ids[0] << " -> "
                          << vertex_ids[1] << " l = " << seq.size() << " ~ " << ids[1]);
            if (!edge_mapper_.count(ids[0])) {
                auto id_distributor = id_storage.GetSegmentIdDistributor(ids, ids + 2);
                auto new_id = graph.AddEdge(vertex_mapper_[vertex_ids[0]], vertex_mapper_[vertex_ids[1]], seq, id_distributor);
                edge_mapper_[ids[0]] = new_id;
                edge_mapper_[ids[1]] = graph.conjugate(new_id);
            }
        }
    }

    IdMapper<typename Graph::VertexId> vertex_mapper_;
    IdMapper<typename Graph::EdgeId> edge_mapper_;
};

}
