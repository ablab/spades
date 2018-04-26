//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/id_mapper.hpp"
#include "io_base.hpp"

namespace debruijn_graph {

namespace graphio {

using io::IdMapper;

template<typename Graph>
class GraphIO : public IOBase<Graph> {
public:
    GraphIO():
            IOBase<Graph>("debruijn graph", ".grp") {}

    const IdMapper<typename Graph::EdgeId> &GetEdgeMapper() {
        return edge_mapper_;
    }

protected:
    void SaveImpl(SaveFile &file, const Graph &graph) override {
        file << graph.GetGraphIdDistributor().GetMax();

        auto component = GraphComponent<Graph>::WholeGraph(graph);
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
        restricted::IdSegmentStorage id_storage = graph.GetGraphIdDistributor().ReserveUpTo(max_id);

        size_t vertex_count, edge_count;
        file >> vertex_count >> edge_count;

        while (vertex_count--) {
            size_t ids[2];
            file >> ids[0] >> ids[1];
            TRACE("Vertex " << ids[0] << " ~ " << ids[1] << " .");

            if (!vertex_mapper_.count(ids[0])) {
                auto id_distributor = id_storage.GetSegmentIdDistributor(ids, ids + 2);
                auto new_id = graph.AddVertex(typename Graph::VertexData(), id_distributor);
                vertex_mapper_[ids[0]] = new_id;
                vertex_mapper_[ids[1]] = graph.conjugate(new_id);
            }
        }

        while (edge_count--) {
            size_t ids[2];
            size_t start_id, fin_id;
            Sequence seq;
            file >> ids[0] >> ids[1] >> start_id >> fin_id >> seq;
            TRACE("Edge " << ids[0] << " : " << start_id << " -> "
                          << fin_id << " l = " << seq.size() << " ~ " << ids[1]);
            if (!edge_mapper_.count(ids[0])) {
                auto id_distributor = id_storage.GetSegmentIdDistributor(ids, ids + 2);
                auto new_id = graph.AddEdge(vertex_mapper_[start_id], vertex_mapper_[fin_id], seq, id_distributor);
                edge_mapper_[ids[0]] = new_id;
                edge_mapper_[ids[1]] = graph.conjugate(new_id);
            }
        }
    }
private:
    IdMapper<typename Graph::VertexId> vertex_mapper_;
    IdMapper<typename Graph::EdgeId> edge_mapper_;
};

}

}
