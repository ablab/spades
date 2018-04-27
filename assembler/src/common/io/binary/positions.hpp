//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"

namespace io {

template<typename Graph>
class EdgePositionsIO : IOBase<typename omnigraph::EdgesPositionHandler<Graph>> {
public:
    typedef omnigraph::EdgesPositionHandler<Graph> Type;
    typedef IdMapper<typename Graph::EdgeId> Mapper;
    EdgePositionsIO(const Mapper &mapper):
            IOBase<Type>("edge positions", ".pos"), mapper_(mapper) {}

private:
    void SaveImpl(SaveFile &file, const Type &edge_pos) override {
        for (auto it = edge_pos.g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            auto pos_it = edge_pos.GetEdgePositions(*it);
            file << it->int_id() << pos_it.size();
            for (const auto &i : pos_it) {
                file << i.contigId;
                file << i.mr.initial_range.start_pos << i.mr.initial_range.end_pos
                     << i.mr.mapped_range.start_pos  << i.mr.initial_range.end_pos;
            }
        }
    }

    void LoadImpl(LoadFile &file, Type &edge_pos) override {
        VERIFY(!edge_pos.IsAttached());
        edge_pos.Attach();
        while (file) { //Read until the end
            auto eid = mapper_[file.Read<size_t>()];
            auto info_count = file.Read<size_t>();
            while (info_count--) {
                auto contig = file.Read<std::string>();
                size_t start_pos, end_pos, m_start_pos, m_end_pos;
                file >> start_pos >> end_pos >> m_start_pos >> m_end_pos;
                edge_pos.AddEdgePosition(eid, contig, start_pos, end_pos, m_start_pos, m_end_pos);
            }
        }
    };

    const Mapper &mapper_;
};

}
