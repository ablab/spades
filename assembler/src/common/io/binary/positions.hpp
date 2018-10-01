//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "io/id_mapper.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"

namespace io {

namespace binary {

inline BinSaveFile &operator<<(BinSaveFile &file, const Range &range) {
    return file << range.start_pos << range.end_pos;
}

template<typename Graph>
class EdgePositionsIO : public IOSingle<typename omnigraph::EdgesPositionHandler<Graph>, EdgeMapper<Graph>> {
public:
    typedef omnigraph::EdgesPositionHandler<Graph> Type;
    typedef EdgeMapper<Graph> Mapper;
    EdgePositionsIO()
            : IOSingle<Type, Mapper>("edge positions", ".pos") {
    }

private:
    void SaveImpl(BinSaveFile &file, const Type &edge_pos) override {
        for (auto it = edge_pos.g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            auto pos_it = edge_pos.GetEdgePositions(*it);
            file << (*it).int_id() << pos_it.size();
            for (const auto &i : pos_it) {
                file << i.contigId << i.mr.initial_range << i.mr.mapped_range;
            }
        }
    }

    void LoadImpl(BinLoadFile &file, Type &edge_pos, const Mapper &mapper) override {
        edge_pos.clear();
        size_t e;
        while (file >> e) { //Read until the end
            auto eid = mapper[e];
            auto info_count = file.Read<size_t>();
            while (info_count--) {
                auto contig = file.Read<std::string>();
                size_t pos[4];
                file >> pos;
                edge_pos.AddEdgePosition(eid, contig, pos[0], pos[1], pos[2], pos[3]);
            }
        }
    };
};

template<typename Graph>
struct IOTraits<omnigraph::EdgesPositionHandler<Graph>> {
    typedef EdgePositionsIO<Graph> Type;
};

} // namespace binary

} //namespace io
