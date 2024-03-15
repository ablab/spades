//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"

namespace io {

namespace binary {

inline BinOStream &operator<<(BinOStream &str, const Range &range) {
    return str << range.start_pos << range.end_pos;
}

template<typename Graph>
class EdgePositionsIO : public IOSingle<typename omnigraph::EdgesPositionHandler<Graph>> {
public:
    typedef omnigraph::EdgesPositionHandler<Graph> Type;
    EdgePositionsIO()
            : IOSingle<Type>("edge positions", ".pos") {
    }

    void SaveImpl(BinOStream &str, const Type &edge_pos) override {
        size_t edges_cnt = edge_pos.g().e_size();
        str << edges_cnt;

        for (auto it = edge_pos.g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            auto pos_it = edge_pos.GetEdgePositions(*it);
            str << (*it).int_id() << pos_it.size();
            for (const auto &i : pos_it) {
                str << i.contigId << i.mr.initial_range << i.mr.mapped_range;
            }
        }
    }

    void LoadImpl(BinIStream &str, Type &edge_pos) override {
        edge_pos.clear();
        size_t edge_cnt;
        str >> edge_cnt;

        for (size_t i = 0; i < edge_cnt; ++i) {
            uint64_t eid;
            str >> eid;
            size_t info_count;
            str >> info_count;
            while (info_count--) {
                auto contig = str.Read<std::string>();
                size_t pos[4];
                str >> pos;
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
