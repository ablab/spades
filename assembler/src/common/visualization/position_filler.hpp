//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/alignment/sequence_mapper.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "io/reads/io_helper.hpp"

namespace visualization {

namespace position_filler {

class PosFiller {
    using Graph = debruijn_graph::Graph;
    using EdgeId = Graph::EdgeId;
    using MapperPtr = std::shared_ptr<debruijn_graph::SequenceMapper<Graph>>;
    const Graph &g_;
    MapperPtr mapper_;
    omnigraph::EdgesPositionHandler<Graph> &edge_pos_;

public:
    typedef omnigraph::MappingPath<EdgeId> MappingPath;

    PosFiller(const Graph &g, MapperPtr mapper, omnigraph::EdgesPositionHandler<Graph> &edge_pos)
            : g_(g), mapper_(mapper), edge_pos_(edge_pos) {}

    PosFiller(debruijn_graph::GraphPack &gp)
            : g_(gp.get<Graph>()), mapper_(debruijn_graph::MapperInstance(gp)),
              edge_pos_(gp.get_mutable<omnigraph::EdgesPositionHandler<Graph>>()) {}

    MappingPath Process(const std::string &s, const std::string &name) const {
        return Process(io::SingleRead(name, s));
    }

    MappingPath Process(const Sequence &s, const std::string &name) const {
        return Process(s.str(), name);
    }

    MappingPath Process(const io::SingleRead &read) const {
        MappingPath path = mapper_->MapRead(read);
        const auto &name = read.name();
        int cur_pos = 0;
        TRACE("Contig " << name << " mapped on " << path.size() << " fragments.");
        for (size_t i = 0; i < path.size(); i++) {
            EdgeId ei = path[i].first;
            omnigraph::MappingRange mr = path[i].second;
            int len = (int) (mr.mapped_range.end_pos - mr.mapped_range.start_pos);
            if (i > 0 &&
                path[i - 1].first != ei &&
                g_.EdgeStart(ei) != g_.EdgeEnd(path[i - 1].first)) {
                TRACE("Contig " << name
                          << " mapped on not adjacent edge. Position in contig is "
                          << path[i - 1].second.initial_range.start_pos
                             + 1
                          << "--"
                          << path[i - 1].second.initial_range.end_pos
                          << " and "
                          << mr.initial_range.start_pos + 1
                          << "--" << mr.initial_range.end_pos);
            }
            edge_pos_.AddEdgePosition(ei, name, mr.initial_range.start_pos,
                                      mr.initial_range.end_pos,
                                      mr.mapped_range.start_pos,
                                      mr.mapped_range.end_pos);
            cur_pos += len;
        }
        return path;
    }

    void Process(io::SingleStream stream) const {
        io::SingleRead read;
        while (!stream.eof()) {
            stream >> read;
            Process(read);
        }
    }

private:
    DECL_LOGGER("PosFiller");
};

inline void FillPos(debruijn_graph::GraphPack &gp, const std::string &contig_file, const std::string &prefix, bool with_rc) {
    PosFiller pos_filler(gp);
    pos_filler.Process(io::PrefixAddingReaderWrapper(io::EasyStream(contig_file, with_rc, false), prefix));
}

inline void FillPos(debruijn_graph::GraphPack &gp, const std::string &s, const std::string &name) {
    PosFiller pos_filler(gp);
    pos_filler.Process(s, name);
}

}
}
