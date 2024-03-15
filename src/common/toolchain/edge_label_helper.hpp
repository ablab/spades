//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/handlers/id_track_handler.hpp"
#include "io/utils/id_mapper.hpp"
#include "io/utils/edge_namer.hpp"

namespace io {

template<class Graph>
class EdgeLabelHelper {
    typedef typename Graph::EdgeId EdgeId;

    const omnigraph::GraphElementFinder<Graph> &element_finder_;
    std::unique_ptr<io::IdMapper<std::string>> id_mapper_;
    EdgeNamingF<Graph> edge_naming_f_;

    mutable bool label_map_filled_;
    mutable std::unordered_map<std::string, size_t> label2id_;

public:

    explicit EdgeLabelHelper(const omnigraph::GraphElementFinder<Graph> &element_finder,
                             io::IdMapper<std::string> *id_mapper) :
            element_finder_(element_finder),
            id_mapper_(id_mapper),
            edge_naming_f_(id_mapper ? io::MapNamingF<Graph>(*id_mapper_) : io::IdNamingF<Graph>()),
            label_map_filled_(false) {
    }

//    void reset_id_mapper(const io::IdMapper<std::string> *id_mapper) {
//        if (id_mapper) {
//            id_mapper_.reset(id_mapper);
//            edge_naming_f_ = io::MapNamingF<Graph>(*id_mapper_);
//        }
//    }

    const EdgeNamingF<Graph> &edge_naming_f() const {
        return edge_naming_f_;
    }

    std::string label(EdgeId e) const {
        return edge_naming_f_(element_finder_.g(), e);
    }

    EdgeId edge(const std::string &label) const {
        if (id_mapper_ && !label_map_filled_) {
            DEBUG("Creating label -> int_id mapping");
            for (auto it = element_finder_.g().ConstEdgeBegin(/*canonical only*/true); !it.IsEnd(); ++it) {
                size_t id = element_finder_.g().int_id(*it);
                DEBUG("Mapping " << (*id_mapper_)[id] << " to " << id);
                label2id_[(*id_mapper_)[id]] = id;
            }
            label_map_filled_ = true;
        }
        return element_finder_.ReturnEdgeId(id_mapper_ ?
                                            utils::get(label2id_, label) :
                                            std::stoi(label));
    }

private:
    DECL_LOGGER("EdgeLabelHelper");
};
}
