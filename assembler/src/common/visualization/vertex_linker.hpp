#pragma once

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/standard_base.hpp"
#include "printing_parameter_storage.hpp"

namespace visualization {

namespace vertex_linker {

template<class Graph>
class VertexLinker
        : public virtual printing_parameter_storage::ParameterStorage<typename Graph::VertexId, string> {
};

template<class Graph>
class MapVertexLinker : public VertexLinker<Graph>,
                        public printing_parameter_storage::MapParameterStorage<typename Graph::VertexId, string> {
public:
    MapVertexLinker() : printing_parameter_storage::MapParameterStorage<typename Graph::VertexId, string>("") {
    }

    MapVertexLinker(const map<typename Graph::VertexId, string> &link_map) :
            printing_parameter_storage::MapParameterStorage<typename Graph::VertexId, string>(link_map, "") {
    }

    virtual ~MapVertexLinker() {
    }
};

template<class Graph>
class EmptyGraphLinker : public MapVertexLinker<Graph> {
public:
    EmptyGraphLinker() {
    }
};

}

}