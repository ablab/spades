#pragma once

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "printing_parameter_storage.hpp"

namespace visualization {

namespace vertex_linker {

template<class Graph>
class VertexLinker
        : public virtual printing_parameter_storage::ParameterStorage<typename Graph::VertexId, std::string> {
};

template<class Graph>
class MapVertexLinker : public VertexLinker<Graph>,
                        public printing_parameter_storage::MapParameterStorage<typename Graph::VertexId, std::string> {
public:
    MapVertexLinker() : printing_parameter_storage::MapParameterStorage<typename Graph::VertexId,std:: string>("") {
    }

    MapVertexLinker(const std::map<typename Graph::VertexId, std::string> &link_map) :
            printing_parameter_storage::MapParameterStorage<typename Graph::VertexId, std::string>(link_map, "") {
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