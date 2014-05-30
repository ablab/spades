#pragma once

//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard_base.hpp"
#include "printing_parameter_storage.hpp"

namespace omnigraph {
namespace visualization {

template<class Graph>
class VertexLinker : public virtual ParameterStorage<typename Graph::VertexId, string> {
};

template<class Graph>
class MapVertexLinker : public VertexLinker<Graph>, public MapParameterStorage<typename Graph::VertexId, string> {
public:
	MapVertexLinker() : MapParameterStorage<typename Graph::VertexId, string>("") {
	}

	MapVertexLinker(const map<typename Graph::VertexId, string> &link_map) : MapParameterStorage<typename Graph::VertexId, string>(link_map, "") {
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
