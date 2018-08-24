//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "graph_pack.hpp"

//Some backward-compatibility functions mimicking old graphio.hpp interface
namespace debruijn_graph {

namespace graphio {

template<typename Graph>
void PrintGraphPack(const string &basename, const graph_pack<Graph> &gp) {
    io::binary::GraphPackIO<Graph> io;
    io.Save(basename, gp);
}

template<typename Graph>
void PrintAll(const std::string &basename, const graph_pack<Graph> &gp) {
    io::binary::FullPackIO<Graph> io;
    io.Save(basename, gp);
}

template<typename Graph>
void ScanGraphPack(const string &basename, graph_pack<Graph> &gp) {
    io::binary::GraphPackIO<Graph> io;
    io.Load(basename, gp);
}

template<typename Graph>
void ScanAll(const std::string &basename, graph_pack<Graph> &gp) {
    io::binary::FullPackIO<Graph> io;
    io.Load(basename, gp);
}

}

}
