//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "graph.hpp"
#include "coverage.hpp"

namespace io {

namespace binary {

/**
 * @brief  This IOer processes the graph with its coverage index.
 */
template<typename Graph>
class BasicGraphIO : public GraphIO<Graph> {
    typedef GraphIO<Graph> Base;

public:
    void Save(const std::string &basename, const Graph &graph) override {
        Base::Save(basename, graph);
        io::binary::Save(basename, graph.coverage_index());
    }

    bool Load(const std::string &basename, Graph &graph) override {
        bool loaded = Base::Load(basename, graph);
        VERIFY(loaded);
        loaded = io::binary::Load(basename, graph.coverage_index());
        VERIFY(loaded);
        return true;
    }

    void BinWrite(std::ostream &os, const Graph &graph) override {
        Base::BinWrite(os, graph);
        io::binary::Write(os, graph.coverage_index());
    }

    bool BinRead(std::istream &is, Graph &graph) override {
        bool loaded = Base::BinRead(is, graph);
        VERIFY(loaded);
        loaded = io::binary::Read(is, graph.coverage_index());
        VERIFY(loaded);
        return true;
    }
};

} // namespace binary

} // namespace io
