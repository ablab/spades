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
};

} // namespace binary

} // namespace io
