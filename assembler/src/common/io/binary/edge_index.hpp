//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "modules/alignment/edge_index.hpp"

namespace io {

namespace binary {

template<typename Graph>
class EdgeIndexIO : public IOSingle<debruijn_graph::EdgeIndex<Graph>> {
public:
    typedef debruijn_graph::EdgeIndex<Graph> Type;
    EdgeIndexIO()
            : IOSingle<Type>("edge index", ".kmidx") {
    }

private:
    void SaveImpl(BinSaveFile &file, const Type &value) override {
        const auto &index = value.inner_index();
        file << (uint32_t)index.k() << index;
    }

    void LoadImpl(BinLoadFile &file, Type &value) override {
        auto &index = value.inner_index();
        uint32_t k_;
        file >> k_;
        VERIFY_MSG(k_ == index.k(), "Cannot read " << this->name_ << ", different Ks");
        index.clear();
        file >> index;
        value.Update();
    }
};

} // namespace binary

} // namespace io
