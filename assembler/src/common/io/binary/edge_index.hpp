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

    void Write(BinOStream &str, const Type &value) override {
        const auto &index = value.inner_index();
        str << (uint32_t)index.k() << index;
    }

    void Read(BinIStream &str, Type &value) override {
        auto &index = value.inner_index();
        uint32_t k_;
        str >> k_;
        VERIFY_MSG(k_ == index.k(), "Cannot read edge index, different Ks");
        index.clear();
        str >> index;
        value.Update();
    }
};

template<typename Graph>
struct IOTraits<debruijn_graph::EdgeIndex<Graph>> {
    typedef EdgeIndexIO<Graph> Type;
};

} // namespace binary

} // namespace io
