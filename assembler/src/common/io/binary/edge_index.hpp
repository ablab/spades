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

    void SaveImpl(BinOStream &str, const Type &value) override {
        str << (uint32_t)value.k() << value;
    }

    void LoadImpl(BinIStream &str, Type &value) override {
        uint32_t k_;
        str >> k_;
        CHECK_FATAL_ERROR(k_ == value.k(), "Cannot read edge index, different Ks");
        value.clear();
        str >> value;
    }
};

template<typename Graph>
struct IOTraits<debruijn_graph::EdgeIndex<Graph>> {
    typedef EdgeIndexIO<Graph> Type;
};

} // namespace binary

} // namespace io
