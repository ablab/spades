//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "assembly_graph/index/edge_position_index.hpp"

namespace io {

template<typename Graph>
class EdgeIndexIO : public IOSingle<debruijn_graph::KmerFreeEdgeIndex<Graph>> {
public:
    typedef debruijn_graph::KmerFreeEdgeIndex<Graph> Type;
    EdgeIndexIO()
            : IOSingle<Type>("edge index", ".kmidx") {
    }

private:
    void SaveImpl(SaveFile &file, const Type &index) override {
        file << (uint32_t)index.k() << index;
    }

    void LoadImpl(LoadFile &file, Type &index) override {
        uint32_t k_;
        file >> k_;
        VERIFY_MSG(k_ == index.k(), "Cannot read " << this->name_ << ", different Ks");
        index.clear();
        file >> index;
    }
};

}
