//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "modules/alignment/kmer_mapper.hpp"

namespace io {

template<typename Graph>
class KmerMapperIO : IOSingle<debruijn_graph::KmerMapper<Graph>> {
public:
    typedef debruijn_graph::KmerMapper<Graph> Type;
    KmerMapperIO()
            : IOSingle<Type>("kmer mapper", ".kmm") {
    }

private:
    void SaveImpl(SaveFile &file, const Type &mapper) override {
        file << (uint32_t)mapper.k() << mapper;
    }

    void LoadImpl(LoadFile &file, Type &mapper) override {
        uint32_t k_;
        file >> k_;
        VERIFY_MSG(k_ == mapper.k(), "Cannot read " << this->name_ << ", different Ks");
        mapper.clear();
        file >> mapper;
    }
};

}
