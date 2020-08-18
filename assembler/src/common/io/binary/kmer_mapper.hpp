//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "modules/alignment/kmer_mapper.hpp"

namespace io {

namespace binary {

template<typename Graph>
class KmerMapperIO : public IOSingle<debruijn_graph::KmerMapper<Graph>> {
public:
    typedef debruijn_graph::KmerMapper<Graph> Type;
    KmerMapperIO()
            : IOSingle<Type>("kmer mapper", ".kmm") {
    }

    void SaveImpl(BinOStream &str, const Type &mapper) override {
        str << (uint32_t)mapper.k() << mapper;
    }

    void LoadImpl(BinIStream &str, Type &mapper) override {
        uint32_t k_;
        str >> k_;
        CHECK_FATAL_ERROR(k_ == mapper.k(), "Cannot read kmer mapper, different Ks");
        mapper.clear();
        str >> mapper;
    }
};

template<typename Graph>
struct IOTraits<debruijn_graph::KmerMapper<Graph>> {
    typedef KmerMapperIO<Graph> Type;
};

} // namespace binary

} // namespace io
