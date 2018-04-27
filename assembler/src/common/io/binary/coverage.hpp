//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "io/id_mapper.hpp"

namespace io {

template<template<typename> typename IndexT, typename Graph>
class CoverageIO : public IOBase<IndexT<Graph>> {
public:
    typedef IndexT<Graph> Type;
    typedef IdMapper<typename Graph::EdgeId> Mapper;
    CoverageIO(const Mapper &mapper):
            IOBase<Type>("coverage", ".cvr"), mapper_(mapper) {}
private:
    void SaveImpl(SaveFile &file, const Type &index) {
        for (auto it = index.g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            auto e = *it;
            file << e.int_id() << index.RawCoverage(e);
        }
    }

    void LoadImpl(LoadFile &file, Type &index) {
        while (file) { //Read until the end
            auto eid = mapper_[file.Read<size_t>()];
            auto cov = file.Read<unsigned>();
            index.SetRawCoverage(eid, cov);
        }
    }

    const Mapper &mapper_;
};

}
