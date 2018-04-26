//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "io/id_mapper.hpp"
#include "paired_info/paired_info.hpp"

namespace io {

template<typename Index>
class PairedIndexIO : public IOBase<Index> {
public:
    typedef IdMapper<typename Index::EdgeId> Mapper;
    PairedIndexIO(const Mapper &mapper):
            IOBase<Index>("paired index", ".prd"), mapper_(mapper) {}
private:
    void SaveImpl(SaveFile &file, const Index &index) override {
        file << index;
    }

    void LoadImpl(LoadFile &file, Index &index) override {
        index.BinRead(file.str(), mapper_);
    }

    const Mapper &mapper_;
};

template<typename Index>
class PairedIndicesIO {
public:
    typedef omnigraph::de::PairedIndices<Index> Type;
    typedef typename PairedIndexIO<Index>::Mapper Mapper;

    PairedIndicesIO(const Mapper &mapper):
        mapper_(mapper) {}

    void Save(const std::string &basename, const Type &indices) {
        PairedIndexIO<Index> io(mapper_);
        for (size_t i = 0; i < indices.size(); ++i) {
            io.Save(basename + "_" + std::to_string(i), indices[i]);
        }
    }

    void Load(const std::string &basename, Type &indices) {
        PairedIndexIO<Index> io(mapper_);
        for (size_t i = 0; i < indices.size(); ++i) {
            io.Load(basename + "_" + std::to_string(i), indices[i]);
        }
    }
private:
    const Mapper &mapper_;
};

}
