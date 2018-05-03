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
class PairedIndexIO : public IOSingle<Index> {
public:
    typedef IdMapper<typename Index::EdgeId> Mapper;
    PairedIndexIO(const Mapper &mapper)
            : IOSingle<Index>("paired index", ".prd"), mapper_(mapper) {
    }

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
class PairedIndicesIO : public IOCollection<Index> {
public:
    PairedIndicesIO(const typename PairedIndexIO<Index>::Mapper &mapper)
            : IOCollection<Index>(new PairedIndexIO<Index>(mapper)) {
    }
};

}
