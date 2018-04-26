//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/id_mapper.hpp"
#include "io_base.hpp"

namespace io {

template<typename Index>
class PairedIndexIO : public IOBase<Index> {
public:
    PairedIndexIO(const IdMapper<typename Index::EdgeId> &mapper):
            IOBase<Index>("paired index", ".prd"), mapper_(mapper) {}
protected:
    void SaveImpl(SaveFile &file, const Index &index) override {
        file << index;
    }

    void LoadImpl(LoadFile &file, Index &index) override {
        index.BinRead(file.str(), mapper_);
    }
private:
    const IdMapper<typename Index::EdgeId> &mapper_;
};

}
