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
class PairedIndexIO : public IOSingle<Index, EdgeMapper<Index>> {
public:
    typedef EdgeMapper<Graph> Mapper;
    PairedIndexIO()
            : IOSingle<Index, Mapper>("paired index", ".prd") {
    }

protected:
    void SaveImpl(SaveFile &file, const Index &index) override {
        file << index;
    }

    void LoadImpl(LoadFile &file, Index &index, const Mapper &mapper) override {
        index.BinRead(file.stream(), mapper);
    }
};

template<typename Index>
class PairedIndicesIO : public IOCollection<omnigraph::de::PairedIndices<Index>, EdgeMapper<Index>> {
public:
    typedef EdgeMapper<Index> Mapper;
    PairedIndicesIO()
            : IOCollection<omnigraph::de::PairedIndices<Index>, Mapper>(
                    std::unique_ptr<IOSingle<Index, Mapper>>(new PairedIndexIO<Index>())) {
    }
};

}
