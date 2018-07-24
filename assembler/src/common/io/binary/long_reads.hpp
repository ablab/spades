//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "io/id_mapper.hpp"
#include "modules/alignment/long_read_storage.hpp"

namespace io {

template<typename Graph>
class LongReadsIO : public IOSingle<debruijn_graph::LongReadContainer<Graph>, EdgeMapper<Graph>> {

public:
    typedef typename debruijn_graph::LongReadContainer<Graph> Type;
    typedef EdgeMapper<Graph> Mapper;

    LongReadsIO()
            : IOSingle<Type, Mapper>("long reads storage", ".mpr") {
    }

private:
    void SaveImpl(SaveFile &file, const Type &container) override {
        for (size_t i = 0; i < container.size(); ++i)
            file << container[i];
    }

    void LoadImpl(LoadFile &file, Type &container, const Mapper &mapper) override {
        //TODO: why and how is the container size set before loading?
        for (size_t i = 0; i < container.size(); ++i)
            container[i].BinRead(file.stream(), mapper);
    }
};

}
