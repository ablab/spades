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

namespace binary {

template<typename Graph>
class LongReadsIO : public IOCollection<debruijn_graph::LongReadContainer<Graph>, EdgeMapper<Graph>> {
public:
    typedef debruijn_graph::LongReadContainer<Graph> Type;
    typedef typename Type::value_type SingleType;
    typedef EdgeMapper<Graph> Mapper;
    LongReadsIO()
            : IOCollection<Type, Mapper>(std::unique_ptr<IOSingle<SingleType, Mapper>>(
                    new IOSingleDefault<SingleType, Mapper>("path storage", ".mpr"))) {
    }
};

template<typename Graph>
struct IOTraits<debruijn_graph::LongReadContainer<Graph>> {
    typedef LongReadsIO<Graph> Type;
};

} // namespace binary

} // namespace io
