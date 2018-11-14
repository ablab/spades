//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "modules/alignment/long_read_storage.hpp"

namespace io {

namespace binary {

template<typename Graph>
class LongReadsIO : public IOCollection<debruijn_graph::LongReadContainer<Graph>> {
public:
    typedef debruijn_graph::LongReadContainer<Graph> Type;
    typedef typename Type::value_type SingleType;
    LongReadsIO()
            : IOCollection<Type>(std::unique_ptr<IOSingle<SingleType>>(
                    new IOSingleDefault<SingleType>("path storage", ".mpr"))) {
    }
};

template<typename Graph>
struct IOTraits<debruijn_graph::LongReadContainer<Graph>> {
    typedef LongReadsIO<Graph> Type;
};

} // namespace binary

} // namespace io
