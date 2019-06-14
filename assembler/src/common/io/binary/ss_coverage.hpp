//***************************************************************************
//* Copyright (c) 2018-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "modules/alignment/rna/ss_coverage.hpp"

namespace io {

namespace binary {

class SSCoverageIO : public IOCollection<debruijn_graph::SSCoverageContainer> {
public:
    typedef debruijn_graph::SSCoverageContainer Type;
    typedef typename Type::value_type SingleType;
    SSCoverageIO()
            : IOCollection<Type>(std::unique_ptr<IOSingle<SingleType>>(
                    new IOSingleDefault<SingleType>("ss coverage", ".sscvr"))) {
    }
};

template<>
struct IOTraits<debruijn_graph::SSCoverageContainer> {
    typedef SSCoverageIO Type;
};

} // namespace binary

} // namespace io
