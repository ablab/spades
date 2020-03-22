//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "common/assembly_graph/paths/bidirectional_path.hpp"

namespace io {

namespace binary {


class TrustedPathsIO : public IOCollection<std::vector<path_extend::BidirectionalPathStorage>> {
public:
    using Type = std::vector<path_extend::BidirectionalPathStorage>;
    using SingleType = path_extend::BidirectionalPathStorage;
    TrustedPathsIO()
            : IOCollection<Type>(std::unique_ptr<IOSingle<SingleType>>(
                    new IOSingleDefault<SingleType>("trusted paths storage", ".bdp"))) {
    }
};

template<>
struct IOTraits<std::vector<path_extend::BidirectionalPathStorage>> {
    using Type = TrustedPathsIO;
};

} // namespace binary

} // namespace io
