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

namespace binary {

template<typename Index>
class PairedIndexIO : public IOSingleDefault<Index, EdgeMapper<Index>> {
public:
    typedef EdgeMapper<Index> Mapper;
    PairedIndexIO()
            : IOSingleDefault<Index, Mapper>("paired index", ".prd") {
    }
};

template<typename G, typename Traits, template<typename, typename> class Container>
struct IOTraits<omnigraph::de::PairedIndex<G, Traits, Container>> {
    typedef PairedIndexIO<omnigraph::de::PairedIndex<G, Traits, Container>> Type;
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

template<typename Index>
struct IOTraits<omnigraph::de::PairedIndices<Index>> {
    typedef PairedIndicesIO<Index> Type;
};

} // namespace binary

} // namespace io
