//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "paired_info/paired_info.hpp"

namespace io {

namespace binary {

template<typename Index>
class PairedIndexIO : public IOSingleDefault<Index> {
public:
    PairedIndexIO()
            : IOSingleDefault<Index>("paired index", ".prd") {
    }
};

template<typename G, typename Traits, template<typename, typename> class Container>
struct IOTraits<omnigraph::de::PairedIndex<G, Traits, Container>> {
    typedef PairedIndexIO<omnigraph::de::PairedIndex<G, Traits, Container>> Type;
};

template<typename Index>
class PairedIndicesIO : public IOCollection<omnigraph::de::PairedIndices<Index>> {
public:
    PairedIndicesIO()
            : IOCollection<omnigraph::de::PairedIndices<Index>>(
                    std::unique_ptr<IOSingle<Index>>(new PairedIndexIO<Index>())) {
    }
};

template<typename Index>
struct IOTraits<omnigraph::de::PairedIndices<Index>> {
    typedef PairedIndicesIO<Index> Type;
};

} // namespace binary

} // namespace io
