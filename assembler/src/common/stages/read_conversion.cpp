//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "read_conversion.hpp"
#include "io/dataset_support/read_converter.hpp"

namespace spades {

void ReadConversion::run(debruijn_graph::conj_graph_pack &, const char *) {
    for (auto &lib : cfg::get_writable().ds.reads) {
        if (!io::ReadConverter::LoadLibIfExists(lib))
            io::ReadConverter::ConvertToBinary(lib);
    }
}

} // namespace spades
