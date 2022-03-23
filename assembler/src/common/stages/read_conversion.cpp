//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "read_conversion.hpp"
#include "io/dataset_support/read_converter.hpp"

namespace spades {

void ReadConversion::run(graph_pack::GraphPack &, const char *) {
    io::ConvertIfNeeded(cfg::get_writable().ds.reads,
                        cfg::get().max_threads);
}

void ReadConversion::load(graph_pack::GraphPack &,
                         const std::string &load_from,
                         const char* prefix) {
    std::string p = fs::append_path(load_from, prefix == NULL ? id() : prefix);
    INFO("Loading current state from " << p);

    debruijn_graph::config::load_lib_data(p);

    io::ConvertIfNeeded(cfg::get_writable().ds.reads,
                        cfg::get().max_threads);
}

void ReadConversion::save(const graph_pack::GraphPack &,
                         const std::string &save_to,
                         const char* prefix) const {
    std::string p = fs::append_path(save_to, prefix == NULL ? id() : prefix);
    INFO("Saving current state to " << p);

    debruijn_graph::config::write_lib_data(p);
}

} // namespace spades
