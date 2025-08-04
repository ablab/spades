//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "adt/pack.hpp"

#include <boost/noncopyable.hpp>

#include <filesystem>

namespace graph_pack {

class GraphPack: public adt::pack, private boost::noncopyable {
public:
    GraphPack(unsigned k,
              const std::filesystem::path &workdir_, size_t lib_count,
              const std::vector<std::string> &genome = std::vector<std::string>(0),
              size_t flanking_range = 50,
              size_t max_mapping_gap = 0,
              size_t max_gap_diff = 0,
              bool detach_indices = true);

    unsigned k() const { return k_; }
    const std::filesystem::path &workdir() const { return workdir_; }

    void DetachAll();

private:
    unsigned k_;
    std::filesystem::path workdir_;
};

} // namespace graph_pack
