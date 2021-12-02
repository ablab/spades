//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace spades {

class ReadConversion : public AssemblyStage {
public:
    ReadConversion()
            : AssemblyStage("Binary Read Conversion", "read_conversion") { }

    void run(graph_pack::GraphPack &, const char *) override;
    void load(graph_pack::GraphPack &, const std::filesystem::path &load_from, const char *prefix = nullptr) override;
    void save(const graph_pack::GraphPack &, const std::filesystem::path &save_to, const char *prefix = nullptr) const override;
    bool run_on_load() const override { return true; }
};

}
