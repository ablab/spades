//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/pipeline/stage.hpp"
#include "common/sequence/range.hpp"
#include <set>
struct Variation;
struct Stats;

namespace debruijn_graph {

    class WastewaterDisentangle : public spades::AssemblyStage {
    public:
        WastewaterDisentangle()
                : AssemblyStage("Wastewater Disentangle", "wastewater_disentangle") {}

        void run(graph_pack::GraphPack &gp, const char*) override;
        void SelectLinaiges(const std::vector<Stats> &stats, std::vector<Stats> &curated_linaiges, std::set<Variation> &unused_variations);
        std::map<std::string, double> AssignRelativeCoverages(const std::vector<Stats> &curated_linaiges, const std::vector<double> &alpha_coverages, const std::vector<std::map<char, double>> &meaningful_coverages);
    };

}
