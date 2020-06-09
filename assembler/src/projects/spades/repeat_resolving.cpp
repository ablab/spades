//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/logger.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "modules/path_extend/pipeline/launcher.hpp"

#include "repeat_resolving.hpp"

namespace debruijn_graph {

static void PEResolving(GraphPack& gp) {
    path_extend::PathExtendParamsContainer params(cfg::get().ds,
                                                  cfg::get().pe_params,
                                                  cfg::get().ss,
                                                  cfg::get().output_dir,
                                                  cfg::get().mode,
                                                  cfg::get().uneven_depth,
                                                  cfg::get().avoid_rc_connections,
                                                  cfg::get().use_scaffolder);

    path_extend::PathExtendLauncher exspander(cfg::get().ds, params, gp);
    exspander.Launch();
}

static bool HasValidLibs() {
    for (const auto& lib : cfg::get().ds.reads) {
        if (!lib.is_repeat_resolvable())
            continue;

        if (!lib.is_paired() ||
            !math::eq(lib.data().mean_insert_size, 0.0)) {
            return true;
        }
    }
    
    return false;
}


void RepeatResolution::load(GraphPack &gp, const std::string &load_from, const char *prefix) {
    // Do nothing in final mode, otherwise - produce saves
    if (!preliminary_)
        return;

    AssemblyStage::load(gp, load_from, prefix);
}


void RepeatResolution::save(const GraphPack &gp, const std::string &save_to, const char *prefix) const {
    // Do nothing in final mode, otherwise - produce saves
    if (!preliminary_)
        return;

    AssemblyStage::save(gp, save_to, prefix);
}

void RepeatResolution::run(GraphPack &gp, const char*) {
    if (cfg::get().developer_mode)
        stats::PrepareForDrawing(gp);

    visualization::graph_labeler::DefaultLabeler<Graph> labeler(gp.get<Graph>(), gp.get<EdgesPositionHandler<Graph>>());
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);
    printer(config::info_printer_pos::before_repeat_resolution);

    //todo awful hack to get around PE using cfg::get everywhere...
    //Is it possible to fix this problem now or still too soon?
    auto tmp_params_storage = cfg::get().pe_params;
    if (preliminary_) {
        INFO("Setting up preliminary path extend settings")
        cfg::get_writable().pe_params = *cfg::get().prelim_pe_params;
    }

    bool no_valid_libs = !HasValidLibs();
    bool use_single_reads = cfg::get().use_single_reads;
    if (cfg::get().rr_enable && no_valid_libs && !use_single_reads)
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");

    if ((no_valid_libs || cfg::get().rm == config::resolving_mode::none) && !use_single_reads) {
        return;
    }
    if (cfg::get().rm == config::resolving_mode::path_extend) {
        INFO("Using Path-Extend repeat resolving");
        PEResolving(gp);
    } else {
        INFO("Unsupported repeat resolver");
    }
    if (preliminary_) {
        INFO("Restoring initial path extend settings")
        cfg::get_writable().pe_params = tmp_params_storage;
    }
}



} // debruijn_graph
