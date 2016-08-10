//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dipspades_config.hpp"
#include "pipeline/config_common.hpp"
#include "utils/files_utils.hpp"
#include "utils/path_helper.hpp"

using namespace dipspades;

void load(dipspades_config::base_params &bp,
        boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(bp.K                    , pt,    "K"                        );
    load(bp.max_memory            , pt,    "max_memory"            );
    load(bp.max_threads            , pt,     "max_threads"            );
    load(bp.read_buffer_size     , pt,     "read_buffer_size"        );
}

void load(dipspades_config::io_params &io,
        boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(io.haplocontigs        , pt,     "haplocontigs"            );
    io.num_libraries = GetAllLinesFromFile(io.haplocontigs).size();

    load(io.log_filename        , pt,     "log_filename"            );

    load(io.output_base            , pt,     "output_base"            );
    if (io.output_base[io.output_base.length() - 1] != '/')
        io.output_base += '/';

    load(io.output_dir            , pt,     "output_dir"            );
    if (io.output_dir[io.output_dir.length() - 1] != '/')
        io.output_dir += '/';

    load(io.tmp_dir            , pt,   "tmp_dir"                );
    if (io.tmp_dir[io.tmp_dir.length() - 1] != '/')
        io.tmp_dir += '/';

    load(io.saves        , pt, "saves"        );
    if(io.saves[io.saves.length() - 1] != '/')
        io.saves += '/';
}

void load(dipspades_config::run_params &rp,
        boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(rp.entry_point         , pt,         "entry_point"    );
    load(rp.developer_mode      , pt,        "developer_mode");
}

void edit_io_params(bool developer_mode, dipspades_config::io_params &io){
    if(developer_mode){
        io.dataset_name = io.output_dir.substr(0, io.output_dir.length() - 1);
        io.output_dir = io.output_base + io.output_dir + "/";
        io.output_root = io.output_dir;
        io.output_suffix = path::MakeLaunchTimeDirName() + "/";
        io.output_dir = io.output_root + io.output_suffix;
        io.output_saves = io.output_dir + "saves/";
//        io.load_from = io.output_root + io.load_from;
        if (io.tmp_dir[0] != '/') { // relative path
            io.tmp_dir = io.output_dir + io.tmp_dir;
        }
        return;
    }

    // no developer mode
    io.dataset_name = io.output_dir;
    io.output_root = io.output_dir;
    io.output_suffix = "";
    io.output_base = "";
    io.output_saves = io.output_dir;
    io.saves = "";
    if (io.tmp_dir[0] != '/') { // relative path
        io.tmp_dir = io.output_dir + io.tmp_dir;
    }
}

inline void load(dipspades_config::polymorphic_br &pbr,
        boost::property_tree::ptree const& pt, bool){
    using config_common::load;
    load(pbr.enabled                        , pt, "enabled"                        );
    load(pbr.rel_bulge_length                , pt, "rel_bulge_length"            );
    load(pbr.rel_bulge_align                , pt, "rel_bulge_align"                );
    load(pbr.paired_vert_abs_threshold        , pt, "paired_vert_abs_threshold"    );
    load(pbr.paired_vert_rel_threshold        , pt, "paired_vert_rel_threshold"    );
    load(pbr.max_bulge_nucls_len            , pt, "max_bulge_nucls_len"            );
    load(pbr.max_neigh_number                , pt, "max_neigh_number"            );
    load(pbr.num_iters_lbr                    , pt, "num_iters_lbr"                );
}

inline void load(dipspades_config::consensus_constructor &cc,
        boost::property_tree::ptree const& pt, bool /*complete*/){
    using config_common::load;
    load(cc.enabled                        , pt, "enabled"                     );
    load(cc.bulge_len_quantile            , pt, "bulge_len_quantile"            );
    load(cc.tails_lie_on_bulges            , pt, "tails_lie_on_bulges"            );
    load(cc.estimate_tails                , pt, "estimate_tails"                );
    load(cc.align_bulge_sides            , pt, "align_bulge_sides"            );
    load(cc.min_overlap_size            , pt, "min_overlap_size"            );
    load(cc.min_lcs_size                , pt, "min_lcs_size"                );
    load(cc.max_loop_length             , pt, "max_loop_length"                );
}

inline void load(dipspades_config::haplotype_assembly &ha,
        boost::property_tree::ptree const& pt, bool /*complete*/){
    using config_common::load;
    load(ha.ha_enabled        , pt,     "ha_enabled"    );
}

void load(dipspades_config &cfg,
        boost::property_tree::ptree const &pt, bool complete){
      using config_common::load;
      load(cfg.bp                , pt,        "bp", complete);
      load(cfg.io                , pt,        "io", complete);
      load(cfg.rp                 , pt,        "rp", complete);
      load(cfg.cc                , pt,        "cc", complete);
      load(cfg.ha                , pt,         "ha", complete);
      load(cfg.pbr                 , pt,         "pbr", complete);
}

void load(dipspades_config &cfg, std::string const &filename) {
      boost::property_tree::ptree pt;
      boost::property_tree::read_info(filename, pt);
      load(cfg, pt, true);
      edit_io_params(cfg.rp.developer_mode, cfg.io);

}
