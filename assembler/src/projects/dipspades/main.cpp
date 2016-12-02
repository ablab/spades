//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * Assembler Main
 */
#include "utils/logger/log_writers.hpp"

#include "utils/segfault_handler.hpp"
#include "utils/memory_limit.hpp"
#include "utils/copy_file.hpp"

#include "pipeline/graph_pack.hpp"
#include "stages/construction.hpp"

#include "version.hpp"

#include "dipspades.hpp"

void make_dirs(){
    make_dir(dsp_cfg::get().io.output_base);
    make_dir(dsp_cfg::get().io.output_root);
    make_dir(dsp_cfg::get().io.output_dir);
    make_dir(dsp_cfg::get().io.output_saves);
    make_dir(dsp_cfg::get().io.tmp_dir);
}

void copy_configs(string cfg_filename, string to) {
  using namespace debruijn_graph;

  if (!make_dir(to)) {
    WARN("Could not create files use in /tmp directory");
  }
  path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".info", true);
}

void load_config(string cfg_filename) {
  path::CheckFileExistenceFATAL(cfg_filename);
  dsp_cfg::create_instance(cfg_filename);
//  string path_to_copy = path::append_path(dsp_cfg::get().io.output_dir, "configs");
//  copy_configs(cfg_filename, path_to_copy);
}

void create_console_logger(string cfg_filename) {
  using namespace logging;

  string log_props_file = dsp_cfg::get().io.log_filename;

  if (!path::FileExists(log_props_file)){
    log_props_file = path::append_path(path::parent_path(cfg_filename), dsp_cfg::get().io.log_filename);
  }

  logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

int main(int /*argc*/, char** argv) {
  perf_counter pc;
  const size_t GB = 1 << 30;

  srand(42);
  srandom(42);

  segfault_handler sh;

  try {
    using namespace debruijn_graph;
    string cfg_filename = argv[1];
    load_config          (cfg_filename);
    make_dirs();
    if(dsp_cfg::get().rp.developer_mode)
        copy_configs(cfg_filename, path::append_path(dsp_cfg::get().io.output_dir, "configs"));
    create_console_logger(cfg_filename);

    INFO("Loaded config from " << cfg_filename);
    
    VERIFY(dsp_cfg::get().bp.K >= runtime_k::MIN_K && dsp_cfg::get().bp.K < runtime_k::MAX_K);
    VERIFY(dsp_cfg::get().bp.K % 2 != 0);

    limit_memory(dsp_cfg::get().bp.max_memory * GB);

    INFO("Starting dipSPAdes, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);
    INFO("Assembling dataset (" << dsp_cfg::get().io.dataset_name << ") with K=" << dsp_cfg::get().bp.K);
    dipspades::run_dipspades();
//    link_output("latest_success");
  } catch (std::bad_alloc const& e) {
    std::cerr << "Not enough memory to run SPAdes. " << e.what() << std::endl;
    return EINTR;
  } catch (std::exception const& e) {
    std::cerr << "Exception caught " << e.what() << std::endl;
    return EINTR;
  } catch (...) {
    std::cerr << "Unknown exception caught " << std::endl;
    return EINTR;
  }

  unsigned ms = (unsigned)pc.time_ms();
  unsigned secs = (ms / 1000) % 60;
  unsigned mins = (ms / 1000 / 60) % 60;
  unsigned hours = (ms / 1000 / 60 / 60);
  INFO("Assembling time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

  // OK
  return 0;
}
