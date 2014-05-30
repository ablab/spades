//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * Assembler Main
 */
#include "standard.hpp"
#include "logger/log_writers.hpp"

#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "launch.hpp"
#include "memory_limit.hpp"
#include "copy_file.hpp"
#include "perfcounter.hpp"
#include "runtime_k.hpp"

#include "config_struct.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void link_output(std::string const& link_name) {
  if (!cfg::get().run_mode)
    return;

  std::string link = cfg::get().output_root + link_name;
  unlink(link.c_str());
  if (symlink(cfg::get().output_suffix.c_str(), link.c_str()) != 0)
    WARN( "Symlink to \"" << link << "\" launch failed");
}

void link_previous_run(std::string const& previous_link_name, std::string const& link_name) {
  if (!cfg::get().run_mode)
    return;

  char buf[255];

  std::string link = cfg::get().output_dir + previous_link_name;
  unlink(link.c_str());
  ssize_t count = readlink((cfg::get().output_root + link_name).c_str(), buf, sizeof(buf) - 1);
  if (count >= 0){
    buf[count] = '\0';
    std::string previous_run("../");
    previous_run = previous_run + buf;
    if (symlink(previous_run.c_str(), link.c_str()) != 0) {
      DEBUG( "Symlink to \"" << link << "\" launch failed : " << previous_run);
    }
  } else {
    DEBUG( "Symlink to \"" << link << "\" launch failed");
  }
}

struct on_exit_output_linker {
  on_exit_output_linker(std::string const& link_name) :
      link_name_(link_name) { }

  ~on_exit_output_linker() {
    link_previous_run("previous", link_name_);
    link_output(link_name_);
  }

 private:
  std::string link_name_;
};

void copy_configs(string cfg_filename, string to) {
  if (!cfg::get().run_mode)
    return;

  using namespace debruijn_graph;

  if (!make_dir(to)) {
    WARN("Could not create files use in /tmp directory");
  }
  path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".info", true);
}

void load_config(string cfg_filename) {
  path::CheckFileExistenceFATAL(cfg_filename);

  cfg::create_instance(cfg_filename);

  if (!cfg::get().project_name.empty()) {
    make_dir(cfg::get().output_base + cfg::get().project_name);
  }

  make_dir(cfg::get().output_root);
  make_dir(cfg::get().tmp_dir);

  make_dir(cfg::get().output_dir);
  if (cfg::get().developer_mode)
    make_dir(cfg::get().output_saves);

  make_dir(cfg::get().temp_bin_reads_path);

  string path_to_copy = path::append_path(cfg::get().output_dir, "configs");
  copy_configs(cfg_filename, path_to_copy);
}

void create_console_logger(string cfg_filename) {
  using namespace logging;

  string log_props_file = cfg::get().log_filename;

  if (!path::FileExists(log_props_file))
    log_props_file = path::append_path(path::parent_path(cfg_filename), cfg::get().log_filename);

  logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

int main(int /*argc*/, char** argv) {
  perf_counter pc;

  const size_t GB = 1 << 30;

  segfault_handler sh(bind(link_output, "latest"));

  srand(42);
  srandom(42);

  try {
    using namespace debruijn_graph;

    string cfg_filename = argv[1];

    load_config          (cfg_filename);
    create_console_logger(cfg_filename);

    on_exit_output_linker try_linker("latest");

    VERIFY(cfg::get().K >= runtime_k::MIN_K && cfg::get().K < runtime_k::MAX_K);
    VERIFY(cfg::get().K % 2 != 0);

    // read configuration file (dataset path etc.)

    limit_memory(cfg::get().max_memory * GB);

    // assemble it!
    INFO("Assembling dataset (" << cfg::get().dataset_file << ") with K=" << cfg::get().K);

    spades::assemble_genome();

    link_output("latest_success");
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
