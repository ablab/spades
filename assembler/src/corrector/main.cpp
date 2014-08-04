//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

// WTF: What is 'assembler' ?
/*
 * Assembler Main
 */
// WTF: get rid of include-all-you-can-ever-imagine-header
#include "include.hpp"
//#include "config_struct.hpp"


// WTF: graph_pack? construction? stage?
#include "graph_pack.hpp"
#include "construction.hpp"
#include "stage.hpp"

// WTF: make proper order of includes
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "read.hpp"
// WTF: never use relative stuff
#include "../include/segfault_handler.hpp"
#include "sam_reader.hpp"
#include "dataset_processor.hpp"
#include "config_struct.hpp"

// FIXME: EVERYWHERE: USE SPACES, NOT TABS! FIX ALL THE CODING STYLE PROBLEMS EVERYWHERE

void create_console_logger() {
  using namespace logging;


  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}


int main(int /*argc*/, char** argv) {
  perf_counter pc;

  srand(42);
  srandom(42);

  create_console_logger();

  // WTF: Check argc
  
  string contig_name (argv[2]);
  string cfg_file(argv[1]);
  corr_cfg::create_instance(cfg_file);
  string work_dir = corr_cfg::get().work_dir;
  if (! path::check_existence(corr_cfg::get().output_dir))
	  path::make_dir(corr_cfg::get().output_dir);
  if (! path::check_existence(corr_cfg::get().work_dir))
	  path::make_dir(corr_cfg::get().work_dir);

  corrector::DatasetProcessor dp(contig_name);
  dp.ProcessDataset();
  unsigned ms = (unsigned)pc.time_ms();
  unsigned secs = (ms / 1000) % 60;
  unsigned mins = (ms / 1000 / 60) % 60;
  unsigned hours = (ms / 1000 / 60 / 60);
  
  INFO("Correcting time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

  // OK
  return 0;
}
