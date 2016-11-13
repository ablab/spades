//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <iostream>
#include <string>
#include <map>
#include <exception>

#include "sequence/seq.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils/memory_limit.hpp"
#include "running_modes.hpp"
#include "config_struct_cclean.hpp"
#include "utils/simple_tools.hpp"
#include "adapter_index.hpp"
#include "utils.hpp"

#include "valid_kmer_generator.hpp"
#include "io/read_processor.hpp"
#include "modules/ssw_cpp.h"
#include "additional.cpp"

#include "job_wrappers.hpp"
#include "brute_force_clean.hpp"

using logging::logger;
using logging::create_logger;
using logging::console_writer;
using std::string;

constexpr int CONFIG_FILE_ARG = 1;

void usage() {
  std::cout << "usage: cclean [program config file]" << std::endl;
}

void create_console_logger() {
  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

int main(int argc, char *argv[]) {

  create_console_logger();

  if (argc < 2) {
    usage();
    return EXIT_FAILURE;
  }

  std::string config_file = argv[CONFIG_FILE_ARG];
  INFO("Loading config from " << config_file.c_str());
  if (!path::FileExists(config_file)) {
      ERROR("File " + config_file + " doesn't exists.");
      return EXIT_FAILURE;
  }
  cfg::create_instance(config_file);

  const std::string &database = cfg::get().database;
  if (!path::FileExists(database)) {
      ERROR("File " + database + " doesn't exists.");
      return EXIT_FAILURE;
  }
  const std::string &dataset = cfg::get().dataset_file_name;
  if (!path::FileExists(dataset)) {
      ERROR("File " + dataset + " doesn't exists.");
      return EXIT_FAILURE;
  }

  clock_t start = clock();

  Cleaner::ProcessDataset();  // Main work here

  INFO("DONE");
  clock_t ends = clock();
  INFO("Processor Time Spent: " << (double) (ends - start) / CLOCKS_PER_SEC
       << " seconds.");
  INFO("Goodbye!");

  return EXIT_SUCCESS;
}
