#include <iostream>
#include <string>
#include <map>
#include <exception>

#include "sequence/seq.hpp"
#include "logger/log_writers.hpp"
#include "memory_limit.hpp"
#include "running_modes.hpp"
#include "config_struct_cclean.hpp"
#include "simple_tools.hpp"
#include "adapter_index.hpp"

#include "valid_kmer_generator.hpp"
#include "io/read_processor.hpp"
#include "ssw/ssw_cpp.h"
#include "additional.cpp"

using logging::logger;
using logging::create_logger;
using logging::console_writer;
using std::string;
using additional::WORK_MODE_TYPE;
using additional::SIMPLE;
using additional::BRUTE_SIMPLE;
using additional::BRUTE_DEEP;

// Terminal args constants, say no magic numbers
constexpr int ARG_CONFIG_FILE = 2;
constexpr int ARG_DB_FILE = 4;
constexpr int ARG_INPUT_FILE = 3;
constexpr int ARG_MODE = 1;
constexpr int ARGS_MIN = 4;
constexpr int ARGS_MAX = 5;

void usage() {
  std::cout << "usage: cclean <mode> <config-file> <database> <input-file>"
            << std::endl;
  std::cout << "\n<config-file>\t Path to config file" << std::endl;
  std::cout << "<database>\t Path to database file" << std::endl;
  std::cout << "<input-file>\t Path to input file. Currently only .gz supported"
            << std::endl;
  std::cout << "<mode>:\t\t<simple> - simple not pair-end fast algorithm."
            << std::endl;
  std::cout << "\t\t<bruteforce> - simple \"fast\" bruteforce mode."
            << std::endl;
  std::cout << "\t\t<bruteforce:deep> - slow mode with analyze of read quality."
            << std::endl;
  std::cout << "\t\tNotice, that in both modes used threeshold from config."
            << std::endl;
  std::cout << "\t\tBut in simple mode it is threeshold of best align score."
            << std::endl;
  std::cout << "\t\tAnd in deep mode it is threeshold of score based on formula: "
            << std::endl;
  std::cout << "\t\tMismatch -Q/10 score, match +0.6 score, where Q - read quality."
            << std::endl;
}

void create_console_logger() {
  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

int main(int argc, char *argv[]) {
  try {
    if (argc < ARGS_MIN || argc > ARGS_MAX) {
      usage();
      return EXIT_SUCCESS;
    }

    string mode_string = argv[ARG_MODE];
    WORK_MODE_TYPE mode;
    // Brutforce mode enabled?
    if (mode_string == "simple")
      mode = SIMPLE;
    else if (mode_string == "bruteforce")
      mode = BRUTE_SIMPLE;
    else if (mode_string == "bruteforce:deep")
      mode = BRUTE_DEEP;
    else {
      usage();
      return EXIT_SUCCESS;
    }

    create_console_logger();
    clock_t start = clock();

    std::string config_file = argv[ARG_CONFIG_FILE];
    INFO("Loading config from " << config_file.c_str());
    if (!FileExists(config_file)) {
        ERROR("File " + config_file + " doesn't exists.");
        return EXIT_SUCCESS;
    }
    cfg::create_instance(config_file);

    // hard memory limit
    // const size_t GB = 1 << 30;
    // limit_memory(cfg::get().general_hard_memory_limit * GB);

    const std::string data_base(argv[ARG_DB_FILE]);
    if (!FileExists(data_base)) {
        ERROR("File " + data_base + " doesn't exists.");
        return EXIT_SUCCESS;
    }
    const std::string input_file(argv[ARG_INPUT_FILE]);
    if (!FileExists(input_file)) {
        ERROR("File " + input_file + " doesn't exists.");
        return EXIT_SUCCESS;
    }

    cclean::AdapterIndex index;
    cclean::AdapterIndexBuilder().FillAdapterIndex(data_base, index);

    ireadstream *input;

    INFO("Init file with reads-to-clean from " << input_file << " ... ");
    input = new ireadstream(input_file);

    INFO("Start matching reads against database ...");
    std::ofstream output(cfg::get().output_file);
    std::ofstream bed(cfg::get().output_bed);
    if (!output.is_open() || !bed.is_open()) {
      ERROR("Cannot open output file: " << cfg::get().output_file << " or "
                                        << cfg::get().output_bed);
      return EXIT_SUCCESS;
    }
    // Main work wrapper
    exactAndAlign(output, bed, input, data_base, index, mode);

    output.close();
    bed.close();

    delete input;

    clock_t ends = clock();
    INFO("Processor Time Spent: " << (double) (ends - start) / CLOCKS_PER_SEC
         << " seconds.");
    INFO("Goodbye!");
  }
  catch (std::exception& e) {
    ERROR("Error: " << e.what());
    return EXIT_SUCCESS;
  }
  catch (...) {
    ERROR("Something unexpected happened!");
    return EXIT_SUCCESS;
  }

  return EXIT_SUCCESS;
}
