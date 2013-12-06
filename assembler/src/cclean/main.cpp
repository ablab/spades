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
using additional::WorkModeType;
using additional::SIMPLE;
using additional::BRUTE_SIMPLE;
using additional::BRUTE_DEEP;

// Terminal args constants, say no magic numbers
constexpr int ArgConfigFile = 2;
constexpr int ArgDbFile = 5;
constexpr int ArgInputFile = 3;
constexpr int ArgOutputFile = 4;
constexpr int ArgMode = 1;
constexpr int ArgsMin = 4;
constexpr int ArgsMax = 6;

void usage() {
  std::cout << "usage: cclean <mode> <config> <input> <output> <database>"
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
    if (argc < ArgsMin || argc > ArgsMax) {
      usage();
      return EXIT_SUCCESS;
    }

    string mode_string = argv[ArgMode];
    WorkModeType mode;
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

    std::string config_file = argv[ArgConfigFile];
    INFO("Loading config from " << config_file.c_str());
    if (!FileExists(config_file)) {
        ERROR("File " + config_file + " doesn't exists.");
        return EXIT_SUCCESS;
    }
    cfg::create_instance(config_file);

    // hard memory limit
    // const size_t GB = 1 << 30;
    // limit_memory(cfg::get().general_hard_memory_limit * GB);

    const std::string data_base(argv[ArgDbFile]);
    if (!FileExists(data_base)) {
        ERROR("File " + data_base + " doesn't exists.");
        return EXIT_SUCCESS;
    }
    const std::string input_file(argv[ArgInputFile]);
    if (!FileExists(input_file)) {
        ERROR("File " + input_file + " doesn't exists.");
        return EXIT_SUCCESS;
    }
    ireadstream *input;
    INFO("Init file with reads-to-clean from " << input_file << " ... ");
    input = new ireadstream(input_file);

    std::string output_file = argv[ArgOutputFile];
    std::ofstream output(output_file);

    INFO("Start matching reads against database ...");
    std::ofstream aligned_output(cfg::get().output_file);
    std::ofstream bed(cfg::get().output_bed);
    if (!aligned_output.is_open() || !bed.is_open() || !output.is_open()) {
      ERROR("Cannot open output file: " << cfg::get().output_file << " or "
            << cfg::get().output_bed << "or " << output_file);
      return EXIT_SUCCESS;
    }

    cclean::AdapterIndex index;
    cclean::AdapterIndexBuilder().FillAdapterIndex(data_base, index);
    // Main work wrapper
    ExactAndAlign(aligned_output, bed, input, output, data_base, index, mode);

    aligned_output.close();
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
