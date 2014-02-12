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
#include "utils.hpp"

#include "valid_kmer_generator.hpp"
#include "io/read_processor.hpp"
#include "ssw/ssw_cpp.h"
#include "additional.cpp"

using logging::logger;
using logging::create_logger;
using logging::console_writer;
using std::string;
using additional::WorkModeType;
using additional::SINGLE_END;
using additional::SINGLE_END_Q;
using additional::BRUTE_SIMPLE;
using additional::BRUTE_WITH_Q;

void usage() {
  std::cout << "usage: cclean [--mode] [--config] [--input] [--output] [--database]"
            << std::endl;
  std::cout << "Commands description:" << std::endl;
  std::cout << "\tmode [or m]\tScaning reads mode. Variants:\n\t\t\tBSE - Brute " <<
               "Single End mode\n\t\t\tBSEQ - Brute Single End mode with Quality" <<
               "\n\t\t\tSE - Fast Single End mode\n\t\t\tSEQ - Fast Single End mode with "<<
               "Quality\n\t\t\tPE - Pair End mode\n\t\t\tPEQ - Pair End mode with Quality." <<
               std::endl;
  std::cout << "\tconfig [or c]\tPath to cclean config file." << std::endl;
  std::cout << "\tinput [or i]\tPath to cclean input fastq file." << std::endl;
  std::cout << "\toutput [or o]\tPath to cclean output file." << std::endl;
  std::cout << "\tdatabase [or d]\tPath to adapters database file." << std::endl;
  std::cout << "\tmlen [or ml]\tMinimum lenght for cuted read for output." << std::endl;
  std::cout << "\tinform [or in]\tFull output with aligned and bed reads in specified files." << std::endl;
  std::cout << "\t\t\tNONE - none information will be outputed. (default)" << std::endl;
  std::cout << "\t\t\tFULL - full information will be outputed." << std::endl;
  std::cout << "example:\ncclean --m=BSE --c=config.info.template --i=dataset_500.fastq "
            << "--o=output_bruteforce.fastq --d=TruSeq2-SE.fa --in=FULL" << std::endl;
}

void create_console_logger() {
  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

int main(int argc, char *argv[]) {

  create_console_logger();

  bool correct_args = false;
  std::string error;
  // Processing arguments in map of args - options
  auto options = cclean_utils::ProcessArgs(argc, argv, &correct_args, &error);

  if (!correct_args) {
    ERROR(error);
    usage();
    return EXIT_FAILURE;
  }

  // check min len arg
  if (options.find("mlen") == options.end())
    options["mlen"] = "0";
  if (atoi(options["mlen"].c_str()) < 0) {
    ERROR("Bad value mlen");
    usage();
    return EXIT_FAILURE;
  }

  string mode_string = options["mode"];
  WorkModeType mode;
  // Type mode for processing
  if (mode_string == "SE") // Single end
    mode = SINGLE_END;
  else if (mode_string == "SEQ") // Bruteforce single end
    mode = SINGLE_END_Q;
  else if (mode_string == "BSE") // Bruteforce single end
    mode = BRUTE_SIMPLE;
  else if (mode_string == "BSEQ") // Bruteforce single end with quality
    mode = BRUTE_WITH_Q;
  else {
    ERROR("Bad mode");
    usage();
    return EXIT_FAILURE;
  }

  clock_t start = clock();

  std::string config_file = options["config"];
  INFO("Loading config from " << config_file.c_str());
  if (!FileExists(config_file)) {
      ERROR("File " + config_file + " doesn't exists.");
      return EXIT_FAILURE;
  }
  cfg::create_instance(config_file);

  // hard memory limit
  // const size_t GB = 1 << 30;
  // limit_memory(cfg::get().general_hard_memory_limit * GB);

  const std::string data_base(options["database"]);
  if (!FileExists(data_base)) {
      ERROR("File " + data_base + " doesn't exists.");
      return EXIT_FAILURE;
  }
  const std::string input_file(options["input"]);
  if (!FileExists(input_file)) {
      ERROR("File " + input_file + " doesn't exists.");
      return EXIT_FAILURE;
  }
  ireadstream *input;
  INFO("Init file with reads-to-clean from " << input_file << " ... ");
  input = new ireadstream(input_file);

  std::string output_file = options["output"];
  std::ofstream output(output_file);

  INFO("Start matching reads against database ...");
  std::ofstream aligned_output(cfg::get().output_file);
  std::ofstream bed(cfg::get().output_bed);
  if (!aligned_output.is_open() || !bed.is_open() || !output.is_open()) {
    ERROR("Cannot open output file: " << cfg::get().output_file << " or "
          << cfg::get().output_bed << "or " << output_file);
    return EXIT_FAILURE;
  }

  cclean::AdapterIndex index;
  cclean::AdapterIndexBuilder().FillAdapterIndex(data_base, index);
  // Main work wrapper
  ExactAndAlign(aligned_output, bed, input, output, data_base, index, mode, options);

  aligned_output.close();
  bed.close();

  delete input;

  clock_t ends = clock();
  INFO("Processor Time Spent: " << (double) (ends - start) / CLOCKS_PER_SEC
       << " seconds.");
  INFO("Goodbye!");

  return EXIT_SUCCESS;
}
