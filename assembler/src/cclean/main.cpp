#include <iostream>
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

void usage() {
  std::cout << "This tool searches contaminations from UniVec db in provided file with reads" << std::endl;
  std::cout << "Usage: QC-pileline config_path mode:{exact, align, both} UniVec_path Fasta/Fastq.gz" << std::endl;
  std::cout << "Currently only .gz files can be read" << std::endl;
}

void create_console_logger() {
  using namespace logging;

  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

int main(int argc, char *argv[]) {
  try {
    create_console_logger();

    clock_t start = clock();

    std::string config_file = "config.inp";
    if (argc > 1) config_file = argv[1];
    INFO("Loading config from " << config_file.c_str());
    CheckFileExistenceFATAL(config_file);
    cfg::create_instance(config_file);

    // hard memory limit
    // const size_t GB = 1 << 30;
    // limit_memory(cfg::get().general_hard_memory_limit * GB);

    if (4 != argc) {
      usage();
      return 0;
    }

    std::string db(argv[2]);
    const std::string dt(argv[3]);

    cclean::AdapterIndex index;
    cclean::AdapterIndexBuilder().FillAdapterIndex(db, index);

    ireadstream * input;

    INFO("Init file with reads-to-clean from " << dt << " ... ");
    input = new ireadstream(dt);

    INFO("Start matching reads against database ...");
    std::ofstream output(cfg::get().output_file);
    std::ofstream bed(cfg::get().output_bed);
    if (!output.is_open() || !bed.is_open()) {
      ERROR("Cannot open output file: " << cfg::get().output_file << " or " << cfg::get().output_bed);
      return 0;
    }

    exactAndAlign(output, bed, input, index);

    output.close();
    bed.close();

    delete input;

    clock_t ends = clock();
    INFO("Processor Time Spent: " << (double) (ends - start) / CLOCKS_PER_SEC << " seconds.");
    INFO("Goodbye!");
  } catch (std::exception& e) {
    ERROR("Error: " << e.what());
    return 0;
  }

  return 0;
}
