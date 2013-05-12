#include <iostream>
#include <map>
#include <exception>

#include "sequence/seq.hpp"
#include "logger/log_writers.hpp"
#include "memory_limit.hpp"
#include "QcException.hpp"
#include "running_modes.hpp"
#include "ssw/ssw_cpp.h"
#include "config_struct_cclean.hpp"
#include "simple_tools.hpp"
#include "adapter_index.hpp"
#include "valid_kmer_generator.hpp"
#include "io/read_processor.hpp"

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

#define USE_INDEX 1

class DummyIndexMatcher {
  const cclean::AdapterIndex &index_;
  size_t found_;

 public:
  DummyIndexMatcher(const cclean::AdapterIndex &index)
      : index_(index), found_(0) {}

  bool operator()(const Read &r) {
    const std::string &seq = r.getSequenceString();
    ValidKMerGenerator<cclean::K> gen(seq.c_str(), NULL, seq.size());
    while (gen.HasMore()) {
      cclean::KMer kmer = gen.kmer();
      if (index_.contains(kmer)) {
#       pragma omp atomic
        found_ += 1;
        // INFO("Contains: " << kmer << " at " << gen.pos() - 1 << " in " << r.getName());
        break;
      }
      gen.Next();
    }

    return false;
  }

  size_t found() const { return found_; }
};


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

#if USE_INDEX
    cclean::AdapterIndex index;
    cclean::AdapterIndexBuilder(16).FillAdapterIndex(index);
#endif

    if (5 != argc || (strcmp(argv[2], "exact") && strcmp(argv[2], "align") && strcmp(argv[2], "both"))) {
      usage();
      return 0;
    }

    const std::string mode(argv[2]);
    std::string db(argv[3]);
    const std::string dt(argv[4]);

    Database * data;
    ireadstream * input;

    INFO("Reading UniVec db at " << db <<  " ... ");
    data = new Database(db);
    INFO("Done");
    INFO("Init file with reads-to-clean at " << dt << " ... ");
    input = new ireadstream(dt);
    INFO("Done");

    INFO("Start matching reads against UniVec ...");
    std::ofstream output(cfg::get().output_file);
    std::ofstream bed(cfg::get().output_bed);
    if (!output.is_open() || !bed.is_open()) {
      ERROR("Cannot open output file: " << cfg::get().output_file << " or " << cfg::get().output_bed);
      return 0;
    }

#if USE_INDEX
    size_t idx = 0, found = 0;
    DummyIndexMatcher matcher(index);
    hammer::ReadProcessor rp(cfg::get().nthreads);
    rp.Run(*input, matcher);
    input->reset();
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

    INFO("Total " << rp.processed() << " reads processed. Matches found in " << matcher.found() << " reads");
#endif

    if (!mode.compare("exact")) {
      exactMatch(output, bed, input, data);
    } else if (!mode.compare("align")) {
      alignment(output, bed, input, data);
    } else if (!mode.compare("both")) {
      exactAndAlign(output, bed, input, data);
    }
    output.close();
    bed.close();

    delete data; //NB Don't delete earlier than AhoCorasick, because otherwise all patterns will be deleted
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
