//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "version.hpp"

#include "pipeline/library.hpp"

#include "io/reads/read_processor.hpp"
#include "io/reads/io_helper.hpp"

#include "utils/ph_map/kmer_maps.hpp"
#include "utils/kmer_mph/kmer_index_builder.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"

#include <clipp/clipp.h>
#include <string>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cerrno>

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;
void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

class SimplePerfectHashMap : public utils::KeyIteratingMap<RtSeq, uint32_t> {
    using base = utils::KeyIteratingMap<RtSeq, uint32_t>;
  public:
    SimplePerfectHashMap(unsigned k)
            : base(k) {}
};

class ParallelSortingSplitter : public kmers::KMerSortingSplitter<RtSeq> {
  using Seq = RtSeq;

  std::vector<std::string> files_;
  size_t read_buffer_size_;

  class BufferFiller {
      size_t processed_;
      ParallelSortingSplitter &splitter_;
      unsigned K_;

    public:
      BufferFiller(ParallelSortingSplitter &splitter, unsigned K)
              : processed_(0), splitter_(splitter), K_(K) {}

      size_t processed() const { return processed_; }

      bool operator()(std::unique_ptr<io::SingleRead> r) {
#         pragma omp atomic
          processed_ += 1;

          const Sequence &seq = r->sequence();

          if (seq.size() < this->K_)
              return false;

          unsigned thread_id = omp_get_thread_num();
          bool stop = false;
          RtSeq kmer = seq.start<RtSeq>(this->K_) >> 'A';
          for (size_t j = this->K_ - 1; j < seq.size(); ++j) {
              kmer <<= seq[j];
              stop |= splitter_.push_back_internal(kmer, thread_id);
          }

          return stop;
      }
  };


  public:
    using kmers::KMerSortingSplitter<RtSeq>::RawKMers;
    ParallelSortingSplitter(const std::string &workdir, unsigned K, size_t read_buffer_size = 0)
            : KMerSortingSplitter<Seq>(workdir, K), read_buffer_size_(read_buffer_size) {}

    void push_back(const std::string &filename) {
        files_.push_back(filename);
    }

    RawKMers Split(size_t num_files, unsigned nthreads) override {
        auto out = PrepareBuffers(num_files, nthreads, read_buffer_size_);

        size_t n = 10;
        BufferFiller filler(*this, K());
        for (const auto &file : files_) {
            INFO("Processing " << file);
            auto irs = io::EasyStream(file, true, true);
            while (!irs.eof()) {
                hammer::ReadProcessor rp(nthreads);
                rp.Run(irs, filler);
                DumpBuffers(out);
                VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

                if (filler.processed() >> n) {
                    INFO("Processed " << filler.processed() << " reads");
                    n += 1;
                }
            }
        }
        INFO("Total " << filler.processed() << " reads processed");

        this->ClearBuffers();

        return out;
    }
};

namespace kmer_count {
struct Args {
    unsigned nthreads = omp_get_max_threads();
    unsigned K = 21;
    std::string workdir, dataset = "";
    size_t read_buffer_size = 536870912;
    std::vector<std::string> input;
};
}

void process_cmdline(int argc, char **argv, kmer_count::Args &args) {
    using namespace clipp;
    bool print_help = false;

    auto cli = (
        (option("-k", "--kmer") & integer("value", args.K)) % "K-mer length",
        (option("-d", "--dataset") & value("dir", args.dataset)) % "Dataset description (in YAML), input files ignored",
        (option("-t", "--threads") & integer("value", args.nthreads)) % "# of threads to use",
        (option("-w", "--workdir") & value("dir", args.workdir)) % "Working directory to use",
        (option("-b", "--bufsize") & integer("value", args.read_buffer_size)) % "Sorting buffer size, per thread",
        (option("-h", "--help").set(print_help)) % "Show help",
        opt_values("input files", args.input)
    );

    auto help_message = make_man_page(cli, argv[0])
        .prepend_section("DESCRIPTION",
                         "SPAdes k-mer counting engine\n\n"
                             "Output: <output_dir>/final_kmers - unordered set of kmers in binary format. "
                             "Kmers from both forward and reverse-complementary reads are taken into account.\n\n"
                             "Output format: All kmers are written sequentially without any separators. "
                             "Each kmer takes the same number of bits. One kmer of length K takes 2*K bits. "
                             "Kmers are aligned by 64 bits. "
                             "For example, one kmer with length=21 takes 8 bytes, with length=33 takes 16 bytes, "
                             "and with length=55 takes 16 bytes. "
                             "Each nucleotide is coded with 2 bits: 00 - A, 01 - C, 10 - G, 11 - T.\n\n"
                             "Example: For kmer: AGCTCT\n"
                             "\tMemory: 6 bits * 2 = 12, 64 bits (8 bytes)\n"
                             "\tLet’s describe bytes: \n"
                             "\tdata[0] = AGCT -> 11 01 10 00 -> 0xd8\n"
                             "\tdata[1] = CT00 -> 00 00 11 01 -> 0x0d\n"
                             "\tdata[2] = 0000 -> 00 00 00 00 -> 0x00\n"
                             "\tdata[3] = 0000 -> 00 00 00 00 -> 0x00\n"
                             "\tdata[4] = 0000 -> 00 00 00 00 -> 0x00\n"
                             "\tdata[5] = 0000 -> 00 00 00 00 -> 0x00\n"
                             "\tdata[6] = 0000 -> 00 00 00 00 -> 0x00\n"
                             "\tdata[7] = 0000 -> 00 00 00 00 -> 0x00\n");
    auto result = parse(argc, argv, cli);
    if (!result || print_help) {
        std::cout << help_message;
        if (print_help) {
            exit(0);
        } else {
            exit(1);
        }
    }

    if (args.input.size() == 0 && args.dataset == "") {
        std::cerr << "ERROR: No input files were specified" << std::endl << std::endl;
        std::cout << help_message << std::endl;
        exit(-1);
    }
}

int main(int argc, char* argv[]) {
    utils::perf_counter pc;

    srand(42);
    srandom(42);
    try {
        kmer_count::Args args;
        process_cmdline(argc, argv, args);

        create_console_logger();

        START_BANNER("SPAdes k-mer counting engine");

        INFO("K-mer length set to " << args.K);
        INFO("# of threads to use: " << args.nthreads);

        SimplePerfectHashMap index(args.K);
        ParallelSortingSplitter splitter(args.workdir, args.K, args.read_buffer_size);
        if (args.dataset != "") {
            io::DataSet<> idataset;
            idataset.load(args.dataset);
            for (const auto &s : idataset.reads())
                splitter.push_back(s);
        } else {
            for (const auto& s : args.input)
                splitter.push_back(s);
        }

        kmers::KMerDiskCounter<RtSeq> counter(args.workdir, std::move(splitter));
        auto res = counter.CountAll(16, args.nthreads, /* merge */ true);
        auto final_kmers = res.final_kmers();

        std::string outputfile_name = fs::append_path(args.workdir, "final_kmers");
        std::rename(final_kmers->file().c_str(), outputfile_name.c_str());

        INFO("K-mer counting done, kmers saved to " << outputfile_name);
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    }

    return 0;
}
