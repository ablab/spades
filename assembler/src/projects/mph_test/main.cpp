//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/indices/perfect_hash_map.hpp"
#include "utils/mph_index/kmer_index_builder.hpp"

#include "io/reads/read_processor.hpp"
#include "io/reads/io_helper.hpp"

#include "version.hpp"

#include <cxxopts/cxxopts.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <string>

using namespace std;
void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

class SimplePerfectHashMap : public debruijn_graph::KeyIteratingMap<RtSeq, uint32_t> {
    using base = debruijn_graph::KeyIteratingMap<RtSeq, uint32_t>;
  public:
    SimplePerfectHashMap(size_t k, const std::string &workdir)
            : base(k, workdir) {}
};

class ParallelSortingSplitter : public KMerSortingSplitter<RtSeq> {
  using Seq = RtSeq;

  std::vector<std::string> files_;
  unsigned nthreads_;
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
    ParallelSortingSplitter(const std::string &workdir, unsigned K, unsigned nthreads, size_t read_buffer_size = 0)
            : KMerSortingSplitter<Seq>(workdir, K), nthreads_(nthreads), read_buffer_size_(read_buffer_size) {}

    void push_back(const std::string &filename) {
        files_.push_back(filename);
    }

    path::files_t Split(size_t num_files) override {
        INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

        path::files_t out = PrepareBuffers(num_files, nthreads_, read_buffer_size_);

        size_t n = 10;
        BufferFiller filler(*this, K());
        for (const auto &file : files_) {
            INFO("Processing " << file);
            auto irs = io::EasyStream(file, true, true);
            while (!irs->eof()) {
                hammer::ReadProcessor rp(nthreads_);
                rp.Run(*irs, filler);
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

int main(int argc, char* argv[]) {
    perf_counter pc;

    srand(42);
    srandom(42);
    try {
        unsigned nthreads;
        unsigned K;
        std::string workdir, dataset;
        std::vector<std::string> input;
        size_t read_buffer_size;

        cxxopts::Options options(argv[0], " <input files> - SPAdes k-mer counting engine");
        options.add_options()
                ("k,kmer", "K-mer length", cxxopts::value<unsigned>(K)->default_value("21"), "K")
                ("d,dataset", "Dataset description (in YAML), input files ignored", cxxopts::value<std::string>(dataset), "file")
                ("t,threads", "# of threads to use", cxxopts::value<unsigned>(nthreads)->default_value(std::to_string(omp_get_max_threads())), "num")
                ("w,workdir", "Working directory to use", cxxopts::value<std::string>(workdir)->default_value("."), "dir")
                ("b,bufsize", "Sorting buffer size, per thread", cxxopts::value<size_t>(read_buffer_size)->default_value("536870912"))
                ("h,help", "Print help");

        options.add_options("Input")
                ("positional", "", cxxopts::value<std::vector<std::string>>(input));

        options.parse_positional("positional");
        options.parse(argc, argv);
        if (options.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        if (!options.count("positional") && !options.count("dataset")) {
            std::cerr << "ERROR: No input files were specified" << std::endl << std::endl;
            std::cout << options.help() << std::endl;
            exit(-1);
        }

        create_console_logger();

        INFO("Starting SPAdes k-mer counting engine, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

        INFO("K-mer length set to " << K);
        INFO("# of threads to use: " << nthreads);

        SimplePerfectHashMap index(K, workdir);
        ParallelSortingSplitter splitter(workdir, K, nthreads, read_buffer_size);
        if (options.count("dataset")) {
            io::DataSet<> idataset;
            idataset.load(dataset);
            for (const auto &s : idataset.reads())
                splitter.push_back(s);
        } else {
            for (const auto& s : input)
                splitter.push_back(s);
        }
        KMerDiskCounter<RtSeq> counter(workdir, splitter);
        counter.CountAll(16, nthreads);
        INFO("K-mer counting done, kmers saved to " << counter.GetFinalKMersFname());
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    } catch (const cxxopts::OptionException &e) {
        std::cerr << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }

    return 0;
}
