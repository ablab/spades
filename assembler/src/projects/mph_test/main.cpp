//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dev_support/logger/log_writers.hpp"
#include "dev_support/segfault_handler.hpp"
#include "modules/data_structures/indices/perfect_hash_map.hpp"

#include "version.hpp"

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

class SimplePerfectHashMap : public debruijn_graph::KeyIteratingMap<runtime_k::RtSeq, uint32_t> {
    using base = debruijn_graph::KeyIteratingMap<runtime_k::RtSeq, uint32_t>;
  public:
    SimplePerfectHashMap(size_t k, const std::string &workdir) 
            : base(k, workdir) {}
};

class SimpleSplitter : public KMerSplitter<runtime_k::RtSeq> {
  using Seq = runtime_k::RtSeq;
    
  public:
    SimpleSplitter(const std::string &workdir, unsigned K)
            : KMerSplitter<Seq>(workdir, K) {}
    path::files_t Split(size_t num_files) override {
        path::files_t out;
        for (unsigned i = 0; i < num_files; ++i)
            out.push_back(GetRawKMersFname(i));
        std::vector<std::unique_ptr<FILE, int(*)(FILE*)> > files;
        for (const auto& file : out)
            files.emplace_back(fopen(file.c_str(), "wb"), fclose);

        for (size_t i = 0; i < 100500; ++i) {
            uint64_t data[Seq::DataSize];
            data[0] = i;
            Seq seq(128, data);
            size_t idx = GetFileNumForSeq(seq, (unsigned)num_files);
            fwrite(seq.data(), sizeof(Seq::DataType), Seq::GetDataSize(K()), files[idx].get());
        }
        
        return out;
    }
};

int main(void) {
    perf_counter pc;

    srand(42);
    srandom(42);
    try {
        create_console_logger();
        unsigned nthreads = 16;
        unsigned K = 128;
        std::string workdir = ".";

        SimplePerfectHashMap index(K, workdir);
        SimpleSplitter splitter(workdir, 128);
        KMerDiskCounter<runtime_k::RtSeq> counter(workdir, splitter);
        index.BuildIndex(counter, 1, nthreads);
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    }

    return 0;
}
