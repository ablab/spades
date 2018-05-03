//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "pipeline/graphio.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "io/reads/io_helper.hpp"

#include "version.hpp"

#include <clipp/clipp.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <string>

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

struct cfg {
    std::string load_from;
    std::string hmmfile;
    size_t k;
    bool acc;
    bool noali;
    double E; double T;
    double domE; double domT;
    double incE; double incT;
    double incdomE; double incdomT;
    bool cut_ga; bool cut_nc; bool cut_tc;
    bool max; double F1; double F2; double F3; bool nobias;

    cfg()
            : load_from(""), hmmfile(""), k(0),
              acc(false), noali(false),
              E(10.0), T(0), domE(10.0), domT(0),
              incE(0.01), incT(0.0), incdomE(0.01), incdomT(0),
              cut_ga(false), cut_nc(false), cut_tc(false),
              max(false), F1(0.02), F2(1e-3), F3(1e-5), nobias(false)
    {}
};



void process_cmdline(int argc, char **argv, cfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.hmmfile    << value("hmm file"),
      cfg.load_from  << value("load from"),
      cfg.k          << value("k-mer size"),
      // Control of output
      cfg.acc     << option("--acc")          % "prefer accessions over names in output",
      cfg.noali   << option("--noali")        % "don't output alignments, so output is smaller",
      // Control of reporting thresholds
      cfg.E       << option("-E")             % "report sequences <= this E-value threshold in output",
      cfg.T       << option("-T")             % "report sequences >= this score threshold in output",
      cfg.domE    << option("--domE")         % "report domains <= this E-value threshold in output",
      cfg.domT    << option("--domT")         % "report domains >= this score cutoff in output",
      // Control of inclusion (significance) thresholds
      cfg.incE    << option("-incE")          % "consider sequences <= this E-value threshold as significant",
      cfg.incT    << option("-incT")          % "consider sequences >= this score threshold as significant",
      cfg.incdomE << option("--incdomE")      % "consider domains <= this E-value threshold as significant",
      cfg.incdomT <<  option("--incdomT")     % "consider domains >= this score threshold as significant",
      // Model-specific thresholding for both reporting and inclusion
      cfg.cut_ga  << option("--cut_ga")       % "use profile's GA gathering cutoffs to set all thresholding",
      cfg.cut_nc  << option("--cut_nc")       % "use profile's NC noise cutoffs to set all thresholding",
      cfg.cut_tc  << option("--cut_tc")       % "use profile's TC trusted cutoffs to set all thresholding",
      // Control of acceleration pipeline
      cfg.max     << option("--max")          % "Turn all heuristic filters off (less speed, more power)",
      cfg.F1      << option("--F1")           % "Stage 1 (MSV) threshold: promote hits w/ P <= F1",
      cfg.F2      << option("--F2")           % "Stage 2 (Vit) threshold: promote hits w/ P <= F2",
      cfg.F3      << option("--F3")           % "Stage 3 (Fwd) threshold: promote hits w/ P <= F3"
  );

  if (!parse(argc, argv, cli)) {
    std::cout << make_man_page(cli, argv[0]);
    exit(1);
  }
}
int main(int argc, char* argv[]) {
    utils::perf_counter pc;

    srand(42);
    srandom(42);

    cfg cfg;
    process_cmdline(argc, argv, cfg);

    create_console_logger();
    INFO("Starting Graph HMM aligning engine, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

    using namespace debruijn_graph;
    ConjugateDeBruijnGraph graph(cfg.k);
    graphio::ScanBasicGraph(cfg.load_from, graph);

    INFO("Graph loaded. Total vertices: " << graph.size());

    for (auto it = graph.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId edge = *it;
        INFO("EdgeId: " << edge << ", length: " << graph.length(edge) << ", seq: " << graph.EdgeNucls(edge));
    }
    
    return 0;
}
