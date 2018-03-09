#include "utils.hpp"

extern "C" {
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
}

#include <clipp/clipp.h>
#include "hmmfile.hpp"
#include "hmmmatcher.hpp"

#include <cstdio>
#include <iostream>

#include "demo.hpp"
#include "cursor.hpp"
#include "graph.hpp"
#include "fees.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"

#include <sys/stat.h>
#include <fstream>

void create_console_logger() {
  using namespace logging;

  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

struct cfg {
  std::string dbfile;
  std::string hmmfile;
  std::string genefile = "";
  std::string log = "";
  std::string output_dir = "";
  bool debug = false;
  size_t top_hmmer_matches = 0;
  size_t k;
  hmmer::hmmer_cfg hcfg;

  cfg()
      : dbfile(""), hmmfile("") {}
};

static int serial_master(const cfg &cfg);

static int output_header(FILE *ofp, const std::string &hmmfile, const std::string &seqfile) {
  if (fprintf(ofp, "# query HMM file:                  %s\n", hmmfile.c_str()) < 0)
    ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# target sequence database:        %s\n", seqfile.c_str()) < 0)
    ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n") < 0)
    ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}

void process_cmdline(int argc, char **argv, cfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.hmmfile << value("hmm file"),
      cfg.dbfile  << value("graph edges file"),
      required("-k") & integer("k", cfg.k)    % "edge length in dB graph",
      required("--output", "-o") & value("output directory", cfg.output_dir)    % "output directory",
      option("--top-hmmer") & integer("N", cfg.top_hmmer_matches)    % "use only N hmmer hmmsearch top matches",
      option("--log") & value("log file", cfg.log)    % "log file",
      option("--genes") & value("gene sequence file", cfg.genefile)    % "gene sequence file",
      cfg.debug   << option("--debug")        % "perform additional debug operations (i.e. found path checking, search in reverted graph, etc)",

      // Control of output
      cfg.hcfg.acc     << option("--acc")          % "prefer accessions over names in output",
      cfg.hcfg.noali   << option("--noali")        % "don't output alignments, so output is smaller",
      // Control of reporting thresholds
      (option("-E") & number("value", cfg.hcfg.E))        % "report sequences <= this E-value threshold in output",
      (option("-T") & number("value", cfg.hcfg.T))        % "report sequences >= this score threshold in output",
      (option("--domE") & number("value", cfg.hcfg.domE)) % "report domains <= this E-value threshold in output",
      (option("--domT") & number("value", cfg.hcfg.domT)) % "report domains >= this score cutoff in output",
      // Control of inclusion (significance) thresholds
      (option("-incE") & number("value", cfg.hcfg.incE))       % "consider sequences <= this E-value threshold as significant",
      (option("-incT") & number("value", cfg.hcfg.incT))       % "consider sequences >= this score threshold as significant",
      (option("-incdomE") & number("value", cfg.hcfg.incdomE)) % "consider domains <= this E-value threshold as significant",
      (option("-incdomT") & number("value", cfg.hcfg.incdomT)) % "consider domains >= this score threshold as significant",
      // Model-specific thresholding for both reporting and inclusion
      cfg.hcfg.cut_ga  << option("--cut_ga")       % "use profile's GA gathering cutoffs to set all thresholding",
      cfg.hcfg.cut_nc  << option("--cut_nc")       % "use profile's NC noise cutoffs to set all thresholding",
      cfg.hcfg.cut_tc  << option("--cut_tc")       % "use profile's TC trusted cutoffs to set all thresholding",
      // Control of acceleration pipeline
      cfg.hcfg.max     << option("--max")             % "Turn all heuristic filters off (less speed, more power)",
      (option("--F1") & number("value", cfg.hcfg.F1)) % "Stage 1 (MSV) threshold: promote hits w/ P <= F1",
      (option("--F2") & number("value", cfg.hcfg.F2)) % "Stage 2 (Vit) threshold: promote hits w/ P <= F2",
      (option("--F3") & number("value", cfg.hcfg.F3)) % "Stage 3 (Fwd) threshold: promote hits w/ P <= F3"
  );

  if (!parse(argc, argv, cli)) {
    std::cout << make_man_page(cli, argv[0]);
    exit(1);
  }
}

static bool is_dir(const std::string &path) {
  struct stat sb;
  return (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode));
}

int main(int argc, char **argv) {
  utils::segfault_handler sh;
  struct cfg cfg;
  int status = eslOK;

  impl_Init();      /* processor specific initialization */
  p7_FLogsumInit(); /* we're going to use table-driven Logsum() approximations at times */

  process_cmdline(argc, argv, cfg);
  if (is_dir(cfg.output_dir)) {
    WARN("Output directory " << cfg.output_dir << " exists");
  } else {
    int status = mkdir(cfg.output_dir.c_str(), 0775);
    if (status != 0) {
      ERROR("Cannot create output directory: " << cfg.output_dir);
      std::exit(1);
    }
  }
  if (cfg.log == "") {
    cfg.log = cfg.output_dir + "/log.log";
  }
  if (cfg.log != "") {
    create_console_logger();
  }
  INFO("Looper started");

  status = serial_master(cfg);

  return status;
}

static int serial_master(const cfg &cfg) {
  FILE *ofp = stdout;
  ESL_STOPWATCH *w = esl_stopwatch_Create();

  /* Open the query profile HMM file */
  hmmer::HMMFile hmmfile(cfg.hmmfile);
  if (!hmmfile.valid()) p7_Fail("Error reading HMM file %s\n", cfg.hmmfile.c_str());

  auto hmmw = hmmfile.read();
  if (hmmw) {
    /* One-time initializations after alphabet <abc> becomes known */
    output_header(ofp, cfg.hmmfile, cfg.dbfile);
  }

  // Read all edges (along with RC)
  auto edges = read_fasta_edges(cfg.dbfile, true);
  INFO("#Edges readed (with RC): " << edges.size());
  remove_duplicates(edges);
  INFO("#Edges remaining: " << edges.size());
  DBGraph graph(cfg.k, edges);
  if (cfg.genefile != "") {
    INFO("Tracing gene paths...");
    auto genes = read_fasta_edges(cfg.genefile, false);
    std::vector<DBGraph::Path> paths;
    INFO(genes.size() << " genes loaded...");
    for (const auto gene : genes) {
      auto path = graph.trace_exact_sequence(gene);
      INFO(path);
      paths.push_back(path);
    }

    std::sort(paths.begin(), paths.end(),
              [](const auto &e1, const auto &e2) { return e1.turns.length() < e2.turns.length(); });
    INFO(paths);

    std::unordered_set<DBGraph::GraphCursor> begins, ends;
    for (const auto &path : paths) {
      begins.insert(path.begin);
      ends.insert(path.end);
    }

    INFO("BEGINS");
    INFO(begins);
    INFO("ENDS");
    INFO(ends);
  }

  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hmmw) {
    P7_HMM *hmm = hmmw->get();

    esl_stopwatch_Start(w);
    hmmer::HMMMatcher matcher(hmmw.get(), cfg.hcfg);

    if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    if (hmm->acc) {
      if (fprintf(ofp, "Accession:   %s\n", hmm->acc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    }
    if (hmm->desc) {
      if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    }

    for (size_t i = 0; i < edges.size(); ++i) {
      matcher.match((std::to_string(i)).c_str(), edges[i].c_str(), "comment");
      if (i % 10000 == 0) {
        TRACE(i << " edges hmmered");
      }
    }
    matcher.summarize();

    esl_stopwatch_Stop(w);

    auto hits = matcher.hits();
    std::vector<size_t> good_edges;
    size_t N = hits->N;
    if (cfg.top_hmmer_matches)
      N = std::min(cfg.top_hmmer_matches, N);
    for (size_t i = 0; i < N; ++i) {
      if (!(hits->hit[i]->flags & p7_IS_REPORTED))
        continue;
      if (!(hits->hit[i]->flags & p7_IS_INCLUDED))
        continue;

      size_t id = std::stoll(hits->hit[i]->name);
      good_edges.push_back(id);
    }
    INFO("Total matched edges: " << good_edges.size());

    auto subgraph = graph.subgraph(good_edges, good_edges, hmm->M * 2);

    INFO("Edges: " << subgraph.n_edges() << ", vertices (k-mers): " << subgraph.n_vertices());
    if (subgraph.loops()) {
      INFO("Subgraph has loops");
    } else {
      INFO("Subgraph has NO loops");
    }

    auto fees = hmm::fees_from_hmm(hmm, hmmw->abc());

    auto general_subgraph = subgraph.general_graph();
    INFO("Before compression");
    INFO("Edges: " << general_subgraph.n_edges() << ", vertices (k-mers): " << general_subgraph.n_vertices()
                   << ", bases: " << general_subgraph.n_bases());
    general_subgraph.check_symmetry();
    general_subgraph.collapse();
    INFO("After compression");
    INFO("Edges: " << general_subgraph.n_edges() << ", vertices (k-mers): " << general_subgraph.n_vertices()
                   << ", bases: " << general_subgraph.n_bases());
    general_subgraph.check_symmetry();
    auto initial = general_subgraph.all();
    auto result = find_best_path(fees, initial);
    INFO(result);

    std::ofstream ofs(cfg.output_dir + "/" + hmm->name + ".fa", std::ios::out);
    ofs.precision(13);
    for (const auto &kv : result) {
      ofs << ">Score=" << kv.second << "\n";
      ofs << kv.first << "\n";
    }

    if (cfg.debug) {
      INFO("Checking path presence in the initial graph");
      for (const auto gene_score : result) {
        auto path = graph.trace_exact_sequence(gene_score.first);
        TRACE(path);
      }

      std::vector<ReversalGraphCursor<Graph::GraphCursor>> rev_initial(CONST_ALL(initial));
      auto reversed_result = find_best_path_rev(fees.reversed(), rev_initial);
      INFO(reversed_result);

      for (const auto &kv : reversed_result) {
        ofs << ">Score=" << kv.second << "\n";
        std::string seq = kv.first;
        std::reverse(ALL(seq));
        ofs << seq << "\n";
      }
    }

    if (cfg.genefile != "") {
      auto genes = read_fasta_edges(cfg.genefile, false);
      Graph geneg(genes);
      auto initial = geneg.begins();
      auto result = find_best_path(fees, initial);

      std::ofstream gene_ofs(cfg.output_dir + "/" + hmm->name + "_genes.fa", std::ios::out);
      gene_ofs.precision(13);
      for (const auto &kv : result) {
        gene_ofs << ">Score=" << kv.second << "\n";
        gene_ofs << kv.first << "\n";
      }
    }

    hmmw = hmmfile.read();
  } /* end outer loop over query HMMs */

#if 0
  switch(hstatus) {
  case eslEOD:       p7_Fail("read failed, HMM file %s may be truncated?", cfg.hmmfile.c_str());      break;
  case eslEFORMAT:   p7_Fail("bad file format in HMM file %s",             cfg.hmmfile.c_str());      break;
  case eslEINCOMPAT: p7_Fail("HMM file %s contains different alphabets",   cfg.hmmfile.c_str());      break;
  case eslEOF:       /* do nothing. EOF is what we want. */                                    break;
  default:           p7_Fail("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg.hmmfile.c_str());
  }
#endif

  /* Cleanup - prepare for exit */
  esl_stopwatch_Destroy(w);

  return eslOK;
}
