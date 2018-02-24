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

#include "hmmfile.hpp"
#include <clipp/clipp.h>

#include <iostream>
#include <cstdio>

#include "demo.hpp"
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
  bool acc;
  bool show_alignments = false;
  double E; double T;
  double domE; double domT;
  double incE; double incT;
  double incdomE; double incdomT;
  bool cut_ga; bool cut_nc; bool cut_tc;
  bool max; double F1; double F2; double F3; bool nobias;

  cfg()
      : dbfile(""), hmmfile(""),
        acc(false),
        E(10.0), T(0), domE(10.0), domT(0),
        incE(0.01), incT(0.0), incdomE(0.01), incdomT(0),
        cut_ga(false), cut_nc(false), cut_tc(false),
        max(false), F1(0.02), F2(1e-3), F3(1e-5), nobias(false)
  {}
};

static int serial_master(const cfg &cfg);

static int
output_header(FILE *ofp, const std::string &hmmfile, const std::string &seqfile) {
  if (fprintf(ofp, "# query HMM file:                  %s\n", hmmfile.c_str()) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# target sequence database:        %s\n", seqfile.c_str()) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
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
      cfg.acc     << option("--acc")          % "prefer accessions over names in output",
      cfg.show_alignments   << option("--show-alignments")        % "output alignments",
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


static bool is_dir(const std::string &path) {
  struct stat sb;
  return (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode));
}



int
main(int argc, char **argv) {
  utils::segfault_handler sh;
  struct cfg cfg;
  int              status   = eslOK;

  impl_Init();          /* processor specific initialization */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */

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


P7_PIPELINE *
pipeline_Create(const cfg &cfg, int M_hint, int L_hint, int long_targets, enum p7_pipemodes_e mode) {
  P7_PIPELINE *pli  = NULL;
  int          seed = 42;

  pli = (P7_PIPELINE*)malloc(sizeof(P7_PIPELINE));

  pli->do_alignment_score_calc = 0;
  pli->long_targets = long_targets;

  if ((pli->fwd = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR_LABEL;
  if ((pli->bck = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR_LABEL;
  if ((pli->oxf = p7_omx_Create(M_hint, 0,      L_hint)) == NULL) goto ERROR_LABEL;
  if ((pli->oxb = p7_omx_Create(M_hint, 0,      L_hint)) == NULL) goto ERROR_LABEL;

  pli->r                  = esl_randomness_CreateFast(seed);
  pli->do_reseeding       = FALSE;
  pli->ddef               = p7_domaindef_Create(pli->r);
  pli->ddef->do_reseeding = pli->do_reseeding;

  /* Configure reporting thresholds */
  pli->by_E            = TRUE;
  pli->E               = cfg.E;
  pli->T               = 0.0;
  pli->dom_by_E        = TRUE;
  pli->domE            = cfg.domE;
  pli->domT            = 0.0;
  pli->use_bit_cutoffs = FALSE;
  if (cfg.T > 0.0) {
    pli->T    = cfg.T;
    pli->by_E = FALSE;
  }
  if (cfg.domT > 0.0) {
    pli->domT     = cfg.domT;
    pli->dom_by_E = FALSE;
  }

  /* Configure inclusion thresholds */
  pli->inc_by_E           = TRUE;
  pli->incE               = cfg.incE;
  pli->incT               = 0.0;
  pli->incdom_by_E        = TRUE;
  pli->incdomE            = cfg.incdomE;
  pli->incdomT            = 0.0;
  if (cfg.incT > 0.0) {
    pli->incT    = cfg.incT;
    pli->inc_by_E = FALSE;
  }
  if (cfg.incdomT > 0.0) {
    pli->incdomT     = cfg.incdomT;
    pli->incdom_by_E = FALSE;
  }

  /* Configure for one of the model-specific thresholding options */
  if (cfg.cut_ga) {
    pli->T        = pli->domT        = 0.0;
    pli->by_E     = pli->dom_by_E    = FALSE;
    pli->incT     = pli->incdomT     = 0.0;
    pli->inc_by_E = pli->incdom_by_E = FALSE;
    pli->use_bit_cutoffs = p7H_GA;
  }
  if (cfg.cut_nc) {
    pli->T        = pli->domT        = 0.0;
    pli->by_E     = pli->dom_by_E    = FALSE;
    pli->incT     = pli->incdomT     = 0.0;
    pli->inc_by_E = pli->incdom_by_E = FALSE;
    pli->use_bit_cutoffs = p7H_NC;
  }
  if (cfg.cut_tc) {
    pli->T        = pli->domT        = 0.0;
    pli->by_E     = pli->dom_by_E    = FALSE;
    pli->incT     = pli->incdomT     = 0.0;
    pli->inc_by_E = pli->incdom_by_E = FALSE;
    pli->use_bit_cutoffs = p7H_TC;
  }

  /* Configure search space sizes for E value calculations  */
  pli->Z       = pli->domZ       = 0.0;
  pli->Z_setby = pli->domZ_setby = p7_ZSETBY_NTARGETS;

  /* Configure acceleration pipeline thresholds */
  pli->do_max        = FALSE;
  pli->do_biasfilter = TRUE;
  pli->do_null2      = TRUE;
  pli->F1     = std::min(1.0, cfg.F1);
  pli->F2     = std::min(1.0, cfg.F2);
  pli->F3     = std::min(1.0, cfg.F2);

  if (long_targets) {
    pli->B1     = 100;
    pli->B2     = 240;
    pli->B3     = 1000;
  } else {
    pli->B1 = pli->B2 = pli->B3 = -1;
  }

  if (cfg.max) {
    pli->do_max        = TRUE;
    pli->do_biasfilter = FALSE;

    pli->F2 = pli->F3 = 1.0;
    pli->F1 = (pli->long_targets ? 0.3 : 1.0); // need to set some threshold for F1 even on long targets. Should this be tighter?
  }

  /* Accounting as we collect results */
  pli->nmodels         = 0;
  pli->nseqs           = 0;
  pli->nres            = 0;
  pli->nnodes          = 0;
  pli->n_past_msv      = 0;
  pli->n_past_bias     = 0;
  pli->n_past_vit      = 0;
  pli->n_past_fwd      = 0;
  pli->pos_past_msv    = 0;
  pli->pos_past_bias   = 0;
  pli->pos_past_vit    = 0;
  pli->pos_past_fwd    = 0;
  pli->mode            = mode;
  pli->show_accessions = cfg.acc;
  pli->show_alignments = cfg.show_alignments;
  pli->hfp             = NULL;
  pli->errbuf[0]       = '\0';

  return pli;

 ERROR_LABEL:
  p7_pipeline_Destroy(pli);
  return NULL;
}

class HMMMatcher {
 public:
  HMMMatcher(const hmmer::HMM &hmmw,
             const cfg &cfg)
      : gm_(NULL, p7_profile_Destroy),
        om_(NULL, p7_oprofile_Destroy),
        bg_(NULL, p7_bg_Destroy),
        pli_(NULL, p7_pipeline_Destroy),
        th_(NULL, p7_tophits_Destroy) {
    P7_HMM *hmm = hmmw.get();

    P7_PROFILE *gm = p7_profile_Create (hmm->M, hmm->abc);
    gm_.reset(gm);

    P7_OPROFILE *om = p7_oprofile_Create(hmm->M, hmm->abc);
    om_.reset(om);

    P7_BG *bg = p7_bg_Create(hmm->abc);
    bg_.reset(bg);

    p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
    p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */

    P7_PIPELINE *pli = pipeline_Create(cfg, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
    pli_.reset(pli);

    p7_pli_NewModel(pli, om, bg);

    th_.reset(p7_tophits_Create());
  }

  void match(const char *name, const char *seq, const char *desc = NULL) {
    ESL_SQ *dbsq = esl_sq_CreateFrom(name, seq, desc, NULL, NULL);
    esl_sq_Digitize(om_->abc, dbsq);

    p7_pli_NewSeq(pli_.get(), dbsq);
    p7_bg_SetLength(bg_.get(), dbsq->n);
    p7_oprofile_ReconfigLength(om_.get(), dbsq->n);

    p7_Pipeline(pli_.get(), om_.get(), bg_.get(), dbsq, nullptr, th_.get());
    p7_pipeline_Reuse(pli_.get());

    esl_sq_Destroy(dbsq);
  }

  void summarize() {
    p7_tophits_SortBySortkey(th_.get());
    p7_tophits_Threshold(th_.get(), pli_.get());
  }

  P7_TOPHITS *hits() const {
    return th_.get();
  }
  P7_PIPELINE *pipeline() const {
    return pli_.get();
  }

 private:
  std::unique_ptr<P7_PROFILE, void(*)(P7_PROFILE*)>  gm_;
  std::unique_ptr<P7_OPROFILE, void(*)(P7_OPROFILE*)> om_;  /* optimized query profile */
  std::unique_ptr<P7_BG, void(*)(P7_BG*)> bg_; /* null model */
  std::unique_ptr<P7_PIPELINE, void(*)(P7_PIPELINE*)> pli_; /* work pipeline */
  std::unique_ptr<P7_TOPHITS, void(*)(P7_TOPHITS*)> th_;
};

static int
serial_master(const cfg &cfg) {
  FILE            *ofp      = stdout;
  ESL_STOPWATCH   *w;
  int              textw    = 120;

  w = esl_stopwatch_Create();

  /* Open the query profile HMM file */
  hmmer::HMMFile hmmfile(cfg.hmmfile);
  if (!hmmfile.valid())
    p7_Fail("Error reading HMM file %s\n", cfg.hmmfile.c_str());

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

    std::unordered_set<DBGraph::GraphPointer> begins, ends;
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
    HMMMatcher matcher(hmmw.get(), cfg);

    if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    if (hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
    if (hmm->desc) { if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

    for (size_t i = 0; i < edges.size(); ++i) {
      matcher.match((std::to_string(i)).c_str(), edges[i].c_str(), "comment");
      if (i % 10000 == 0) {
        TRACE(i << " edges hmmered");
      }
    }
    matcher.summarize();

    esl_stopwatch_Stop(w);

    // p7_tophits_Targets(ofp, matcher.hits(), matcher.pipeline(), textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    // p7_tophits_Domains(ofp, matcher.hits(), matcher.pipeline(), textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    // p7_pli_Statistics(ofp, matcher.pipeline(), w); if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

    auto hits = matcher.hits();
    std::vector<size_t> good_edges;
    size_t N = hits->N;
    if (cfg.top_hmmer_matches) {
      N = std::min(cfg.top_hmmer_matches, N);
    }
    for (size_t i = 0; i < N; ++i) {
      size_t id = std::stoll(hits->hit[i]->name);
      good_edges.push_back(id);
    }

    auto subgraph = graph.subgraph(good_edges, good_edges, hmm->M * 2);

    INFO("Edges: " << subgraph.n_edges() << ", vertices (k-mers): " << subgraph.n_vertices());
    if (subgraph.loops()) {
      INFO("Subgraph has loops");
    } else {
      INFO("Subgraph has NO loops");
    }

    auto fees = fees_from_hmm(hmm, hmmw->abc());

    auto general_subgraph = subgraph.general_graph();
    INFO("Before compression");
    INFO("Edges: " << general_subgraph.n_edges() << ", vertices (k-mers): " << general_subgraph.n_vertices() << ", bases: " << general_subgraph.n_bases());
    general_subgraph.check_symmetry();
    general_subgraph.collapse();
    INFO("After compression");
    INFO("Edges: " << general_subgraph.n_edges() << ", vertices (k-mers): " << general_subgraph.n_vertices() << ", bases: " << general_subgraph.n_bases());
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

      std::vector<ReversalGraphPointer<Graph::GraphPointer>> rev_initial(CONST_ALL(initial));
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
