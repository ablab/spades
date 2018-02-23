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

#include "hmmfile.hpp"

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

extern "C" {
#include "p7_config.h"

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

  P7_PIPELINE *
  pipeline_Create(const cfg &cfg, int M_hint, int L_hint, int long_targets, enum p7_pipemodes_e mode) {
    P7_PIPELINE *pli  = NULL;
    int          seed = 42;

    pli = (P7_PIPELINE*)malloc(sizeof(P7_PIPELINE));

    pli->do_alignment_score_calc = 0;
    pli->long_targets = long_targets;

    if ((pli->fwd = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR;
    if ((pli->bck = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR;
    if ((pli->oxf = p7_omx_Create(M_hint, 0,      L_hint)) == NULL) goto ERROR;
    if ((pli->oxb = p7_omx_Create(M_hint, 0,      L_hint)) == NULL) goto ERROR;

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
    pli->show_alignments = !cfg.noali;
    pli->hfp             = NULL;
    pli->errbuf[0]       = '\0';

    return pli;

 ERROR:
    p7_pipeline_Destroy(pli);
    return NULL;
  }

  std::unique_ptr<P7_PROFILE, void(*)(P7_PROFILE*)>  gm_;
  std::unique_ptr<P7_OPROFILE, void(*)(P7_OPROFILE*)> om_;  /* optimized query profile */
  std::unique_ptr<P7_BG, void(*)(P7_BG*)> bg_; /* null model */
  std::unique_ptr<P7_PIPELINE, void(*)(P7_PIPELINE*)> pli_; /* work pipeline */
  std::unique_ptr<P7_TOPHITS, void(*)(P7_TOPHITS*)> th_;
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
      (option("-E") & number("value", cfg.E))        % "report sequences <= this E-value threshold in output",
      (option("-T") & number("value", cfg.E))        % "report sequences >= this score threshold in output",
      (option("--domE") & number("value", cfg.domE)) % "report domains <= this E-value threshold in output",
      (option("--domT") & number("value", cfg.domT)) % "report domains >= this score cutoff in output",
      // Control of inclusion (significance) thresholds
      (option("-incE") & number("value", cfg.incE))       % "consider sequences <= this E-value threshold as significant",
      (option("-incT") & number("value", cfg.incT))       % "consider sequences >= this score threshold as significant",
      (option("-incdomE") & number("value", cfg.incdomE)) % "consider domains <= this E-value threshold as significant",
      (option("-incdomT") & number("value", cfg.incdomT)) % "consider domains >= this score threshold as significant",
      // Model-specific thresholding for both reporting and inclusion
      cfg.cut_ga  << option("--cut_ga")       % "use profile's GA gathering cutoffs to set all thresholding",
      cfg.cut_nc  << option("--cut_nc")       % "use profile's NC noise cutoffs to set all thresholding",
      cfg.cut_tc  << option("--cut_tc")       % "use profile's TC trusted cutoffs to set all thresholding",
      // Control of acceleration pipeline
      cfg.max     << option("--max")             % "Turn all heuristic filters off (less speed, more power)",
      (option("--F1") & number("value", cfg.F1)) % "Stage 1 (MSV) threshold: promote hits w/ P <= F1",
      (option("--F2") & number("value", cfg.F2)) % "Stage 2 (Vit) threshold: promote hits w/ P <= F2",
      (option("--F3") & number("value", cfg.F3)) % "Stage 3 (Fwd) threshold: promote hits w/ P <= F3"
  );

  if (!parse(argc, argv, cli)) {
    std::cout << make_man_page(cli, argv[0]);
    exit(1);
  }
}

int main(int argc, char* argv[]) {
    utils::perf_counter pc;
    int textw = 120;

    srand(42);
    srandom(42);

    cfg cfg;
    process_cmdline(argc, argv, cfg);

    create_console_logger();
    INFO("Starting Graph HMM aligning engine, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

    /* Open the query profile HMM file */
    hmmer::HMMFile hmmfile(cfg.hmmfile);
    if (!hmmfile.valid())
        FATAL_ERROR("Error reading HMM file "<< cfg.hmmfile);

    using namespace debruijn_graph;
    ConjugateDeBruijnGraph graph(cfg.k);
    graphio::ScanBasicGraph(cfg.load_from, graph);
    INFO("Graph loaded. Total vertices: " << graph.size());

    auto hmmw = hmmfile.read();
    ESL_STOPWATCH *w = esl_stopwatch_Create();
    /* Outer loop: over each query HMM in <hmmfile>. */
    while (hmmw) {
        P7_HMM *hmm = hmmw->get();

        HMMMatcher matcher(hmmw.get(), cfg);

        if (fprintf(stderr, "Query:       %s  [M=%d]\n", hmm->name, hmm->M) < 0) FATAL_ERROR("write failed");
        if (hmm->acc)  { if (fprintf(stderr, "Accession:   %s\n", hmm->acc)  < 0) FATAL_ERROR("write failed"); }
        if (hmm->desc) { if (fprintf(stderr, "Description: %s\n", hmm->desc) < 0) FATAL_ERROR("write failed"); }

        esl_stopwatch_Start(w);
        for (auto it = graph.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            EdgeId edge = *it;
            // FIXME: this conversion is pointless
            std::string ref = std::to_string(graph.int_id(edge));
            std::string seq = graph.EdgeNucls(edge).str();
            //  INFO("EdgeId: " << edge << ", length: " << graph.length(edge) << ", seq: " << graph.EdgeNucls(edge));
            matcher.match(ref.c_str(), seq.c_str());
        }

        matcher.summarize();
        esl_stopwatch_Stop(w);

        p7_tophits_Targets(stderr, matcher.hits(), matcher.pipeline(), textw); if (fprintf(stderr, "\n\n") < 0) FATAL_ERROR("write failed");
        p7_tophits_Domains(stderr, matcher.hits(), matcher.pipeline(), textw); if (fprintf(stderr, "\n\n") < 0) FATAL_ERROR("write failed");
        p7_pli_Statistics(stderr, matcher.pipeline(), w); if (fprintf(stderr, "//\n") < 0) FATAL_ERROR("write failed");

        hmmw = hmmfile.read();
    } /* end outer loop over query HMMs */

    esl_stopwatch_Destroy(w);

    return 0;
}
