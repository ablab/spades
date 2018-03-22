#include "hmmmatcher.hpp"

extern "C" {
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
}

#include "hmmfile.hpp"

using namespace hmmer;

HMMMatcher::HMMMatcher(const HMM &hmmw,
                       const hmmer_cfg &cfg)
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

    P7_PIPELINE *pli = pipeline_create(cfg, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
    pli_.reset(pli);

    p7_pli_NewModel(pli, om, bg);

    th_.reset(p7_tophits_Create());
}

void HMMMatcher::match(const char *name, const char *seq, const char *desc) {
    ESL_SQ *dbsq = esl_sq_CreateFrom(name, seq, desc, NULL, NULL);
    esl_sq_Digitize(om_->abc, dbsq);

    p7_pli_NewSeq(pli_.get(), dbsq);
    p7_bg_SetLength(bg_.get(), int(dbsq->n));
    p7_oprofile_ReconfigLength(om_.get(), int(dbsq->n));

    p7_Pipeline(pli_.get(), om_.get(), bg_.get(), dbsq, nullptr, th_.get());
    p7_pipeline_Reuse(pli_.get());

    esl_sq_Destroy(dbsq);
}

void HMMMatcher::summarize() {
    p7_tophits_SortBySortkey(th_.get());
    p7_tophits_Threshold(th_.get(), pli_.get());
}

P7_TOPHITS *HMMMatcher::top_hits() const {
    return th_.get();
}
P7_PIPELINE *HMMMatcher::pipeline() const {
    return pli_.get();
}

P7_PIPELINE *
HMMMatcher::pipeline_create(const hmmer_cfg &cfg, int M_hint, int L_hint, int long_targets, unsigned mode) {
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
    pli->mode            = (p7_pipemodes_e)mode;
    pli->show_accessions = cfg.acc;
    pli->show_alignments = !cfg.noali;
    pli->hfp             = NULL;
    pli->errbuf[0]       = '\0';

    return pli;

ERROR:
    p7_pipeline_Destroy(pli);
    return NULL;
}

const char *HMMMatcher::Hit::name() const { return h_->name; }
const char *HMMMatcher::Hit::acc() const { return h_->acc; }
const char *HMMMatcher::Hit::desc() const { return h_->desc; }
float HMMMatcher::Hit::score() const { return h_->score; }
double HMMMatcher::Hit::lnP() const { return h_->lnP; }
size_t HMMMatcher::Hit::ndom() const { return h_->ndom; }
uint32_t HMMMatcher::Hit::flags() const { return h_->flags; }
bool HMMMatcher::Hit::reported() const { return flags() & p7_IS_REPORTED; }
bool HMMMatcher::Hit::included() const { return flags() & p7_IS_INCLUDED; }

const HMMMatcher::HitIterator HMMMatcher::hit_begin() const {
    return { th_->hit };
}
const HMMMatcher::HitIterator HMMMatcher::hit_end() const {
    return { th_->hit + th_->N };
}

HMMMatcher::Hit::DomainIterator &HMMMatcher::Hit::DomainIterator::operator++() {
    ++d_;
    return *this;
}

bool HMMMatcher::Domain::reported() const { return d_->is_reported; }
bool HMMMatcher::Domain::included() const { return d_->is_included; }

std::pair<int, int> HMMMatcher::Domain::env() const { return { d_->ienv, d_->jenv }; }
std::pair<int, int> HMMMatcher::Domain::seqpos() const { return { d_->ad->sqfrom, d_->ad->sqto }; }
std::pair<int, int> HMMMatcher::Domain::hmmpos() const { return { d_->ad->hmmfrom, d_->ad->hmmto }; }
long HMMMatcher::Domain::L() const { return d_->ad->L; }
long HMMMatcher::Domain::M() const { return d_->ad->M; }

float HMMMatcher::Domain::bitscore() const { return d_->bitscore; }
double HMMMatcher::Domain::lnP() const { return d_->lnP; }
float HMMMatcher::Domain::oasc() const { return d_->oasc; }

const HMMMatcher::Hit::DomainIterator HMMMatcher::Hit::domain_begin() const {
    return { h_->dcl };
}
const HMMMatcher::Hit::DomainIterator HMMMatcher::Hit::domain_end() const {
    return { h_->dcl + h_->ndom };
}
