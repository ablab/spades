#pragma once

#include <memory>

extern "C" {
    typedef struct p7_tophits_s P7_TOPHITS;
    typedef struct p7_pipeline_s P7_PIPELINE;
    typedef struct p7_profile_s P7_PROFILE;
    typedef struct p7_oprofile_s P7_OPROFILE;
    typedef struct p7_bg_s P7_BG;

}

namespace hmmer {

class HMM;

struct hmmer_cfg {
    bool acc;
    bool noali;
    double E; double T;
    double domE; double domT;
    double incE; double incT;
    double incdomE; double incdomT;
    bool cut_ga; bool cut_nc; bool cut_tc;
    bool max; double F1; double F2; double F3; bool nobias;

    hmmer_cfg()
            : acc(false), noali(false),
              E(10.0), T(0), domE(10.0), domT(0),
              incE(0.01), incT(0.0), incdomE(0.01), incdomT(0),
              cut_ga(false), cut_nc(false), cut_tc(false),
              max(false), F1(0.02), F2(1e-3), F3(1e-5), nobias(false)
    {}
};

class HMMMatcher {
  public:
    HMMMatcher(const HMM &hmmw,
               const hmmer_cfg &cfg);
    void match(const char *name, const char *seq, const char *desc = NULL);
    
    void summarize();
    P7_TOPHITS *hits() const;
    P7_PIPELINE *pipeline() const;

  private:

    P7_PIPELINE*
    pipeline_create(const hmmer_cfg &cfg, int M_hint, int L_hint, int long_targets, unsigned mode);

    std::unique_ptr<P7_PROFILE, void(*)(P7_PROFILE*)>  gm_;
    std::unique_ptr<P7_OPROFILE, void(*)(P7_OPROFILE*)> om_;  /* optimized query profile */
    std::unique_ptr<P7_BG, void(*)(P7_BG*)> bg_; /* null model */
    std::unique_ptr<P7_PIPELINE, void(*)(P7_PIPELINE*)> pli_; /* work pipeline */
    std::unique_ptr<P7_TOPHITS, void(*)(P7_TOPHITS*)> th_;
};


};

