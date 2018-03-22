#pragma once

#include <llvm/ADT/iterator.h>
#include <llvm/ADT/iterator_range.h>

#include <memory>
#include "hmmer_fwd.h"

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
    class HitIterator;
    class DomainIterator;

    class Domain {
      public:
        Domain(const P7_DOMAIN *d)
                : d_(d) {}

        bool reported() const;
        bool included() const;

        float bitscore() const;
        double lnP() const;
        float oasc() const;
        std::pair<int, int> env() const;
        std::pair<int, int> seqpos() const;
        std::pair<int, int> hmmpos() const;
        long L() const;
        long M() const;

        const P7_DOMAIN *domain() const { return d_; }

      private:
        const P7_DOMAIN *d_;
        friend class DomainIterator;
    };

    class Hit {
      public:
        class DomainIterator : public llvm::iterator_facade_base<DomainIterator,
                                                                 std::forward_iterator_tag,
                                                                 Domain> {
          public:
            DomainIterator(P7_DOMAIN *d)
                    : d_(d) {}

            Domain operator*() const { return { d_ }; }
            DomainIterator &operator++();
            bool operator==(const DomainIterator &that) const { return d_ == that.d_; }


          private:
            P7_DOMAIN *d_;
        };

        Hit(const P7_HIT *h)
                : h_(h) {}

        const char *name() const;
        const char *acc() const;
        const char *desc() const;
        float score() const;
        double lnP() const;
        size_t ndom() const;
        uint32_t flags() const;

        bool reported() const;
        bool included() const;

        const DomainIterator domain_begin() const;
        const DomainIterator domain_end() const;

        llvm::iterator_range<DomainIterator> domains() const {
            return llvm::make_range(domain_begin(), domain_end());
        }

        const P7_HIT *hit() const { return h_; }

      private:
        const P7_HIT *h_;
    };

    class HitIterator : public llvm::iterator_facade_base<HitIterator,
                                                          std::forward_iterator_tag,
                                                          Hit> {
      public:
        HitIterator(P7_HIT **hl)
                : hl_(hl) {}

        Hit operator*() const { return *hl_; }
        HitIterator &operator++() { ++hl_; return *this; }
        bool operator==(const HitIterator &that) const { return hl_ == that.hl_; }

      private:
        P7_HIT **hl_;
    };

    HMMMatcher(const HMM &hmmw,
               const hmmer_cfg &cfg);
    void match(const char *name, const char *seq, const char *desc = NULL);

    void summarize();
    P7_TOPHITS *top_hits() const;
    P7_PIPELINE *pipeline() const;

    const HitIterator hit_begin() const;
    const HitIterator hit_end() const;

    llvm::iterator_range<HitIterator> hits() const {
        return llvm::make_range(hit_begin(), hit_end());
    }

  private:
    P7_PIPELINE*
    pipeline_create(const hmmer_cfg &cfg,
                    int M_hint, int L_hint,
                    int long_targets, unsigned mode);

    std::unique_ptr<P7_PROFILE, void(*)(P7_PROFILE*)>  gm_;
    std::unique_ptr<P7_OPROFILE, void(*)(P7_OPROFILE*)> om_;  /* optimized query profile */
    std::unique_ptr<P7_BG, void(*)(P7_BG*)> bg_;              /* null model */
    std::unique_ptr<P7_PIPELINE, void(*)(P7_PIPELINE*)> pli_; /* work pipeline */
    std::unique_ptr<P7_TOPHITS, void(*)(P7_TOPHITS*)> th_;
};


};
