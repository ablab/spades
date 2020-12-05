//***************************************************************************
//* Copyright (c) 2014-2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/gap_info.hpp"

namespace path_extend {

typedef omnigraph::GapDescription<debruijn_graph::Graph> GapDescription;

class GapAnalyzer {

public:
    static const int INVALID_GAP = GapDescription::INVALID_GAP;
    GapAnalyzer(const debruijn_graph::Graph& g)
            : g_(g) { }

    virtual GapDescription FixGap(const GapDescription &gap) const = 0;
    virtual ~GapAnalyzer() = default;
protected:
    const debruijn_graph::Graph& g_;
};

class HammingGapAnalyzer: public GapAnalyzer {
    const double min_gap_score_;
    const size_t short_overlap_threshold_;
    const size_t basic_overlap_length_;

    static constexpr double MIN_OVERLAP_COEFF = 0.05;

    double ScoreGap(const Sequence& s1, const Sequence& s2) const;

public:
    //todo review parameters in usages
    HammingGapAnalyzer(const debruijn_graph::Graph& g,
                       double min_gap_score,
                       size_t short_overlap_threshold,
                       size_t basic_overlap_length):
                GapAnalyzer(g),
                min_gap_score_(min_gap_score),
                short_overlap_threshold_(short_overlap_threshold),
                basic_overlap_length_(basic_overlap_length) {
        DEBUG("HammingGapAnalyzer params: \n min_gap_score " << min_gap_score_ <<
              "\n short_overlap_threshold " << short_overlap_threshold_ <<
              "\n basic_overlap_length " << basic_overlap_length_);
    }

    GapDescription FixGap(const GapDescription &gap) const override;

private:
    DECL_LOGGER("HammingGapAnalyzer");
};

//LA stands for Local Alignment
//TODO if current setting will work -- get rid of flank_*_coefficient params
class LAGapAnalyzer: public GapAnalyzer {
public:
    LAGapAnalyzer(const debruijn_graph::Graph& g, size_t min_la_length,
                  double flank_multiplication_coefficient,
                  int flank_addition_coefficient) :
            GapAnalyzer(g),
            min_la_length_(min_la_length),
            flank_multiplication_coefficient_(flank_multiplication_coefficient),
            flank_addition_coefficient_(flank_addition_coefficient) {
        DEBUG("flank_multiplication_coefficient - " << flank_multiplication_coefficient_);
        DEBUG("flank_addition_coefficient  - " << flank_addition_coefficient_ );
    }

    GapDescription FixGap(const GapDescription &gap) const override;

private:
    DECL_LOGGER("LAGapAnalyzer");
    const size_t min_la_length_;
    const double flank_multiplication_coefficient_;
    const int flank_addition_coefficient_;

    static constexpr double IDENTITY_RATIO = 0.9;
    static constexpr double ESTIMATED_GAP_MULTIPLIER = 2.0;
    static constexpr size_t GAP_ADDITIONAL_COEFFICIENT = 30;
};

class CompositeGapAnalyzer: public GapAnalyzer {
public:

    CompositeGapAnalyzer(const debruijn_graph::Graph& g,
                         const std::vector<std::shared_ptr<GapAnalyzer>> &joiners,
                         size_t may_overlap_threshold, int must_overlap_threshold, size_t artificial_gap) :
            GapAnalyzer(g),
            joiners_(joiners),
            may_overlap_threshold_(may_overlap_threshold),
            must_overlap_threshold_(must_overlap_threshold),
            artificial_gap_(artificial_gap)
    {  }

    GapDescription FixGap(const GapDescription &gap) const override;

private:
    std::vector<std::shared_ptr<GapAnalyzer>> joiners_;
    const size_t may_overlap_threshold_;
    const int must_overlap_threshold_;
    const size_t artificial_gap_;

    DECL_LOGGER("CompositeGapAnalyzer");
};


}
