//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "genomic_info_filler.hpp"
#include "modules/coverage_model/kmer_coverage_model.hpp"
#include "modules/simplification/ec_threshold_finder.hpp"

using namespace debruijn_graph;

static std::vector<size_t> extract(const std::map<size_t, size_t> &hist) {
    std::map<size_t, size_t> tmp = hist;

    size_t maxcov = 0;
    for (auto it = tmp.cbegin(), et = tmp.cend(); it != et; ++it)
        maxcov = std::max(maxcov, it->first);

    // Touch all the values until maxcov to make sure all the values exist in the map
    for (size_t i = 0; i <= maxcov; ++i)
        tmp[i];

    // Extract the values
    std::vector<size_t> res(maxcov);
    for (size_t i = 0; i < maxcov; ++i)
        res[i] = tmp[i + 1];

    return res;
}

void GenomicInfoFiller::run(conj_graph_pack &gp, const char*) {
    if (cfg::get().uneven_depth) {
        ErroneousConnectionThresholdFinder<conj_graph_pack::graph_t> finder(gp.g);
        std::map<size_t, size_t> hist = finder.ConstructHistogram();
        double avg = finder.AvgCoverage();
        double gthr = finder.FindThreshold(hist);
        INFO("Average edge coverage: " << avg);
        INFO("Graph threshold: " << gthr);

        gp.ginfo.set_cov_histogram(extract(hist));
        gp.ginfo.set_ec_bound(std::min(avg, gthr));
    } else {
        // Fit the coverage model and get the threshold
        coverage_model::KMerCoverageModel CovModel(gp.ginfo.cov_histogram(), cfg::get().kcm.probability_threshold, cfg::get().kcm.strong_probability_threshold);
        CovModel.Fit();

        gp.ginfo.set_genome_size(CovModel.GetGenomeSize());
        gp.ginfo.set_ec_bound((double)CovModel.GetErrorThreshold());
        if (CovModel.converged()) {
            gp.ginfo.set_estimated_mean((double)CovModel.GetMeanCoverage());
            INFO("Mean coverage was calculated as " << gp.ginfo.estimated_mean());
        } else
            INFO("Failed to estimate mean coverage");

        if (cfg::get().kcm.use_coverage_threshold) {
            VERIFY(math::gr(cfg::get().ds.aRL, 0.));
            double arl = math::gr(cfg::get().ds.aRL, double(gp.k_value)) ?
                         cfg::get().ds.aRL : double(cfg::get().ds.no_merge_RL);
            double coef = (arl - double(gp.k_value)) / arl;
            gp.ginfo.set_trusted_bound(CovModel.converged() && cfg::get().kcm.coverage_threshold == 0.0 ?
                                       double(CovModel.GetLowThreshold()) :
                                       cfg::get().kcm.coverage_threshold * coef);
        }
    }

    INFO("EC coverage threshold value was calculated as " << gp.ginfo.ec_bound());
    INFO("Trusted kmer low bound: " << gp.ginfo.trusted_bound());
}
