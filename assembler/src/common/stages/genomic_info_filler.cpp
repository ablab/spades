//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "genomic_info_filler.hpp"
#include "pipeline/genomic_info.hpp"
#include "assembly_graph/core/graph.hpp"
#include "modules/coverage_model/kmer_coverage_model.hpp"
#include "modules/simplification/ec_threshold_finder.hpp"

using namespace debruijn_graph;

static std::vector<size_t> extract(const std::map<size_t, size_t> &hist) {
    size_t maxcov = (hist.size() ? hist.rbegin()->first : 0);

    // Extract the values
    std::vector<size_t> res(maxcov);
    for (const auto &entry : hist) {
        if (entry.first < 1 || entry.first > maxcov)
            continue;

        res[entry.first - 1] = entry.second;
    }

    return res;
}

void GenomicInfoFiller::run(GraphPack &gp, const char*) {
    auto &ginfo = gp.get_mutable<GenomicInfo>();

    if (cfg::get().uneven_depth) {
        omnigraph::ErroneousConnectionThresholdFinder<Graph> finder(gp.get<Graph>());
        std::map<size_t, size_t> hist = finder.ConstructHistogram();
        double avg = finder.AvgCoverage();
        double gthr = finder.FindThreshold(hist);
        INFO("Average edge coverage: " << avg);
        INFO("Graph threshold: " << gthr);

        ginfo.set_cov_histogram(extract(hist));
        ginfo.set_ec_bound(std::min(avg, gthr));
    } else {
        // Fit the coverage model and get the threshold
        coverage_model::KMerCoverageModel CovModel(ginfo.cov_histogram(),
                                                   cfg::get().kcm.probability_threshold,
                                                   cfg::get().kcm.strong_probability_threshold);
        CovModel.Fit();

        ginfo.set_genome_size(CovModel.GetGenomeSize());
        ginfo.set_ec_bound((double)CovModel.GetErrorThreshold());
        if (CovModel.converged()) {
            ginfo.set_estimated_mean((double)CovModel.GetMeanCoverage());
            INFO("Mean coverage was calculated as " << ginfo.estimated_mean());
        } else
            INFO("Failed to estimate mean coverage");

        if (cfg::get().kcm.use_coverage_threshold) {
            VERIFY(math::gr(cfg::get().ds.aRL, 0.));
            double arl = math::gr(cfg::get().ds.aRL, double(gp.k())) ?
                         cfg::get().ds.aRL : double(cfg::get().ds.no_merge_RL);
            double coef = (arl - double(gp.k())) / arl;
            ginfo.set_trusted_bound(CovModel.converged() && cfg::get().kcm.coverage_threshold == 0.0 ?
                                       double(CovModel.GetLowThreshold()) :
                                       cfg::get().kcm.coverage_threshold * coef);
        }
    }

    INFO("EC coverage threshold value was calculated as " << ginfo.ec_bound());
    INFO("Trusted kmer low bound: " << ginfo.trusted_bound());
}
