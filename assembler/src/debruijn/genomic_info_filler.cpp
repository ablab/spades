//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "genomic_info_filler.hpp"

#include "graph_pack.hpp"
#include "kmer_coverage_model.hpp"
#include "omni/omni_tools.hpp"

#include <yaml-cpp/yaml.h>

#include <map>
#include <vector>

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

bool GenomicInfo::Load(const std::string &filename) {
    std::ifstream ifs(filename.c_str());
    if (!ifs)
        return false;

    YAML::Node node = YAML::Load(ifs);

    ec_bound_ = node["ec bound"].as<double>(0);
    trusted_bound_ = node["trusted bound"].as<size_t>(0);
    genome_size_ = node["genome size"].as<size_t>(0);
    cov_histogram_ = node["coverage histogram"].as<std::vector<size_t> >(std::vector<size_t>());

    return true;
}

void GenomicInfo::Save(const std::string &filename) const {
    std::ofstream ofs(filename.c_str());

    YAML::Node node;
    node["ec bound"] = ec_bound_;
    node["trusted bound"] = trusted_bound_;
    node["genome size"] = genome_size_;
    node["coverage histogram"] = cov_histogram_;

    ofs << node;
}

void GenomicInfoFiller::run(conj_graph_pack &gp, const char*) {
    if (cfg::get().ds.single_cell) {
        ErroneousConnectionThresholdFinder<decltype(gp.g)> finder(gp.g);
        std::map<size_t, size_t> hist = finder.ConstructHistogram();
        double avg = finder.AvgCoverage();
        double gthr = finder.FindThreshold(hist);
        INFO("Average edge coverage: " << avg);
        INFO("Graph threshold: " << gthr);

        gp.ginfo.set_cov_histogram(extract(hist));
        gp.ginfo.set_ec_bound(std::max(avg, gthr));
    } else {
        // First, get k-mer coverage histogram
        std::map<size_t, size_t> tmp;
        size_t maxcov = 0;
        size_t kmer_per_record = 1;
        if(conj_graph_pack::index_t::InnerIndexT::storing_type::IsInvertable()) {
            kmer_per_record = 2;
        }
        for (auto I = gp.index.inner_index().value_cbegin(), E = gp.index.inner_index().value_cend(); I != E;  ++I) {
            size_t ccov = I->count;
            maxcov = std::max(ccov, maxcov);
            tmp[ccov] += kmer_per_record;
        }

        gp.ginfo.set_cov_histogram(extract(tmp));

        // Fit the coverage model and get the threshold
        cov_model::KMerCoverageModel CovModel(gp.ginfo.cov_histogram(), cfg::get().kcm.probability_threshold, cfg::get().kcm.strong_probability_threshold);
        CovModel.Fit();

        gp.ginfo.set_genome_size(CovModel.GetGenomeSize());
        gp.ginfo.set_ec_bound((double)CovModel.GetErrorThreshold());
        // ginfo.set.trusted_bound(CovModel.GetLowThreshold());
    }
    INFO("EC coverage threshold value was calculated as " << gp.ginfo.ec_bound());
}
