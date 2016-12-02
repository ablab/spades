//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "genomic_info_filler.hpp"

#include "utils/coverage_model/kmer_coverage_model.hpp"
#include "modules/simplification/ec_threshold_finder.hpp"

#include "llvm/Support/YAMLTraits.h"
#include "llvm/Support/Errc.h"
#include "llvm/Support/FileSystem.h"

#include <string>

#include <map>
#include <vector>

using namespace llvm;
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

namespace llvm { namespace yaml {
template <>
struct MappingTraits<GenomicInfo> {
    static void mapping(yaml::IO &io, GenomicInfo &info) {
      info.yamlize(io);
    }
};


template <>
struct SequenceTraits<std::vector<std::size_t>> {
    static size_t size(IO &, std::vector<std::size_t> &seq) {
        return seq.size();
    }
    static size_t&
    element(IO &, std::vector<std::size_t> &seq, size_t index) {
        if (index >= seq.size())
            seq.resize(index+1);
        return seq[index];
    }
    static const bool flow = true;
};
}}

void GenomicInfo::yamlize(yaml::IO &io) {
    io.mapOptional("ec bound", ec_bound_, 0.0);
    io.mapOptional("estimated mean", estimated_mean_, 0.0);
    io.mapOptional("trusted bound", trusted_bound_, 0.0);
    io.mapOptional("genome size", genome_size_, size_t(0));
    io.mapOptional("coverage histogram", cov_histogram_);
}


bool GenomicInfo::Load(const std::string &filename) {
    ErrorOr<std::unique_ptr<MemoryBuffer>> Buf = MemoryBuffer::getFile(filename);
    if (!Buf)
        return false;

    yaml::Input yin(*Buf.get());
    yin >> *this;

    if (yin.error())
        return false;
    
    return true;
}

void GenomicInfo::Save(const std::string &filename) const {
    std::error_code EC;
    llvm::raw_fd_ostream ofs(filename, EC, llvm::sys::fs::OpenFlags::F_Text);
    llvm::yaml::Output yout(ofs);
    yout << const_cast<GenomicInfo&>(*this);
}

void GenomicInfoFiller::run(conj_graph_pack &gp, const char*) {
    if (cfg::get().uneven_depth) {
        ErroneousConnectionThresholdFinder<decltype(gp.g)> finder(gp.g);
        std::map<size_t, size_t> hist = finder.ConstructHistogram();
        double avg = finder.AvgCoverage();
        double gthr = finder.FindThreshold(hist);
        INFO("Average edge coverage: " << avg);
        INFO("Graph threshold: " << gthr);

        gp.ginfo.set_cov_histogram(extract(hist));
        gp.ginfo.set_ec_bound(std::min(avg, gthr));
    } else {
        // First, get k-mer coverage histogram
        std::map<size_t, size_t> tmp;
        size_t maxcov = 0;
        size_t kmer_per_record = 1;
        if (conj_graph_pack::index_t::InnerIndex::storing_type::IsInvertable())
            kmer_per_record = 2;

        for (auto I = gp.index.inner_index().value_cbegin(), E = gp.index.inner_index().value_cend(); I != E;  ++I) {
            size_t ccov = I->count;
            maxcov = std::max(ccov, maxcov);
            tmp[ccov] += kmer_per_record;
        }

        gp.ginfo.set_cov_histogram(extract(tmp));

        // Fit the coverage model and get the threshold
        utils::coverage_model::KMerCoverageModel CovModel(gp.ginfo.cov_histogram(), cfg::get().kcm.probability_threshold, cfg::get().kcm.strong_probability_threshold);
        CovModel.Fit();

        gp.ginfo.set_genome_size(CovModel.GetGenomeSize());
        gp.ginfo.set_ec_bound((double)CovModel.GetErrorThreshold());
        if (CovModel.converged()) {
            gp.ginfo.set_estimated_mean((double)CovModel.GetMeanCoverage());
            INFO("Mean coverage was calculated as " << gp.ginfo.estimated_mean());
        } else
            INFO("Failed to estimate mean coverage");

        if (cfg::get().kcm.use_coverage_threshold) {
            double coef = (cfg::get().ds.aRL() - double(cfg::get().K) + 1) / cfg::get().ds.aRL();
            if (coef < 0)
                coef = double(cfg::get().ds.RL() - cfg::get().K + 1) / double(cfg::get().ds.RL());
            gp.ginfo.set_trusted_bound(CovModel.converged() && cfg::get().kcm.coverage_threshold == 0.0 ?
                                       double(CovModel.GetLowThreshold()) :
                                       cfg::get().kcm.coverage_threshold * coef);
        }
    }

    INFO("EC coverage threshold value was calculated as " << gp.ginfo.ec_bound());
    INFO("Trusted kmer low bound: " << gp.ginfo.trusted_bound());
}
