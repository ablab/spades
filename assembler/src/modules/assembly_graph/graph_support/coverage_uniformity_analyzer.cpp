//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#include "coverage_uniformity_analyzer.hpp"
namespace debruijn_graph {
double CoverageUniformityAnalyzer::CountMedianCoverage() const{
    vector <pair<double, size_t> > coverages;
    size_t total_len = 0, short_len = 0, cur_len = 0;
    for (auto iter = g_.ConstEdgeBegin(); ! iter.IsEnd(); ++iter){
        if (g_.length(*iter) > length_bound_) {
            coverages.push_back(make_pair(g_.coverage(*iter), g_.length(*iter)));
            total_len += g_.length(*iter);
        } else {
            short_len += g_.length(*iter);
        }
    }
    if (total_len == 0){
        INFO("Median coverage detection failed, not enough long edges");
        return -1.0;
    }
    std::sort(coverages.begin(), coverages.end());
    size_t i = 0;
    while (cur_len < total_len/2 && i <coverages.size()) {
        cur_len += coverages[i].second;
        i++;
    }
    INFO ("genomic coverage is "<< coverages[i - 1].first << " calculated of length " << size_t (double(total_len) * 0.5));
    return coverages[i - 1].first;
}

std::pair<size_t, size_t> CoverageUniformityAnalyzer::TotalLengthsNearMedian(double allowed_variation, double median_coverage) const{
    std::pair<size_t, size_t> res(0,0);
    for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        if (g_.length(*iter) > length_bound_) {
            if (g_.coverage(*iter) < median_coverage * (1 + allowed_variation) &&
                g_.coverage(*iter) > median_coverage * (1 - allowed_variation)) {
                res.first += g_.length(*iter);
            } else {
                res.second += g_.length(*iter);
            }
        }
    }
    return res;
}

size_t CoverageUniformityAnalyzer::TotalLongEdgeLength() const {
    size_t res = 0;
    for (auto iter = g_.ConstEdgeBegin(); ! iter.IsEnd(); ++iter){
        if (g_.length(*iter) > length_bound_) {
            res += g_.length(*iter);
        }
    }
    return res;
}

double CoverageUniformityAnalyzer::UniformityFraction(double allowed_variation, double median_coverage) const {
    std::pair<size_t, size_t> lengths = TotalLengthsNearMedian(allowed_variation, median_coverage);
    size_t total_len = lengths.first + lengths.second;
    if (total_len == 0) {
        WARN(" No edges longer than length bound(" << length_bound_ <<" )");
        return 0;
    }
    return double(lengths.first) / double(total_len);
}

}