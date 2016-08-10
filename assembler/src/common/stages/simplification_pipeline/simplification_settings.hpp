//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "pipeline/config_struct.hpp"

namespace debruijn {

namespace simplification {

class LengthThresholdFinder {
public:
    static size_t MaxTipLength(size_t read_length, size_t k, double coeff) {
        return std::max((size_t) math::round((double)std::min(k, read_length / 2) * coeff),
                        read_length);
    }

    static size_t MaxBulgeLength(size_t k, double coeff,
                                 size_t additive_coeff) {
        return std::max((size_t) math::round((double)k * coeff), k + additive_coeff);
    }

    static size_t MaxErroneousConnectionLength(size_t k, size_t param) {
        return k + param;
    }

    static size_t MaxTipOriginatedECLength(size_t read_length, size_t k,
                                           double coeff) {
        return 2 * MaxTipLength(read_length, k, coeff) - 1;
    }
};

//todo use GenomicInfo as field!
class SimplifInfoContainer {
    size_t read_length_;
    double detected_mean_coverage_;
    double detected_coverage_bound_;
    bool main_iteration_;
    size_t chunk_cnt_;
    debruijn_graph::config::pipeline_type mode_;

public: 
    SimplifInfoContainer(debruijn_graph::config::pipeline_type mode) : 
        read_length_(-1ul),
        detected_mean_coverage_(-1.0),
        detected_coverage_bound_(-1.0),
        main_iteration_(false),
        chunk_cnt_(-1ul),
        mode_(mode) {
    }

    size_t read_length() const {
        VERIFY(read_length_ != -1ul);
        return read_length_;
    }

    double detected_mean_coverage() const {
        VERIFY(math::ge(detected_mean_coverage_, 0.));
        return detected_mean_coverage_;
    }

    double detected_coverage_bound() const {
        VERIFY(math::ge(detected_coverage_bound_, 0.));
        return detected_coverage_bound_;
    }

    bool main_iteration() const {
        return main_iteration_;
    }

    size_t chunk_cnt() const {
        VERIFY(chunk_cnt_ != -1ul);
        return chunk_cnt_;
    }

    debruijn_graph::config::pipeline_type mode() const {
        return mode_;
    }

    SimplifInfoContainer& set_read_length(size_t read_length) {
        read_length_ = read_length;
        return *this;
    }

    SimplifInfoContainer& set_detected_coverage_bound(double detected_coverage_bound) {
        detected_coverage_bound_ = detected_coverage_bound;
        return *this;
    }

    SimplifInfoContainer& set_detected_mean_coverage(double detected_mean_coverage) {
        detected_mean_coverage_ = detected_mean_coverage;
        return *this;
    }

    SimplifInfoContainer& set_main_iteration(bool main_iteration) {
        main_iteration_ = main_iteration;
        return *this;
    }

    SimplifInfoContainer& set_chunk_cnt(size_t chunk_cnt) {
        chunk_cnt_ = chunk_cnt;
        return *this;
    }
};

}

}
