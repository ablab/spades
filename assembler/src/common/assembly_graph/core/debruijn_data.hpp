//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <set>
#include <cstring>
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"
#include "sequence/sequence_tools.hpp"
#include "utils/standard_base.hpp"

namespace debruijn_graph {
class DeBruijnMaster;

class DeBruijnVertexData {
    friend class DeBruinMaster;
public:
    DeBruijnVertexData() {

    }
};

class CoverageData {
 private:
    unsigned coverage_;

 public:
    CoverageData()
            : coverage_(0) {
    }

    void inc_coverage(int value) {
        VERIFY(value >= 0 || coverage_ > unsigned(-value));
        coverage_ += value;
    }

    void set_coverage(unsigned coverage) {
        coverage_ = coverage;
    }

    //not length normalized
    unsigned coverage() const {
        return coverage_;
    }
};

class DeBruijnEdgeData {
    friend class DeBruinMaster;
    CoverageData coverage_;
    CoverageData flanking_cov_;
    Sequence nucls_;
public:

    DeBruijnEdgeData(const Sequence &nucls) :
            nucls_(nucls) {
    }

    const Sequence& nucls() const {
        return nucls_;
    }

    void inc_raw_coverage(int value) {
        coverage_.inc_coverage(value);
    }

    void set_raw_coverage(unsigned coverage) {
        coverage_.set_coverage(coverage);
    }

    unsigned raw_coverage() const {
        return coverage_.coverage();
    }

    void inc_flanking_coverage(int value) {
        flanking_cov_.inc_coverage(value);
    }

    void set_flanking_coverage(unsigned flanking_coverage) {
        flanking_cov_.set_coverage(flanking_coverage);
    }

    //not length normalized
    unsigned flanking_coverage() const {
        return flanking_cov_.coverage();
    }

    size_t size() const {
        return nucls_.size();
    }
};

class DeBruijnDataMaster {
private:
    const size_t k_;

public:
    typedef DeBruijnVertexData VertexData;
    typedef DeBruijnEdgeData EdgeData;

    DeBruijnDataMaster(size_t k) :
            k_(k) {
    }

    const EdgeData MergeData(const std::vector<const EdgeData*>& to_merge, bool safe_merging = true) const;

    std::pair<VertexData, std::pair<EdgeData, EdgeData>> SplitData(const EdgeData& edge, size_t position, bool is_self_conj = false) const;

    EdgeData GlueData(const EdgeData&, const EdgeData& data2) const;

    bool isSelfConjugate(const EdgeData &data) const {
        return data.nucls() == !(data.nucls());
    }

    EdgeData conjugate(const EdgeData &data) const {
        return EdgeData(!(data.nucls()));
    }

    VertexData conjugate(const VertexData & /*data*/) const {
        return VertexData();
    }

    size_t length(const EdgeData& data) const {
        return data.nucls().size() - k_;
    }

    size_t length(const VertexData& ) const {
        return k_;
    }

    size_t k() const {
        return k_;
    }

};

//typedef DeBruijnVertexData VertexData;
//typedef DeBruijnEdgeData EdgeData;
//typedef DeBruijnDataMaster DataMaster;

inline const DeBruijnEdgeData DeBruijnDataMaster::MergeData(const std::vector<const DeBruijnEdgeData*>& to_merge, bool safe_merging) const {
    std::vector<Sequence> ss;
    ss.reserve(to_merge.size());
    for (auto it = to_merge.begin(); it != to_merge.end(); ++it) {
        ss.push_back((*it)->nucls());
    }
    return EdgeData(MergeOverlappingSequences(ss, k_, safe_merging));
}

inline std::pair<DeBruijnVertexData, std::pair<DeBruijnEdgeData, DeBruijnEdgeData>> DeBruijnDataMaster::SplitData(const EdgeData& edge,
                                                                                                                  size_t position, 
                                                                                                                  bool is_self_conj) const {
    const Sequence& nucls = edge.nucls();
    size_t end = nucls.size();
    if (is_self_conj) {
        VERIFY(position < end);
        end -= position;
    }
    return std::make_pair(VertexData(), std::make_pair(EdgeData(edge.nucls().Subseq(0, position + k_)), EdgeData(nucls.Subseq(position, end))));
}

inline DeBruijnEdgeData DeBruijnDataMaster::GlueData(const DeBruijnEdgeData&, const DeBruijnEdgeData& data2) const {
    return data2;
}

}
