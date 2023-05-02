//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "graph_core.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"
#include "sequence/sequence_tools.hpp"

#include <llvm/ADT/PointerSumType.h>
#include <llvm/ADT/PointerEmbeddedInt.h>

#include <utility>
#include <vector>
#include <set>
#include <cstring>
#include <cstdint>

namespace debruijn_graph {
class DeBruijnDataMaster;

class DeBruijnVertexData {
    friend class DeBruijnDataMaster;
    typedef size_t LinkId;

    enum OverlapKind {
        ComplexOverlap,
        ExplicitOverlap
    };

    struct OverlapStorage {
        OverlapStorage() = default;

        OverlapStorage(const std::vector<LinkId> &other_links)
                : links_(other_links) {}

        ~OverlapStorage() {}

        template<class It>
        void add_links(It from, It to) {
            while (from != to)
                links_.push_back(*from++);
        }

        void add_link(LinkId added_link) {
            links_.push_back(added_link);
        }

        const auto &links() const {
          return links_;
        }

        std::vector<LinkId> move() {
            auto links_copy = links_;
            links_.clear();
            return links_copy;
        }

        void clear() {
            links_.clear();
        }

        std::vector<LinkId> links_;
    };

    typedef llvm::PointerEmbeddedInt<uint32_t, 32> SimpleOverlap;
    typedef llvm::PointerSumType<OverlapKind,
                                 llvm::PointerSumTypeMember<ComplexOverlap, OverlapStorage*>,
                                 llvm::PointerSumTypeMember<ExplicitOverlap, SimpleOverlap>> Overlap;

    Overlap overlap_;

public:
    explicit DeBruijnVertexData(const std::vector<LinkId> &links)
            : overlap_(Overlap::create<ComplexOverlap>(new OverlapStorage(links))) {}

    explicit DeBruijnVertexData(unsigned overlap)
            : overlap_(Overlap::create<ExplicitOverlap>(overlap)) {}

    DeBruijnVertexData(DeBruijnVertexData&& that) {
        overlap_ = that.overlap_;
        that.overlap_.clear();
    }

    DeBruijnVertexData(const DeBruijnVertexData&) = delete;

    ~DeBruijnVertexData() {
        if (has_complex_overlap()) {
            delete complex_overlap();
            overlap_.clear();
        }
    }

    void set_overlap(unsigned overlap) {
        if (has_complex_overlap())
            delete complex_overlap();
        overlap_.set<ExplicitOverlap>(overlap);
    }

    unsigned overlap() const {
        return overlap_.get<ExplicitOverlap>();
    }

    auto links() const {
        return overlap_.get<ComplexOverlap>()->links();
    }

    auto move_links() {
        return overlap_.get<ComplexOverlap>()->move();
    }

    auto clear_links() {
        return overlap_.get<ComplexOverlap>()->clear();
    }

    void add_link(LinkId link) {
        overlap_.get<ComplexOverlap>()->add_link(link);
    }

    void add_links(const std::vector<LinkId> &links) {
        overlap_.get<ComplexOverlap>()->add_links(links.begin(), links.end());
    }

    bool has_complex_overlap() const {
        return overlap_.is<ComplexOverlap>();
    }

    OverlapStorage *complex_overlap() {
        return overlap_.get<ComplexOverlap>();
    }

    const OverlapStorage *complex_overlap() const {
        return overlap_.get<ComplexOverlap>();
    }

    DeBruijnVertexData clone() const {
        return has_complex_overlap() ?
                DeBruijnVertexData(links()) : DeBruijnVertexData(overlap());
    }
};

class CoverageData {
 private:
    uint32_t coverage_;

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
    friend class DeBruijnDataMaster;
    CoverageData coverage_;
    CoverageData flanking_cov_;
    Sequence nucls_;
public:

    explicit DeBruijnEdgeData(const Sequence &nucls) :
            nucls_(nucls) {}

    DeBruijnEdgeData(DeBruijnEdgeData&&) = default;
    DeBruijnEdgeData(const DeBruijnEdgeData&) = delete;

    const Sequence& nucls() const {
        return nucls_;
    }

    DeBruijnEdgeData clone() const {
        return DeBruijnEdgeData(nucls_);
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
    unsigned k_;

public:
    typedef DeBruijnVertexData VertexData;
    typedef DeBruijnEdgeData EdgeData;
    typedef DeBruijnVertexData::LinkId LinkId;
    typedef DeBruijnVertexData::OverlapStorage OverlapStorage;

    DeBruijnDataMaster(unsigned k)
            : k_(k) {}

    const EdgeData MergeData(const std::vector<const EdgeData *> &to_merge, const std::vector<uint32_t> &overlaps,
                             bool safe_merging = true) const;

    std::tuple<VertexData, EdgeData, EdgeData> SplitData(const EdgeData& edge, size_t position, bool is_self_conj = false) const;

    EdgeData GlueData(const EdgeData&, const EdgeData& data2) const;

    bool isSelfConjugate(const EdgeData &data) const {
        return data.nucls() == !(data.nucls());
    }

    EdgeData conjugate(const EdgeData &data) const {
        return EdgeData(!(data.nucls()));
    }

    VertexData conjugate(const VertexData &data) const {
        return data.clone();
    }

    size_t length(const EdgeData& data) const {
        return data.nucls().size() - k_;
    }

    // FIXME: make use of it!
    size_t length(const VertexData &data) const {
        return data.overlap();
    }

    unsigned k() const {
        return k_;
    }

    void set_k(unsigned k) {
        k_ = k;
    }
};

//typedef DeBruijnVertexData VertexData;
//typedef DeBruijnEdgeData EdgeData;
//typedef DeBruijnDataMaster DataMaster;

inline const DeBruijnEdgeData DeBruijnDataMaster::MergeData(const std::vector<const EdgeData *> &to_merge,
                                                            const std::vector<uint32_t> &overlaps,
                                                            bool safe_merging) const {
    std::vector<Sequence> ss;
    ss.reserve(to_merge.size());
    for (auto it = to_merge.begin(); it != to_merge.end(); ++it) {
        ss.push_back((*it)->nucls());
    }
    return EdgeData(MergeOverlappingSequences(ss, overlaps, safe_merging));
}

inline std::tuple<DeBruijnVertexData, DeBruijnEdgeData, DeBruijnEdgeData> DeBruijnDataMaster::SplitData(const EdgeData& edge,
                                                                                                        size_t position,
                                                                                                        bool is_self_conj) const {
    const Sequence& nucls = edge.nucls();
    size_t end = nucls.size();
    if (is_self_conj) {
        VERIFY(position < end);
        end -= position;
    }
    return { VertexData(k_), EdgeData(edge.nucls().Subseq(0, position + k_)), EdgeData(nucls.Subseq(position, end)) };
}

inline DeBruijnEdgeData DeBruijnDataMaster::GlueData(const DeBruijnEdgeData&, const DeBruijnEdgeData& data2) const {
    return data2.clone();
}

}
