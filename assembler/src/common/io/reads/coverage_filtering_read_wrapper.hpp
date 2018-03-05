//***************************************************************************
//* Copyright (c) 2017 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <memory>
#include "io/reads/ireader.hpp"

#include "adt/cqf.hpp"
#include "adt/cyclichash.hpp"
#include "utils/kmer_counting.hpp"

namespace io {

class KmerMultiplicityFiller {
    std::vector<unsigned> &mlts_;
    const utils::CQFKmerFilter &kmer_mlt_index_;

public:
    KmerMultiplicityFiller(std::vector<unsigned> &mlts,
                           const utils::CQFKmerFilter &kmer_mlt_index) :
            mlts_(mlts), kmer_mlt_index_(kmer_mlt_index) {}

    void ProcessKmer(const RtSeq &/*kmer*/, uint64_t hash) {
        mlts_.push_back(unsigned(kmer_mlt_index_.lookup(hash)));
    }
};

template<class Hasher>
unsigned CountMedianMlt(const Sequence &s, unsigned k, const Hasher &hasher,
                        const utils::CQFKmerFilter &kmer_mlt_index) {
    if (s.size() < k)
        return 0;

    std::vector<unsigned> mlts;
    KmerMultiplicityFiller mlt_filler(mlts, kmer_mlt_index);
    utils::KmerSequenceProcessor<Hasher, KmerMultiplicityFiller> processor(hasher, mlt_filler);
    processor.ProcessSequence(s, k);

    size_t n = mlts.size() / 2;
    std::nth_element(mlts.begin(), mlts.begin() + n, mlts.end());
    return mlts[n];
}

template<class Hasher>
class CoverageFilterBase {
    const unsigned k_;
    const Hasher hasher_;
    const utils::CQFKmerFilter &kmer_mlt_index_;
    const unsigned thr_;
public:
    CoverageFilterBase(unsigned k, const Hasher &hasher,
                       const utils::CQFKmerFilter &kmer_mlt_index,
                       unsigned threshold) :
            k_(k), hasher_(hasher),
            kmer_mlt_index_(kmer_mlt_index), thr_(threshold) {
    }

    bool CheckMedianMlt(const Sequence &s) const {
        return CountMedianMlt(s, k_, hasher_, kmer_mlt_index_) >= thr_;
    }
};

template<class SingleReadType, class Hasher>
class CoverageFilter : public CoverageFilterBase<Hasher> {
    typedef CoverageFilterBase<Hasher> base;
public:
    CoverageFilter(unsigned k, const Hasher &hasher,
                   const utils::CQFKmerFilter &kmer_mlt_index,
                   unsigned thr) :
            base(k, hasher, kmer_mlt_index, thr) {}

    bool operator()(const SingleReadType &r) const {
        return this->CheckMedianMlt(r.sequence());
    }
};

template<class SingleReadType, class Hasher>
class CoverageFilter<UniversalPairedRead<SingleReadType>, Hasher> : public CoverageFilterBase<Hasher> {
    typedef CoverageFilterBase<Hasher> base;
    typedef UniversalPairedRead<SingleReadType> PairedReadType;
public:
    CoverageFilter(unsigned k, const Hasher &hasher,
                   const utils::CQFKmerFilter &kmer_mlt_index,
                   unsigned thr) :
            base(k, hasher, kmer_mlt_index, thr) {}

    bool operator()(const PairedReadType& r) const {
        return this->CheckMedianMlt(r.first().sequence()) || this->CheckMedianMlt(r.second().sequence());
    }

};

template<class ReadType, class Hasher>
inline std::shared_ptr<ReadStream<ReadType>> CovFilteringWrap(std::shared_ptr<ReadStream<ReadType>> reader_ptr,
                                                              unsigned k, const Hasher &hasher,
                                                              const utils::CQFKmerFilter &cqf, unsigned thr) {
    CoverageFilter<ReadType, Hasher> filter(k, hasher, cqf, thr);
    return io::FilteringWrap<ReadType>(reader_ptr, [=](const ReadType &r) { return filter(r); });
}

template<class ReadType, class Hasher>
inline ReadStreamList<ReadType> CovFilteringWrap(ReadStreamList<ReadType> &readers,
                                                 unsigned k, const Hasher &hasher,
                                                 const utils::CQFKmerFilter &filter, unsigned thr) {
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers.size(); ++i)
        answer.push_back(CovFilteringWrap(readers.ptr_at(i),
                                          k, hasher, filter, thr));

    return answer;
}

}
