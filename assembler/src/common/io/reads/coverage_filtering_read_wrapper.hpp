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

template<class Hasher>
unsigned CountMedianMlt(const Sequence &s, unsigned k, const Hasher &hasher,
                        const utils::CQFKmerFilter &kmer_mlt_index) {
    if (s.size() < k)
        return 0;

    std::vector<unsigned> mlts;

    auto process_f = [&] (const RtSeq& /*kmer*/, uint64_t hash) {
        mlts.push_back(unsigned(kmer_mlt_index.lookup(hash)));
    };

    utils::KmerHashProcessor<Hasher> processor(hasher, process_f);
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

template<class ReadType, class Hasher>
class CoverageFilter;

//FIXME reduce code duplication, improve dispatching
template<class Hasher>
class CoverageFilter<SingleRead, Hasher> : public CoverageFilterBase<Hasher> {
    typedef CoverageFilterBase<Hasher> base;
public:
    CoverageFilter(unsigned k, const Hasher &hasher,
                   const utils::CQFKmerFilter &kmer_mlt_index,
                   unsigned thr) :
            base(k, hasher, kmer_mlt_index, thr) {}

    bool operator()(const SingleRead& r) const {
        this->CheckMedianMlt(r.sequence());
    }
};

template<class Hasher>
class CoverageFilter<SingleReadSeq, Hasher> : public CoverageFilterBase<Hasher> {
    typedef CoverageFilterBase<Hasher> base;
public:
    CoverageFilter(unsigned k, const Hasher &hasher,
                   const utils::CQFKmerFilter &kmer_mlt_index,
                   unsigned thr) :
            base(k, hasher, kmer_mlt_index, thr) {}

    bool operator()(const SingleReadSeq& r) const {
        this->CheckMedianMlt(r.sequence());
    }
};

template<class Hasher>
class CoverageFilter<PairedRead, Hasher> : public CoverageFilterBase<Hasher> {
    typedef CoverageFilterBase<Hasher> base;
public:
    CoverageFilter(unsigned k, const Hasher &hasher,
                   const utils::CQFKmerFilter &kmer_mlt_index,
                   unsigned thr) :
            base(k, hasher, kmer_mlt_index, thr) {}

    bool operator()(const PairedRead& r) const {
        this->CheckMedianMlt(r.first().sequence()) || this->CheckMedianMlt(r.second().sequence());
    }
};

template<class Hasher>
class CoverageFilter<PairedReadSeq, Hasher> : public CoverageFilterBase<Hasher> {
    typedef CoverageFilterBase<Hasher> base;
public:
    CoverageFilter(unsigned k, const Hasher &hasher,
                   const utils::CQFKmerFilter &kmer_mlt_index,
                   unsigned thr) :
            base(k, hasher, kmer_mlt_index, thr) {}

    bool operator()(const PairedReadSeq& r) const {
        this->CheckMedianMlt(r.first().sequence()) || this->CheckMedianMlt(r.second().sequence());
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
