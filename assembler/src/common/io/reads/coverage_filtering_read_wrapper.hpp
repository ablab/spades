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

typedef qf::cqf<RtSeq> CQFKmerFilter;
//typedef CyclicHash<64, uint8_t, NDNASeqHash<uint8_t>> SeqHasher;
typedef SymmetricCyclicHash<uint8_t, uint64_t> SeqHasher;

template<class ReadType>
class CoverageFilteringReaderWrapper: public DelegatingWrapper<ReadType> {
    typedef std::shared_ptr<ReadStream<ReadType>> ReadStreamPtr;
public:
    explicit CoverageFilteringReaderWrapper(ReadStreamPtr reader,
                                            unsigned k,
                                            const CQFKmerFilter &filter, unsigned thr)
            : DelegatingWrapper<ReadType>(reader),
              k_(k), filter_(filter), thr_(thr),
              hasher_(k) {
        StepForward();
    }

    explicit CoverageFilteringReaderWrapper(const CoverageFilteringReaderWrapper &reader) = delete;
    void operator=(const CoverageFilteringReaderWrapper &reader) = delete;

    CoverageFilteringReaderWrapper& operator>>(ReadType& read) override {
        read = next_read_;
        StepForward();
        return *this;
    }

    void reset() override {
        this->reader().reset();
        StepForward();
    }

private:
    unsigned k_;
    const CQFKmerFilter &filter_;
    unsigned thr_;
    SeqHasher hasher_;

    ReadType next_read_;
    std::vector<unsigned> cov_;

    void StepForward() {
        while (!this->eof()) {
            this->reader() >> next_read_;

            // Calculate median coverage
            const Sequence &seq = next_read_.sequence();
            if (seq.size() < k_)
                continue;

            cov_.clear();
            RtSeq kmer = seq.start<RtSeq>(k_) >> 'A';
            SeqHasher::digest d;
            for (size_t j = k_ - 1; j < seq.size(); ++j) {
                kmer <<= seq[j];
                if (!kmer.IsMinimal())
                    d = hasher_.hash(!kmer);
                else
                    d = hasher_.hash(kmer);
                cov_.push_back(unsigned(filter_.lookup(d)));
            }

            size_t n = cov_.size() / 2;
            std::nth_element(cov_.begin(), cov_.begin() + n, cov_.end());
             if (cov_[n] >= thr_)
                return;
        }
    }
};

template<class ReadType>
inline std::shared_ptr<ReadStream<ReadType>> CovFilteringWrap(std::shared_ptr<ReadStream<ReadType>> reader_ptr,
                                                              unsigned k,
                                                              const CQFKmerFilter &filter, unsigned thr) {
    return std::make_shared<CoverageFilteringReaderWrapper<ReadType>>(reader_ptr, k, filter, thr);
}

template<class ReadType>
inline ReadStreamList<ReadType> CovFilteringWrap(ReadStreamList<ReadType> &readers,
                                                 unsigned k,
                                                 const CQFKmerFilter &filter, unsigned thr) {
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers.size(); ++i)
        answer.push_back(CovFilteringWrap(readers.ptr_at(i),
                                          k, filter, thr));

    return answer;
}

template<class Hasher>
unsigned CountMedianMlt(const Sequence &s, unsigned k, const CQFKmerFilter &kmer_mlt_index) {
    std::vector<unsigned> mlts;

    auto process_f = [&] (const RtSeq& /*kmer*/, uint64_t hash) {
        mlts.push_back(unsigned(kmer_mlt_index.lookup(hash)));
    };

    utils::KmerHashProcessor<Hasher> processor(process_f);
    processor.ProcessSequence(s, k);

    size_t n = mlts.size() / 2;
    std::nth_element(mlts.begin(), mlts.begin() + n, mlts.end());
    return mlts[n];
}

//template<class PairedReadType>
//inline ReadStreamList<PairedReadType> PairedCovFilteringWrap(const ReadStreamList<PairedReadType> &readers,
//                                                  unsigned k,
//                                                  unsigned thr) {
//    utils::StoringTypeFilter<utils::InvertableStoring> filter;
//    SeqHasher hasher(k);
//    auto single_readers = io::SquashingWrap<PairedReadType>(readers);
//
//    size_t kmers_cnt_est = EstimateCardinality(k, single_readers, filter);
//    auto cqf = make_shared<CQFKmerFilter>([&](const RtSeq &s) { return hasher.hash(s); },
//                                          kmers_cnt_est);
//
//    utils::FillCoverageHistogram(*cqf, k, single_readers, filter, thr + 1);
//
//    auto filter_f = [=] (io::PairedRead& p_r) { return CountMedianMlt(p_r.first().sequence(), k, hasher, *cqf) > thr ||
//                    CountMedianMlt(p_r.second().sequence(), k, hasher, *cqf) > thr; };
//    return io::FilteringWrap<PairedReadType>(readers, filter_f);
//}

}
