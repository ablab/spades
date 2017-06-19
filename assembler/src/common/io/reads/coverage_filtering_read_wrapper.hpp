//***************************************************************************
//* Copyright (c) 2017 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef COMMON_IO_COVFILTERINGREADERWRAPPER_HPP_
#define COMMON_IO_COVFILTERINGREADERWRAPPER_HPP_

#include "ireader.hpp"

#include "adt/cqf.hpp"
#include "adt/cyclichash.hpp"

#include <memory>

namespace io {

using CQFKmerFilter = qf::cqf<RtSeq>;

template<class ReadType>
class CoverageFilteringReaderWrapper: public DelegatingWrapper<ReadType> {
    typedef std::shared_ptr<ReadStream<ReadType>> ReadStreamPtr;
    using SeqHasher = CyclicHash<64, uint8_t, NDNASeqHash<uint8_t>>;
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
                cov_.push_back(unsigned(filter_.lookup(d, /* lock */ false)));
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

}

#endif /* COMMON_IO_COVFILTERINGREADERWRAPPER_HPP_ */
