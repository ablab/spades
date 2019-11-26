#pragma once

#include "adt/cyclichash.hpp"
#include "adt/hll.hpp"
#include "adt/cqf.hpp"
#include "ph_map/storing_traits.hpp"
#include "io/reads/read_processor.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/logger/logger.hpp"

namespace utils {

typedef qf::cqf CQFKmerFilter;

template<class Hasher, class KmerProcessor, class KmerFilter = StoringTypeFilter<SimpleStoring>>
class KmerSequenceProcessor {
    typedef uint64_t HashT;
    typedef typename rolling_hash::chartype CharT;
    Hasher hasher_;
    KmerProcessor &processor_;
    const KmerFilter filter_;

public:
    KmerSequenceProcessor(const Hasher &hasher, KmerProcessor &processor,
                      const KmerFilter &filter = StoringTypeFilter<SimpleStoring>()) :
            hasher_(hasher), processor_(processor), filter_(filter) {
    }

    void ProcessSequence(const Sequence &s, unsigned k) {
        RtSeq kmer = s.start<RtSeq>(k) >> 'A';
        auto hash = hasher_.hash(kmer);
        for (size_t j = k - 1; j < s.size(); ++j) {
            CharT inchar = (CharT) s[j];
            hash = hasher_.hash_update(hash, (CharT) kmer[0], inchar);
            //TODO can be optimized, no need to shift the kmer
            kmer <<= inchar;
            if (!filter_.filter(kmer))
                continue;
            processor_.ProcessKmer(kmer, (HashT) hash);
        }
    }

};

class HllProcessor {
    hll::hll<> &hll_;
public:
    HllProcessor(hll::hll<> &hll)
            : hll_(hll) { }

    void ProcessKmer(const RtSeq &/*kmer*/, uint64_t hash) {
        hll_.add(hash);
    }

};



template<class ReadStream, class SeqHasher, class Processor, class KmerFilter>
size_t FillFromStream(ReadStream &stream, const SeqHasher &hasher,
                      Processor &processor, unsigned k,
                      size_t max_read_cnt = std::numeric_limits<size_t>::max(),
                      const KmerFilter &filter = KmerFilter()) {
    size_t reads = 0;

    KmerSequenceProcessor<SeqHasher, Processor, KmerFilter>
            kmer_hash_processor(hasher, processor, filter);

    typename ReadStream::ReadT r;
    while (!stream.eof()) {
        stream >> r;
        reads += 1;

        const Sequence &seq = r.sequence();
        if (seq.size() < k)
            continue;

        kmer_hash_processor.ProcessSequence(seq, k);
        if (reads >= max_read_cnt)
            break;
    }

    return reads;
}

class CQFProcessor {
    CQFKmerFilter &cqf_;
    CQFKmerFilter &local_cqf_;
    const unsigned thr_;
public:
    CQFProcessor(CQFKmerFilter &cqf,
                 CQFKmerFilter &local_cqf,
                 unsigned thr) :
            cqf_(cqf), local_cqf_(local_cqf), thr_(thr) {
    }

    void ProcessKmer(const RtSeq &/*kmer*/, uint64_t hash) {
        // First try and insert in the main QF. If lock can't be
        // acquired in the first attempt then insert the item in the
        // local QF.
        if (cqf_.lookup(hash, /* lock */ true) >= thr_)
            return;

        if (!cqf_.add(hash, /* count */ 1,
                /* lock */ true, /* spin */ false)) {
            local_cqf_.add(hash, /* count */ 1,
                    /* lock */ false, /* spin */ false);
            if (local_cqf_.insertions() > local_cqf_.slots() / 2)
                cqf_.merge(local_cqf_);
        }
    }

};

template<class Hasher, class KMerFilter = utils::StoringTypeFilter<utils::SimpleStoring>>
class HllFiller {
 private:
    std::vector<HllProcessor> processors;
    std::vector<KmerSequenceProcessor<Hasher, HllProcessor, KMerFilter>> kmer_seq_processors;
    unsigned k;
    std::vector<size_t> reads;

    unsigned stop_after_log_read_cnt = 15;
 public:
    HllFiller(std::vector<hll::hll<>>& hlls, const Hasher& hasher, const KMerFilter& filter,
              unsigned k) {
        size_t nthreads = hlls.size();
        reads.resize(nthreads, 0);
        processors.reserve(nthreads);
        for (unsigned i = 0; i < nthreads; ++i) {
            processors.push_back(HllProcessor(hlls[i]));
            kmer_seq_processors.push_back(KmerSequenceProcessor<Hasher, HllProcessor, KMerFilter>
                                              (hasher, processors[i], filter));
            this->k = k;
        }
    }

    //Return value: should we interrupt reads processing
    template <class Read>
    bool operator()(std::unique_ptr<Read> r) {
        unsigned thread_id = (unsigned)omp_get_thread_num();
        reads[thread_id] += 1;
        const Sequence &seq = r->sequence();
        if (seq.size() < k) {
            return false;
        }
        kmer_seq_processors[thread_id].ProcessSequence(seq, k);

        if (reads[thread_id] >> stop_after_log_read_cnt) {
#           pragma omp atomic
            stop_after_log_read_cnt += 1;
            return true;
        }

        return false;
    }

    size_t processed_reads() const {
        size_t sum_reads = 0;
        for (size_t i = 0; i < reads.size(); ++i) {
            sum_reads += reads[i];
        }
        return sum_reads;
    }
};

template<class ReadStream, class Hasher, class KMerFilter = utils::StoringTypeFilter<utils::SimpleStoring>>
size_t EstimateCardinalityForOneStream(unsigned k, ReadStream &streams, const Hasher &hasher,
                           const KMerFilter &filter = utils::StoringTypeFilter<utils::SimpleStoring>()) {
    streams.reset();
    unsigned nthreads = (unsigned)omp_get_max_threads();
    std::vector<hll::hll<>> hlls(nthreads);
    size_t n = 15, reads = 0;
    HllFiller<Hasher, KMerFilter> hll_filler(hlls, hasher, filter, k);

    for (size_t i = 0; i < streams.size(); ++i) {
        while (!streams[i].eof()) {
            hammer::ReadProcessor rp(nthreads);
            rp.Run(streams[i], hll_filler);

            reads = hll_filler.processed_reads();
            if (reads >> n) {
                INFO("Processed " << reads << " reads");
                n += 1;
            }
        }
    }
    INFO("Total " << reads << " reads processed");

    for (size_t i = 1; i < hlls.size(); ++i) {
        hlls[0].merge(hlls[i]);
        hlls[i].clear();
    }

    double res = hlls[0].cardinality();

    INFO("Estimated " << size_t(res) << " distinct kmers");
    return size_t(res);
}

template<class ReadStream, class Hasher, class KMerFilter = utils::StoringTypeFilter<utils::SimpleStoring>>
size_t EstimateCardinalityUpperBound(unsigned k, ReadStream &streams, const Hasher &hasher,
                           const KMerFilter &filter = utils::StoringTypeFilter<utils::SimpleStoring>()) {
    unsigned stream_num = unsigned(streams.size());
    std::vector<hll::hll<>> hlls(stream_num);
    std::vector<HllProcessor> processors;
    for (size_t i = 0; i < hlls.size(); ++i) {
        processors.push_back(HllProcessor(hlls[i]));
    }

    streams.reset();
    size_t reads = 0, n = 15;
    while (!streams.eof()) {
#       pragma omp parallel for reduction(+:reads)
        for (unsigned i = 0; i < stream_num; ++i) {
            reads += FillFromStream(streams[i], hasher, processors[omp_get_thread_num()], k, 1000000, filter);
        }

        if (reads >> n) {
            INFO("Processed " << reads << " reads");
            n += 1;
        }
    }
    INFO("Total " << reads << " reads processed");

    for (size_t i = 1; i < hlls.size(); ++i) {
        hlls[0].merge(hlls[i]);
        hlls[i].clear();
    }

    double res = hlls[0].upper_bound_cardinality();

    INFO("Estimated " << size_t(res) << " distinct kmers");
    return size_t(res);
}

template<class Hasher, class ReadStream, class KMerFilter = utils::StoringTypeFilter<utils::SimpleStoring>>
void FillCoverageHistogram(qf::cqf &cqf, unsigned k, const Hasher &hasher, ReadStream &streams,
                           unsigned thr, const KMerFilter &filter = utils::StoringTypeFilter<utils::SimpleStoring>()) {
    unsigned stream_num = unsigned(streams.size());

    // Create fallback per-thread CQF using same hash_size (important!) but different # of slots
    std::vector<qf::cqf> local_cqfs;
    local_cqfs.reserve(stream_num);
    for (unsigned i = 0; i < stream_num; ++i)
        local_cqfs.emplace_back(1 << 16, cqf.hash_bits());

    INFO("Counting threshold " << thr);
    streams.reset();
    size_t reads = 0, n = 15;
    while (!streams.eof()) {
        #pragma omp parallel for reduction(+:reads)
        for (unsigned i = 0; i < stream_num; ++i) {
            CQFProcessor processor(cqf, local_cqfs[i], thr);
            reads += FillFromStream(streams[i], hasher, processor, k, 1000000, filter);
        }

        if (reads >> n) {
            INFO("Processed " << reads << " reads");
            n += 1;
        }
    }

    INFO("Merging local CQF");
    for (unsigned i = 0; i < stream_num; ++i) {
        cqf.merge(local_cqfs[i]);
    }

    INFO("Total " << reads << " reads processed");
}

}
