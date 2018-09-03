#pragma once

#include "adt/cyclichash.hpp"
#include "adt/hll.hpp"
#include "adt/cqf.hpp"
#include "ph_map/storing_traits.hpp"
#include "common/utils/parallel/openmp_wrapper.h"

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

template<class ReadStream, class Hasher, class KMerFilter = utils::StoringTypeFilter<utils::SimpleStoring>>
size_t EstimateCardinality(unsigned k, ReadStream &streams, const Hasher &hasher,
                           const KMerFilter &filter = utils::StoringTypeFilter<utils::SimpleStoring>()) {
    unsigned stream_num = unsigned(streams.size());
    std::vector<hll::hll<>> hlls(stream_num);

    streams.reset();
    size_t reads = 0, n = 15;
    while (!streams.eof()) {
#       pragma omp parallel for reduction(+:reads)
        for (unsigned i = 0; i < stream_num; ++i) {
            HllProcessor processor(hlls[omp_get_thread_num()]);
            reads += FillFromStream(streams[i], hasher, processor, k, 1000000, filter);
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

    std::pair<double, bool> res = hlls[0].cardinality();
    if (!res.second) {
        INFO("Estimated " << size_t(res.first) << " distinct kmers (rough upper limit)");
        return 256ull * 1024 * 1024;
    }

    INFO("Estimated " << size_t(res.first) << " distinct kmers");
    return size_t(res.first);
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
