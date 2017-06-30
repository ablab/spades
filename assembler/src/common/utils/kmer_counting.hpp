#pragma once

#include "adt/cyclichash.hpp"
#include "adt/hll.hpp"
#include "adt/cqf.hpp"

namespace utils {

typedef qf::cqf<RtSeq> CQFKmerFilter;
typedef CyclicHash<64, uint8_t, NDNASeqHash<uint8_t>> SeqHasher;

template<class Hasher, class KmerFilter>
class KmerHashProcessor {
    typedef typename Hasher::digest HashT;
    typedef typename Hasher::char_t CharT;
    typedef std::function<void (const RtSeq&, HashT)> ProcessF;
    const Hasher hasher_;
    ProcessF process_f_;
    const KmerFilter filter_;
    RtSeq kmer_;
    HashT hash_;
public:
    KmerHashProcessor(const Hasher &hasher,
                      const ProcessF &process_f,
                      const KmerFilter &filter) :
            hasher_(hasher), process_f_(process_f), filter_(filter) {
    }

    void Init(const RtSeq &kmer) {
        kmer_ = kmer;
        hash_ = hasher_.hash(kmer);
    }

    void ProcessKmer(char inchar) {
        kmer_ <<= inchar;
        hash_ = hasher_.hash_update(hash_, (CharT) kmer_[0], (CharT) inchar);
        if (!filter_.filter(kmer_))
            return;
        process_f_(kmer_, hash_);
    }

    //todo use unsigned type for k
    void ProcessSequence(const Sequence &s, size_t k) {
        Init(s.start<RtSeq>(k) >> 'A');
        for (size_t j = k - 1; j < s.size(); ++j) {
            ProcessKmer(s[j]);
        }
    }

};

class HllProcessor {
    hll::hll<RtSeq> &hll_;
public:
    HllProcessor(hll::hll<RtSeq> &hll) : hll_(hll) {
    }

    void ProcessKmer(const RtSeq &/*kmer*/, uint64_t hash) {
        hll_.add(hash);
    }

    void Finalize() {}

};

template<class ReadStream, class Processor, class KmerFilter>
size_t FillFromStream(ReadStream &stream, Processor &processor, unsigned k,
                      const KmerFilter &filter = KmerFilter()) {
    size_t reads = 0;

    KmerHashProcessor<SeqHasher, KmerFilter> kmer_hash_processor(SeqHasher(k),
                                                                 [&](const RtSeq &kmer, uint64_t hash) { processor.ProcessKmer(kmer, hash); },
                                                                 filter);
    typename ReadStream::ReadT r;
    while (!stream.eof()) {
        stream >> r;
        reads += 1;

        const Sequence &seq = r.sequence();
        if (seq.size() < k)
            continue;

        kmer_hash_processor.ProcessSequence(seq, k);
        if (reads >= 1000000)
            break;

    }
    processor.Finalize();

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
        // accuired in the first attempt then insert the item in the
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

    void Finalize() {
        cqf_.merge(local_cqf_);
    }

};

//FIXME further reduce code duplication
template<class ReadStream, class KMerFilter>
size_t EstimateCardinality(unsigned k, ReadStream &streams,
                           const KMerFilter &filter) {
    unsigned nthreads = (unsigned) streams.size();
    SeqHasher hasher(k);

    std::vector<hll::hll<RtSeq>> hlls;
    hlls.reserve(nthreads);
    for (unsigned i = 0; i < nthreads; ++i)
        hlls.emplace_back([&](const RtSeq &s) { return hasher.hash(s); });

    streams.reset();
    size_t reads = 0, n = 15;
    while (!streams.eof()) {
#           pragma omp parallel for num_threads(nthreads) reduction(+:reads)
        for (unsigned i = 0; i < nthreads; ++i) {
            HllProcessor processor(hlls[omp_get_thread_num()]);
            reads += FillFromStream(streams[i], processor, k, filter);
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
    if (!res.second)
        return 256ull * 1024 * 1024;

    INFO("Estimated " << size_t(res.first) << " distinct kmers");
    return size_t(res.first);
}

template<class ReadStream, class KMerFilter>
void FillCoverageHistogram(qf::cqf<RtSeq> &cqf, unsigned k, ReadStream &streams,
                           const KMerFilter &filter, unsigned thr = 2) {
    unsigned nthreads = (unsigned) streams.size();
    SeqHasher hasher(k);

    // Create fallback per-thread CQF using same hash_size (important!) but different # of slots
    std::vector<qf::cqf<RtSeq>> local_cqfs;
    local_cqfs.reserve(nthreads);
    for (unsigned i = 0; i < nthreads; ++i)
        local_cqfs.emplace_back([&](const RtSeq &s) { return hasher.hash(s); },
                                1 << 16, cqf.hash_bits());

    INFO("Filtering threshold " << thr);
    streams.reset();
    size_t reads = 0, n = 15;
    while (!streams.eof()) {
#           pragma omp parallel for num_threads(nthreads) reduction(+:reads)
        for (unsigned i = 0; i < nthreads; ++i) {
            CQFProcessor processor(cqf, local_cqfs[i], thr);
            reads += FillFromStream(streams[i], processor, k, filter);
        }

        if (reads >> n) {
            INFO("Processed " << reads << " reads");
            n += 1;
        }
    }
    INFO("Total " << reads << " reads processed");
}

}
