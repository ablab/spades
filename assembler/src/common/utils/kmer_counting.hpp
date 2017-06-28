#pragma once

#include "adt/cyclichash.hpp"
#include "adt/hll.hpp"
#include "adt/cqf.hpp"

namespace utils {

typedef qf::cqf<RtSeq> CQFKmerFilter;
typedef CyclicHash<64, uint8_t, NDNASeqHash<uint8_t>> SeqHasher;

template<class ReadStream, class KmerFilter>
size_t FillHLLFromStream(ReadStream &stream, unsigned k,
                         hll::hll<RtSeq> &hll, const KmerFilter &filter = KmerFilter()) {
    size_t reads = 0;
    SeqHasher hasher(k);
    typename ReadStream::ReadT r;
    while (!stream.eof()) {
        stream >> r;
        reads += 1;

        const Sequence &seq = r.sequence();
        if (seq.size() < k)
            continue;

        RtSeq kmer = seq.start<RtSeq>(k) >> 'A';
        SeqHasher::digest d = hasher.hash(kmer);
        for (size_t j = k - 1; j < seq.size(); ++j) {
            uint8_t outchar = kmer[0];
            uint8_t inchar = seq[j];
            kmer <<= inchar;
            d = hasher.hash_update(d, outchar, inchar);

            if (!filter.filter(kmer))
                continue;

            hll.add(d);
        }

        if (reads >= 1000000)
            break;
    }

    return reads;
}

template<class ReadStream, class KmerFilter>
size_t FillCQFFromStream(ReadStream &stream, unsigned k, unsigned thr,
                         CQFKmerFilter &cqf, CQFKmerFilter &local_cqf,
                         const KmerFilter &filter = KmerFilter()) {
    size_t reads = 0;
    SeqHasher hasher(k);
    typename ReadStream::ReadT r;
    while (!stream.eof()) {
        stream >> r;
        reads += 1;

        const Sequence &seq = r.sequence();
        if (seq.size() < k)
            continue;

        RtSeq kmer = seq.start<RtSeq>(k) >> 'A';
        SeqHasher::digest d = hasher.hash(kmer);
        for (size_t j = k - 1; j < seq.size(); ++j) {
            uint8_t outchar = kmer[0];
            uint8_t inchar = seq[j];
            kmer <<= inchar;
            d = hasher.hash_update(d, outchar, inchar);
            if (!filter.filter(kmer))
                continue;

            // First try and insert in the main QF. If lock can't be
            // accuired in the first attempt then insert the item in the
            // local QF.
            if (cqf.lookup(d, /* lock */ true) >= thr)
                continue;

            if (!cqf.add(d, /* count */ 1,
                    /* lock */ true, /* spin */ false)) {
                local_cqf.add(d, /* count */ 1,
                        /* lock */ false, /* spin */ false);
                if (local_cqf.insertions() > local_cqf.slots() / 2)
                    cqf.merge(local_cqf);
            }
        }

        if (reads >= 1000000)
            break;
    }
    cqf.merge(local_cqf);

    return reads;
}

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
        for (unsigned i = 0; i < nthreads; ++i)
            reads += FillHLLFromStream(streams[i], k, hlls[omp_get_thread_num()],
                                       filter);

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
        for (unsigned i = 0; i < nthreads; ++i)
            reads += FillCQFFromStream(streams[i], k, thr,
                                       cqf, local_cqfs[i],
                                       filter);

        if (reads >> n) {
            INFO("Processed " << reads << " reads");
            n += 1;
        }
    }
    INFO("Total " << reads << " reads processed");
}

}
