//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "kmer_extension_index.hpp"
#include "utils/kmer_mph/kmer_splitters.hpp"
#include "adt/cyclichash.hpp"
#include "adt/hll.hpp"
#include "adt/cqf.hpp"

namespace utils {

class DeBruijnExtensionIndexBuilder {
private:
    using CQFKmerFilter = qf::cqf<RtSeq>;
    using SeqHasher = CyclicHash<64, uint8_t, NDNASeqHash<uint8_t>>;

    template<class BaseFilter>
    class CQFKmerFilterWrapper {
      public:
        CQFKmerFilterWrapper(BaseFilter filter,
                             const CQFKmerFilter &cqf, unsigned thr)
                : filter_(filter), cqf_(cqf), thr_(thr) {}

        bool filter(const RtSeq &kmer) const {
            if (!filter_.filter(kmer))
                return false;
            return cqf_.lookup(kmer) >= thr_;
        }

      private:
        BaseFilter filter_;
        const CQFKmerFilter &cqf_;
        unsigned thr_;
    };

    template<class ReadStream, class KmerFilter>
    size_t FillHLLFromStream(ReadStream &stream, unsigned k,
                             hll::hll<RtSeq> &hll, const KmerFilter &filter = KmerFilter()) const {
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
                             const KmerFilter &filter = KmerFilter()) const {
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
                if (cqf.lookup(d, /* lock */ true) > thr)
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

public:
    template<class ReadStream, class Index>
    size_t FillExtensionsFromStream(ReadStream &stream, Index &index) const {
        unsigned k = index.k();
        size_t rl = 0;

        while (!stream.eof()) {
            typename ReadStream::ReadT r;
            stream >> r;
            rl = std::max(rl, r.size());

            const Sequence &seq = r.sequence();
            if (seq.size() < k + 1)
                continue;

            typename Index::KeyWithHash kwh = index.ConstructKWH(seq.start<RtSeq>(k));
            for (size_t j = k; j < seq.size(); ++j) {
                char nnucl = seq[j], pnucl = kwh[0];
                index.AddOutgoing(kwh, nnucl);
                kwh <<= nnucl;
                index.AddIncoming(kwh, pnucl);
            }
        }

        return rl;
    }

    template<class Index>
    void FillExtensionsFromIndex(const std::string &KPlusOneMersFilename,
                                 Index &index) const {
        unsigned KPlusOne = index.k() + 1;

        typename Index::kmer_iterator it(KPlusOneMersFilename,
                                         RtSeq::GetDataSize(KPlusOne));
        for (; it.good(); ++it) {
            RtSeq kpomer(KPlusOne, *it);

            char pnucl = kpomer[0], nnucl = kpomer[KPlusOne - 1];
            TRACE("processing k+1-mer " << kpomer);
            index.AddOutgoing(index.ConstructKWH(RtSeq(KPlusOne - 1, kpomer)),
                              nnucl);
            // FIXME: This is extremely ugly. Needs to add start / end methods to extract first / last N symbols...
            index.AddIncoming(index.ConstructKWH(RtSeq(KPlusOne - 1, kpomer << 0)),
                              pnucl);
        }
    }

    template<class ReadStream, class KMerFilter>
    size_t EstimateCardinality(unsigned k, ReadStream &streams,
                               const KMerFilter &filter) const {
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
                               const KMerFilter &filter, unsigned thr = 2) const {
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

public:
    template<class Index, class Streams>
    void BuildExtensionIndexFromStream(fs::TmpDir workdir, Index &index,
                                       Streams &streams, io::SingleStream *contigs_stream = 0,
                                       size_t read_buffer_size = 0) const {
        unsigned nthreads = (unsigned) streams.size();
        using KmerFilter = StoringTypeFilter<typename Index::storing_type>;

        std::vector<std::string> kmerfiles;
        if (1) {
            SeqHasher hasher(index.k() + 1);

            INFO("Estimating k-mers cardinality");
            size_t kmers = EstimateCardinality(index.k() + 1, streams, KmerFilter());

            // Create main CQF using # of slots derived from estimated # of k-mers
            qf::cqf<RtSeq> cqf([&](const RtSeq &s) { return hasher.hash(s); },
                               kmers);

            INFO("Building k-mer coverage histogram");
            unsigned thr = 2;
            FillCoverageHistogram(cqf, index.k() + 1, streams, KmerFilter(), thr);

            // First, build a k+1-mer index
            CQFKmerFilterWrapper<KmerFilter> filter(KmerFilter(), cqf, thr);
            DeBruijnReadKMerSplitter<typename Streams::ReadT, CQFKmerFilterWrapper<KmerFilter> >
                    splitter(workdir, index.k() + 1, 0xDEADBEEF, streams,
                             contigs_stream, read_buffer_size, filter);
            KMerDiskCounter<RtSeq> counter(workdir, splitter);
            counter.CountAll(nthreads, nthreads, /* merge */ false);

            BuildExtensionIndexFromKPOMers(workdir, index, counter,
                                           nthreads, read_buffer_size);
        } else {
            // First, build a k+1-mer index
            DeBruijnReadKMerSplitter<typename Streams::ReadT, KmerFilter >
                    splitter(workdir, index.k() + 1, 0xDEADBEEF, streams,
                             contigs_stream, read_buffer_size);
            KMerDiskCounter<RtSeq> counter(workdir, splitter);
            counter.CountAll(nthreads, nthreads, /* merge */ false);

            BuildExtensionIndexFromKPOMers(workdir, index, counter,
                                           nthreads, read_buffer_size);
        }
    }

    template<class Index, class Counter>
    void BuildExtensionIndexFromKPOMers(fs::TmpDir workdir,
                                        Index &index, Counter &counter,
                                        unsigned nthreads, size_t read_buffer_size = 0) const {
        VERIFY(counter.k() == index.k() + 1);
        std::vector<std::string> kmerfiles;

        // Now, count unique k-mers from k+1-mers
        DeBruijnKMerKMerSplitter<StoringTypeFilter<typename Index::storing_type> >
                splitter(workdir, index.k(),
                         index.k() + 1, Index::storing_type::IsInvertable(), read_buffer_size);
        for (unsigned i = 0; i < counter.num_buckets(); ++i)
            splitter.AddKMers(counter.GetMergedKMersFname(i));
        KMerDiskCounter<RtSeq> counter2(workdir, splitter);

        BuildIndex(index, counter2, 16, nthreads);

        // Build the kmer extensions
        INFO("Building k-mer extensions from k+1-mers");
#       pragma omp parallel for num_threads(nthreads)
        for (unsigned i = 0; i < nthreads; ++i)
            FillExtensionsFromIndex(counter.GetMergedKMersFname(i), index);
        INFO("Building k-mer extensions from k+1-mers finished.");
    }

private:
    DECL_LOGGER("DeBruijnExtensionIndexBuilder");
};

}
