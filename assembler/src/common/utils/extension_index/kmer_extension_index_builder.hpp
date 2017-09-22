//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "kmer_extension_index.hpp"
#include "utils/kmer_mph/kmer_index_builder.hpp"
#include "utils/kmer_mph/kmer_splitters.hpp"
#include "utils/kmer_counting.hpp"
#include "utils/ph_map/perfect_hash_map_builder.hpp"

namespace utils {

class DeBruijnExtensionIndexBuilder {
public:
    template<class ReadStream, class Index>
    void FillExtensionsFromStream(ReadStream &stream, Index &index) const {
        unsigned k = index.k();

        while (!stream.eof()) {
            typename ReadStream::read_type r;
            stream >> r;

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

public:
    template<class Index, class Streams>
    void BuildExtensionIndexFromStream(fs::TmpDir workdir, Index &index,
                                       Streams &streams, io::SingleStream *contigs_stream = 0,
                                       size_t read_buffer_size = 0) const {
        unsigned nthreads = (unsigned) streams.size();
        using KmerFilter = StoringTypeFilter<typename Index::storing_type>;

        // First, build a k+1-mer index
        DeBruijnReadKMerSplitter<typename Streams::ReadT, KmerFilter >
                splitter(workdir, index.k() + 1, 0xDEADBEEF, streams,
                         contigs_stream, read_buffer_size);
        KMerDiskCounter<RtSeq> counter(workdir, splitter);
        counter.CountAll(nthreads, nthreads, /* merge */ false);

        BuildExtensionIndexFromKPOMers(workdir, index, counter,
                                       nthreads, read_buffer_size);

    }

    template<class Index, class Counter>
    void BuildExtensionIndexFromKPOMers(fs::TmpDir workdir,
                                        Index &index, Counter &counter,
                                        unsigned nthreads, size_t read_buffer_size = 0) const {
        VERIFY(counter.k() == index.k() + 1);

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
