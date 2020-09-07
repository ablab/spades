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
#include "io/reads/multifile_reader.hpp"

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
    kmers::KMerDiskStorage<RtSeq>
    BuildExtensionIndexFromStream(fs::TmpDir workdir, Index &index,
                                  Streams &streams,
                                  size_t read_buffer_size = 0) const {
        unsigned nthreads = (unsigned) streams.size();
        using KmerFilter = StoringTypeFilter<typename Index::storing_type>;

        // First, build a k+1-mer index
        using Splitter = DeBruijnReadKMerSplitter<typename Streams::ReadT, KmerFilter>;
        auto counter = std::make_unique<kmers::KMerDiskCounter<RtSeq>>(workdir,
                                                                       Splitter(workdir, index.k() + 1, streams, read_buffer_size));
        auto kmers = counter->CountAll(nthreads, nthreads, /* merge */ false);

        BuildExtensionIndexFromKPOMers(workdir, index, kmers,
                                       nthreads, read_buffer_size);

        return kmers;
    }

    template<class Index, class KMerStorage>
    void BuildExtensionIndexFromKPOMers(fs::TmpDir workdir,
                                        Index &index, const KMerStorage &kmers,
                                        unsigned nthreads, size_t read_buffer_size = 0) const {
        VERIFY(kmers.k() == index.k() + 1);

        // Now, count unique k-mers from k+1-mers
        using Splitter = DeBruijnKMerKMerSplitter<StoringTypeFilter<typename Index::storing_type>>;
        Splitter splitter(workdir, index.k(),
                          index.k() + 1, Index::storing_type::IsInvertable(), read_buffer_size);
        for (unsigned i = 0; i < kmers.num_buckets(); ++i)
            splitter.AddKMers(kmers.kmer_file(i));
        kmers::KMerDiskCounter<RtSeq> counter2(workdir, std::move(splitter));

        BuildIndex(index, counter2, 16, nthreads);

        // Build the kmer extensions
        INFO("Building k-mer extensions from k+1-mers");
#       pragma omp parallel for num_threads(nthreads)
        for (unsigned i = 0; i < nthreads; ++i)
            FillExtensionsFromIndex(kmers.kmer_file(i), index);
        INFO("Building k-mer extensions from k+1-mers finished.");
    }

private:
    DECL_LOGGER("DeBruijnExtensionIndexBuilder");
};

}
