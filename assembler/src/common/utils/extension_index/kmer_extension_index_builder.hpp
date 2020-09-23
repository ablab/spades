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


    template<class Index, class It>
    void FillExtensionsFromIndex(It begin, It end,
                                 Index &index) const {
        unsigned KPlusOne = index.k() + 1;
        for (; begin != end; ++begin) {
            RtSeq kpomer(KPlusOne, begin->first); // FIXME: temporary until boophm refactoring

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
        kmers::KMerDiskCounter<RtSeq> counter(workdir,
                                              Splitter(workdir, index.k() + 1, streams, read_buffer_size));
        auto kmers = counter.Count(10 * nthreads, nthreads);

        BuildExtensionIndexFromKPOMers(workdir, index, kmers,
                                       nthreads, read_buffer_size);

        return kmers;
    }

    template<class Index, class KMerStorage>
    void BuildExtensionIndexFromKPOMers(fs::TmpDir workdir,
                                        Index &index, const KMerStorage &kpomers,
                                        unsigned nthreads, size_t read_buffer_size = 0) const {
        VERIFY(kpomers.k() == index.k() + 1);

        // Now, count unique k-mers from k+1-mers
        using Splitter = DeBruijnKMerKMerSplitter<StoringTypeFilter<typename Index::storing_type>,
                                                  typename KMerStorage::kmer_iterator>;
        Splitter splitter(workdir, index.k(),
                          index.k() + 1, Index::storing_type::IsInvertable(), read_buffer_size);
        for (unsigned i = 0; i < kpomers.num_buckets(); ++i)
            splitter.AddKMers(adt::make_range(kpomers.bucket_begin(i), kpomers.bucket_end(i)));
        kmers::KMerDiskCounter<RtSeq> counter(workdir, std::move(splitter));

        BuildIndex(index, counter, kpomers.num_buckets(), nthreads);

        // Build the kmer extensions
        INFO("Building k-mer extensions from k+1-mers");
#       pragma omp parallel for num_threads(nthreads)
        for (size_t i = 0; i < kpomers.num_buckets(); ++i)
            FillExtensionsFromIndex(kpomers.bucket_begin(i), kpomers.bucket_end(i),
                                    index);
        INFO("Building k-mer extensions from k+1-mers finished.");
    }

private:
    DECL_LOGGER("DeBruijnExtensionIndexBuilder");
};

}
