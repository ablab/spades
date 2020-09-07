#pragma once
//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "perfect_hash_map_builder.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include <cstdlib>

namespace utils {

struct CoverageHashMapBuilder : public utils::PerfectHashMapBuilder {
    template<class ReadStream, class Index>
    void FillCoverageFromStream(ReadStream &stream, Index &index) const {
        typedef typename Index::KeyType Kmer;
        unsigned k = index.k();

        while (!stream.eof()) {
            typename ReadStream::ReadT r;
            stream >> r;

            const Sequence &seq = r.sequence();
            if (seq.size() < k)
                continue;

            typename Index::KeyWithHash kwh = index.ConstructKWH(seq.start<Kmer>(k) >> 'A');
            for (size_t j = k - 1; j < seq.size(); ++j) {
                kwh <<= seq[j];
                if (!kwh.is_minimal() || !index.valid(kwh))
                    continue;

#                   pragma omp atomic
                index.get_raw_value_reference(kwh) += 1;
            }
        }
    }

    template<class Index, class KMerStorage, class Streams>
    void BuildIndex(Index &index,
                    const KMerStorage& storage,
                    Streams &streams) const {
        unsigned nthreads = (unsigned)streams.size();

        utils::PerfectHashMapBuilder::BuildIndex(index, storage, nthreads);
        INFO("Collecting k-mer coverage information from reads, this takes a while.");

        streams.reset();
#       pragma omp parallel for num_threads(nthreads)
        for (size_t i = 0; i < streams.size(); ++i) {
            FillCoverageFromStream(streams[i], index);
        }
    }
};
}
