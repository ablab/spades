//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace debruijn_graph {
class LinkCleaner {
private:
    typedef utils::DeBruijnExtensionIndex<> Index;
    typedef Index::KMer Kmer;
    typedef Index::KeyWithHash KeyWithHash;
    Index &index_;

    void CleanForwardLinks(const KeyWithHash &kh, char ch) {
        if (index_.CheckOutgoing(kh, ch)) {
            KeyWithHash next_kh = index_.GetOutgoing(kh, ch);
            if (!index_.CheckIncoming(next_kh, kh[0])) {
                index_.DeleteOutgoing(kh, ch);
            }
        }
    }

    void CleanBackwardLinks(const KeyWithHash &kh, char ch) {
        if (index_.CheckIncoming(kh, ch)) {
            KeyWithHash prev_kh = index_.GetIncoming(kh, ch);
            if (!index_.CheckOutgoing(prev_kh, kh[index_.k() - 1])) {
                index_.DeleteIncoming(kh, ch);
            }
        }
    }

public:
    LinkCleaner(Index &index) : index_(index) {}

    //TODO make parallel
    void CleanLinks() {
        vector<Index::kmer_iterator> iters = index_.kmer_begin(10 * omp_get_max_threads());
#pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < iters.size(); i++) {
            for (Index::kmer_iterator &it = iters[i]; it.good(); ++it) {
                const KeyWithHash kh = index_.ConstructKWH(RtSeq(index_.k(), *it));
                if (kh.is_minimal()) {
                    for (char ch = 0; ch < 4; ch++) {
                        CleanForwardLinks(kh, ch);
                        CleanBackwardLinks(kh, ch);
                    }
                }
            }
        }
    }
};

}  // namespace debruijn_graph
