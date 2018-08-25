//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/ph_map/perfect_hash_map.hpp"
#include "utils/kmer_mph/kmer_index.hpp"
#include <array>
#include <numeric>

namespace debruijn_graph {

class AlternativeEarlyTipClipper {
public:
    typedef utils::DeBruijnExtensionIndex<> Index;
    typedef Index::KMer Kmer;
    typedef Index::KeyWithHash KeyWithHash;

    AlternativeEarlyTipClipper(Index &index, size_t length_bound) : index_(index), length_bound_(length_bound) {}

    // Methods return the number of removed edges
    size_t ClipTips() {
        INFO("Early tip clipping");
        size_t result = ClipTips(10 * omp_get_max_threads());
        INFO(result << " " << (index_.k() + 1) << "-mers were removed by early tip clipper");
        return result;
    }

    // This interface is for MPI version of EarlyTC stage
    size_t ClipTips(size_t n_chunks, const std::vector<size_t> &chunks) {
        std::vector<Index::kmer_iterator> all_iters = index_.kmer_begin(n_chunks);
        std::vector<Index::kmer_iterator> iters;
        for (size_t chunk : chunks) {
            if (chunk < all_iters.size()) {  // all_iters.size() could be less than required
                iters.push_back(std::move(all_iters[chunk]));
            }
        }

        return ClipTips(iters);
    }

private:
    Index &index_;
    size_t length_bound_;

    /*
     * This method starts from the kmer that is second in the tip counting from the junction vertex
     * It returns all kmers of a tip into tip vector
     * In case it did not end as a tip or if it was too long tip the method returns empty vector
     */
    std::vector<KeyWithHash> FindForward_(KeyWithHash kh) {
        std::vector<KeyWithHash> tip;
        while (tip.size() < length_bound_ && index_.CheckUniqueIncoming(kh) && index_.CheckUniqueOutgoing(kh)) {
            tip.push_back(kh);
            kh = index_.GetUniqueOutgoing(kh);
        }
        tip.push_back(kh);
        if (index_.CheckUniqueIncoming(kh) && index_.IsDeadEnd(kh)) {
            return tip;
        } else {
            return {};
        }
    }

    size_t RemoveTip_(const vector<KeyWithHash> &tip) {
        for (const auto &kh : tip) index_.IsolateVertex(kh);
        return tip.size();
    }

    size_t RemoveTips_(const std::array<vector<KeyWithHash>, 4> &tips, size_t max) {
        size_t result = 0;
        for (const auto &tip : tips) {
            if (tip.size() < max) result += RemoveTip_(tip);
        }
        return result;
    }

    size_t RemoveForward_(const KeyWithHash &kh) {
        std::array<vector<KeyWithHash>, 4> tips;
        size_t max = 0;
        for (char c = 0; c < 4; c++) {
            if (index_.CheckOutgoing(kh, c)) {
                KeyWithHash khc = index_.GetOutgoing(kh, c);
                tips[c] = FindForward_(khc);
                size_t len = tips[c].empty() ? std::numeric_limits<size_t>::max() : tips[c].size();
                if (len > max) max = len;
            }
        }
        return RemoveTips_(tips, max);
    }

    size_t ClipTips(std::vector<Index::kmer_iterator> &iters) {
        std::vector<size_t> result(iters.size());
        std::vector<std::vector<KeyWithHash>> tipped_junctions(iters.size());
#pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < iters.size(); i++) {
            for (Index::kmer_iterator &it = iters[i]; it.good(); ++it) {
                auto seq = RtSeq(index_.k(), *it);
                for (const auto &s : {seq, !seq}) {  // The file stores only canonical (i.e., s > !s) k-mers
                    KeyWithHash kh = index_.ConstructKWH(s);
                    if (index_.OutgoingEdgeCount(kh) >= 2) {
                        size_t removed = RemoveForward_(kh);
                        result[i] += removed;
                        if (removed) {
                            tipped_junctions[i].push_back(kh);
                        }
                    }
                }
            }
        }

        size_t n_tipped_junctions = std::accumulate(tipped_junctions.cbegin(),
                                                    tipped_junctions.cend(),
                                                    size_t(0),
                                                    [](size_t sum, const auto &v) { return sum + v.size(); });
        INFO("#tipped junctions: " << n_tipped_junctions);

        // Remove links leading to tips
        size_t clipped_tips = 0;
#pragma omp parallel for schedule(guided) reduction(+:clipped_tips)
        for (size_t i = 0; i < tipped_junctions.size(); i++) {
            for (const KeyWithHash &kh : tipped_junctions[i]) {
                clipped_tips += RemoveUnpairedForwardLinks_(kh);
            }
        }
        INFO("Clipped tips: " << clipped_tips);

        size_t sum = std::accumulate(result.cbegin(), result.cend(), size_t(0));
        return sum;
    }

    size_t ClipTips(size_t n_chunks) {
        auto iters = index_.kmer_begin(n_chunks);
        return ClipTips(iters);
    }

    bool RemoveUnpairedForwardLink_(const KeyWithHash &kh, char ch) {
        if (index_.CheckOutgoing(kh, ch)) {
            KeyWithHash next_kh = index_.GetOutgoing(kh, ch);
            if (!index_.CheckIncoming(next_kh, kh[0])) {
                index_.DeleteOutgoing(kh, ch);
                return true;
            }
        }
        return false;
    }

    size_t RemoveUnpairedForwardLinks_(const KeyWithHash &kh) {
        size_t count = 0;
        for (char ch = 0; ch < 4; ++ch) {
            count += RemoveUnpairedForwardLink_(kh, ch);
        }
        return count;
    }

protected:
    DECL_LOGGER("Early tip clipping");
};

}  // namespace debruijn_graph
