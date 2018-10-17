//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
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

class EarlyTipClipperProcessor {
public:
    typedef utils::DeBruijnExtensionIndex<> Index;
    typedef Index::KMer Kmer;
    typedef Index::KeyWithHash KeyWithHash;

    EarlyTipClipperProcessor(Index &index, size_t length_bound) : index_(index), length_bound_(length_bound) {}

    // Methods return the number of removed edges
    size_t ClipTips() {
        INFO("Early tip clipping");
        size_t result = ClipTips(10 * omp_get_max_threads());
        INFO(result << " " << (index_.k() + 1) << "-mers were removed by early tip clipper");
        return result;
    }

    size_t ClipTips(std::vector<Index::kmer_iterator> &iters) {
        std::vector<size_t> removed_kmers(iters.size());
        std::vector<std::vector<KeyWithHash>> tipped_junctions(iters.size());
#pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < iters.size(); i++) {
            TipsArray tips;
            for (auto &tip : tips) {
                tip.reserve(length_bound_);
            }

            for (Index::kmer_iterator &it = iters[i]; it.good(); ++it) {
                auto seq = RtSeq(index_.k(), *it);
                for (const auto &s : {seq, !seq}) {  // The file stores only canonical (i.e., s > !s) k-mers
                    KeyWithHash kh = index_.ConstructKWH(s);
                    if (index_.OutgoingEdgeCount(kh) >= 2) {
                        size_t removed = RemoveForward(kh, tips);
                        removed_kmers[i] += removed;
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

        // Remove inconsistent links leading to removed tips
        // Tips is already removed but their roots still store phantom extensions
        size_t clipped_tips = 0;
#pragma omp parallel for schedule(guided) reduction(+:clipped_tips)
        for (size_t i = 0; i < tipped_junctions.size(); i++) {
            for (const KeyWithHash &kh : tipped_junctions[i]) {
                clipped_tips += RemoveInconsistentForwardLinks(kh);
            }
        }
        INFO("Clipped tips: " << clipped_tips);

        size_t n_removed_kmers = std::accumulate(removed_kmers.cbegin(), removed_kmers.cend(), size_t(0));
        return n_removed_kmers;
    }


private:
    typedef std::array<std::vector<KeyWithHash>, 4> TipsArray;
    Index &index_;
    size_t length_bound_;

    /*
     * This method starts from the kmer that is second in the tip counting from the junction vertex
     * It records all kmers of a tip into tip vector
     * In case it did not end as a tip or if it was too long tip the method leaves the tip vector empty
     */
    void FindForward(KeyWithHash kh, std::vector<KeyWithHash> &tip) const {
        while (tip.size() < length_bound_ && index_.CheckUniqueIncoming(kh) && index_.CheckUniqueOutgoing(kh)) {
            tip.push_back(kh);
            kh = index_.GetUniqueOutgoing(kh);
        }
        tip.push_back(kh);
        if (!index_.CheckUniqueIncoming(kh) || !index_.IsDeadEnd(kh)) {
            // Branching or too long tip -> clear output vector
            tip.clear();
        }
    }

    size_t RemoveTip(const std::vector<KeyWithHash> &tip) {
        for (const auto &kh : tip) {
            index_.IsolateVertex(kh);
        }
        return tip.size();
    }

    size_t RemoveTips(const TipsArray &tips, size_t max) {
        size_t result = 0;
        for (const auto &tip : tips) {
            if (tip.size() < max) {
                result += RemoveTip(tip);
            }
        }
        return result;
    }

    size_t RemoveForward(const KeyWithHash &kh, TipsArray &tips) {
        size_t max = 0;
        auto mask = index_.get_value(kh);
        for (char c = 0; c < 4; c++) {
            tips[c].clear();
            if (mask.CheckOutgoing(c)) {
                KeyWithHash khc = index_.GetOutgoing(kh, c);
                FindForward(khc, tips[c]);
                size_t len = tips[c].empty() ? std::numeric_limits<size_t>::max() : tips[c].size();
                if (len > max) max = len;
            }
        }
        return RemoveTips(tips, max);
    }

    size_t ClipTips(size_t n_chunks) {
        auto iters = index_.kmer_begin(n_chunks);
        return ClipTips(iters);
    }

    size_t RemoveInconsistentForwardLinks(const KeyWithHash &kh) {
        size_t count = 0;
        auto mask = index_.get_value(kh);
        for (char c = 0; c < 4; ++c) {
            if (mask.CheckOutgoing(c)) {
                KeyWithHash next_kh = index_.GetOutgoing(kh, c);
                if (!index_.CheckIncoming(next_kh, kh[0])) {
                    index_.DeleteOutgoing(kh, c);
                    ++count;
                }
            }
        }
        return count;
    }

protected:
    DECL_LOGGER("Early tip clipping");
};

}  // namespace debruijn_graph
