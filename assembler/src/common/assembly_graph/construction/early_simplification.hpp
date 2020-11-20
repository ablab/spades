//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/extension_index/kmer_extension_index.hpp"
#include "utils/ph_map/perfect_hash_map.hpp"
#include "utils/kmer_mph/kmer_index.hpp"
#include "math/xmath.h"
#include <array>
#include <numeric>

namespace debruijn_graph {

static size_t RemoveInconsistentForwardLinks(utils::DeBruijnExtensionIndex<> &index,
                                             const utils::DeBruijnExtensionIndex<>::KeyWithHash &kh) {
    size_t count = 0;
    auto mask = index.get_value(kh);
    for (char c = 0; c < 4; ++c) {
        if (!mask.CheckOutgoing(c))
            continue;

        auto next_kh = index.GetOutgoing(kh, c);
        if (!index.CheckIncoming(next_kh, kh[0])) {
            index.DeleteOutgoing(kh, c);
            count += 1;
        }
    }
    return count;
}

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
                clipped_tips += RemoveInconsistentForwardLinks(index_, kh);
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
protected:
    DECL_LOGGER("Early tip clipping");
};


class EarlyLowComplexityClipperProcessor {
public:
    typedef utils::DeBruijnExtensionIndex<> Index;
    typedef Index::KMer Kmer;
    typedef Index::KeyWithHash KeyWithHash;

    EarlyLowComplexityClipperProcessor(Index &index,
                                       double at_ratio, size_t min_length, size_t max_length)
            : index_(index), ratio_(at_ratio), min_len_(min_length), max_len_(max_length) {}

    // @brief Removes low-complexity edges of length 1 (so effectively junction k-mers without unique extension).
    // @return Number of removed edges
    size_t RemoveATEdges() {
        INFO("Remove short poly A/T edges");
        auto iters = index_.kmer_begin(10 * omp_get_max_threads());
        size_t result = RemoveATEdges(iters);
        INFO(result << " " << (index_.k() + 1) << "-mers were removed by early poly A/T remover");
        return result;
    }

    size_t RemoveATEdges(std::vector<Index::kmer_iterator> &iters) {
        std::vector<std::vector<std::pair<KeyWithHash, char>>> at_edges(iters.size());
        double thr = index_.k() * ratio_;
#pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < iters.size(); i++) {
            for (Index::kmer_iterator &it = iters[i]; it.good(); ++it) {
                auto seq = RtSeq(index_.k(), *it);
                for (const auto &s : {seq, !seq}) {  // The file stores only canonical (i.e., s > !s) k-mers
                    KeyWithHash kh = index_.ConstructKWH(s);
                    auto extensions = index_.get_value(kh);

                    // Start from junction
                    if (!extensions.IsJunction())
                        continue;

                    // See, if this is a low-complexity k-mer
                    std::array<size_t, 4> counts = std::array<size_t, 4>();
                    for (size_t i = 0; i < index_.k(); ++i)
                        counts[s[i]] += 1;

                    size_t curm = *std::max_element(counts.begin(), counts.end());
                    if (math::ls(double(curm), thr))
                        continue;

                    for (char c = 0; c < 4; ++c) {
                        if (!extensions.CheckOutgoing(c))
                            continue;

                        auto next = index_.get_value(index_.GetOutgoing(kh, c));
                        // See, if this edge of length 1: the next should be junction or dead-end
                        if (!next.IsJunction() && !next.IsDeadEnd())
                            continue;

                        DEBUG("Removing edge: " << kh << "," << nucl(c) << "cnt: " << curm << ", thr: " << thr);

                        at_edges[i].emplace_back(kh, c);
                    }
                }
            }

        }

        size_t n_edges = std::accumulate(at_edges.cbegin(), at_edges.cend(),
                                         size_t(0),
                                         [](size_t sum, const auto &v) { return sum + v.size(); });

        size_t removed_links = 0;
        // Now, iterate over all A/T vertices removing incoming / outgoing links
#pragma omp parallel for schedule(guided) reduction(+:removed_links)
        for (size_t i = 0; i < at_edges.size(); ++i) {
            for (auto &edge : at_edges[i]) {
                if (!index_.CheckOutgoing(edge.first, edge.second))
                    continue;

                auto next = index_.GetOutgoing(edge.first, edge.second);
                index_.DeleteOutgoing(edge.first, edge.second);
                index_.DeleteIncoming(next, edge.first[0]);
                removed_links += 2;
            }
        }

#if 0
        for (size_t i = 0; i < at_edges.size(); ++i) {
            for (auto &edge : at_edges[i]) {
                VERIFY(!index_.CheckOutgoing(edge.first, edge.second));
                auto next = index_.GetOutgoing(edge.first, edge.second);
                VERIFY(!index_.CheckIncoming(next, edge.first[0]));
            }
        }
#endif

        INFO("Links removed: " << removed_links);
        return n_edges;
    }


    // @brief Removes low-complexity tips.
    // @return Number of removed edges
    size_t RemoveATTips() {
        INFO("Remove poly A/T tips");
        auto iters = index_.kmer_begin(10 * omp_get_max_threads());
        size_t result = RemoveATTips(iters);
        INFO(result << " " << (index_.k() + 1) << "-mers were removed by early poly A/T tip clipper");
        return result;
    }

    size_t RemoveATTips(std::vector<Index::kmer_iterator> &iters) {
        std::vector<std::vector<KeyWithHash>> at_edges(iters.size());
        size_t removed_kmers = 0;
        #pragma omp parallel for schedule(guided) reduction(+:removed_kmers)
        for (size_t i = 0; i < iters.size(); i++) {
            std::vector<KeyWithHash> tip;
            tip.reserve(max_len_);
            for (Index::kmer_iterator &it = iters[i]; it.good(); ++it) {
                auto seq = RtSeq(index_.k(), *it);
                for (const auto &s : {seq, !seq}) {  // The file stores only canonical (i.e., s > !s) k-mers
                    KeyWithHash kh = index_.ConstructKWH(s);
                    auto extensions = index_.get_value(kh);

                    // Start from tip ends
                    if (!extensions.IsDeadEnd() || !extensions.CheckUniqueIncoming())
                        continue;

                    // Go backward until the junction point counting the complexity
                    std::array<size_t, 4> counts = { 0, 0, 0, 0 };
                    tip.clear();

                    do {
                        tip.push_back(kh);
                        counts[kh[index_.k() - 1]] += 1;
                        kh = index_.GetUniqueIncoming(kh);
                    } while (tip.size() < max_len_ && !index_.IsJunction(kh));

                    // Now kh is a junction point (possible dead start) or still unique extension
                    // Bail out in the second case, as this tip is too long
                    if (!index_.IsDeadEnd(kh) && !index_.IsJunction(kh))
                        continue;

                    // If the tip is short, calculate the complexity up to min_len bp
                    for (size_t i = tip.size() - 1; i < min_len_; ++i)
                        counts[kh[index_.k() - 1 - i]] += 1;

                    // See, if this is a low-complexity tip
                    size_t curm = *std::max_element(counts.begin(), counts.end());
                    double thr = double(std::max(tip.size(), min_len_)) * ratio_;
                    if (math::ls(double(curm), thr))
                        continue;

                    at_edges[i].push_back(kh);
                    removed_kmers += tip.size();
                    for (const KeyWithHash &tip_kh : tip) {
                        index_.IsolateVertex(tip_kh);
                    }
                }
            }

        }

        // Remove inconsistent links leading to removed tips
        // Tips are already removed but their roots still store phantom extensions
        size_t clipped_tips = 0;
        #pragma omp parallel for schedule(guided) reduction(+:clipped_tips)
        for (size_t i = 0; i < at_edges.size(); i++) {
            for (const KeyWithHash &kh : at_edges[i]) {
                clipped_tips += RemoveInconsistentForwardLinks(index_, kh);
            }
        }
        INFO("Clipped tips: " << clipped_tips);

        return removed_kmers;
    }

private:
    Index &index_;
    double ratio_;
    size_t min_len_;
    size_t max_len_;

protected:
    DECL_LOGGER("Early A/T remover");
};

}  // namespace debruijn_graph
