//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"

namespace utils {

typedef debruijn_graph::Graph Graph;
typedef debruijn_graph::EdgeId EdgeId;
typedef debruijn_graph::VertexId VertexId;

inline bool IsPathRegionCorrect(const std::pair<size_t, size_t> &region, size_t path_size) {
    return region.first < path_size && region.second < path_size;
}

inline Sequence GetSequenceOfPathRegion(const Graph &g, size_t k_value, const std::vector<EdgeId> &path,
                                        const std::pair<size_t, size_t> &region) {
    VERIFY(IsPathRegionCorrect(region, path.size()));

    if (region.first > region.second)
        return Sequence();

    EdgeId cur_edge = path[region.first];
    Sequence seq = g.EdgeNucls(cur_edge);

    for (auto i = region.first + 1; i <= region.second; ++i) {
        Sequence next_seq = g.EdgeNucls(path[i]);
        seq = seq + next_seq.Subseq(k_value, next_seq.size());
    }

    return seq;
}

inline Sequence GetSequenceByPath(const Graph &g, size_t k_value, const std::vector<EdgeId> &path) {
    if (path.size() == 0)
        return Sequence();
    return GetSequenceOfPathRegion(g, k_value, path, std::pair<size_t, size_t>(0, path.size() - 1));
}

inline size_t AlignmentOfSequencesByParts(const Sequence &seq1, const Sequence &seq2) {
    size_t max_length = 10000;
    if (std::min<size_t>(seq1.size(), seq2.size()) > max_length) {
        size_t shrink1 = max_length;
        size_t num_full_iter = seq1.size() / shrink1;

        size_t summary_dist = 0;
        size_t shrink2 = size_t((double(shrink1) / double(seq1.size())) * double(seq2.size()));
        for (size_t i = 0; i < num_full_iter; i++) {
            Sequence cur_seq1 = seq1.Subseq(shrink1 * i, shrink1 * (i + 1));
            Sequence cur_seq2 = seq2.Subseq(shrink2 * i, shrink2 * (i + 1));
            summary_dist += EditDistance(cur_seq1, cur_seq2);
        }

        if (seq1.size() % shrink1 != 0) {
            Sequence cur_seq1 = seq1.Subseq(shrink1 * num_full_iter, seq1.size());
            Sequence cur_seq2 = seq2.Subseq(shrink2 * num_full_iter, seq2.size());
            summary_dist += EditDistance(cur_seq1, cur_seq2);
        }

        return summary_dist;
    }
    return EditDistance(seq1, seq2);
}

inline double RelAlignmentOfSequences(const Sequence &seq1, const Sequence &seq2) {
    return double(AlignmentOfSequencesByParts(seq1, seq2)) / double(std::min<size_t>(seq1.size(), seq2.size()));
}

inline double RelativeLengthEquality(size_t len1, size_t len2) {
    return double(std::min<size_t>(len1, len2)) / double(std::max<size_t>(len1, len2));
}

}
