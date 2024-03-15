#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "common/utils/verify.hpp"

#include "aa_cursor.hpp"
#include "assembly_graph/core/graph.hpp"
#include "debruijn_graph_cursor.hpp"
#include "sequence/aa.hpp"

using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;
using Graph = debruijn_graph::ConjugateDeBruijnGraph;

class NaiveEdgeIndex {
    const static size_t SIZE = 1 << 20;
    static size_t hash(const std::string &kmer) { return std::hash<std::string>()(kmer) % SIZE; }

public:
    void build_nt(const Graph &graph, size_t k) {
        index_.resize(SIZE);
        k_ = k;
        size_t count = 0;
        for (auto it = graph.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            if (count % 100000 == 0) {
                INFO(count << " edges processed");
            }
            ++count;
            EdgeId edge = *it;
            std::string seq = graph.EdgeNucls(edge).str();
            size_t len = seq.length();
            VERIFY(len >= k);
            for (size_t i = 0; i < len - k + 1; ++i) {
                std::string kmer = seq.substr(i, k);
                index_[hash(kmer)].push_back({edge, i});
            }
        }
    }


    void build_aa(const Graph &graph, size_t k) {
        index_.resize(SIZE);
        k_ = k;
        size_t count = 0;
        for (auto it = graph.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            if (count % 100000 == 0) {
                INFO(count << " edges processed");
            }
            ++count;
            EdgeId edge = *it;
            std::string seq = graph.EdgeNucls(edge).str();
            for (int shift : {0, 1, 2}) {
                std::string seq_aa = aa::translate(seq.c_str() + shift);
                size_t len = seq_aa.length();
                for (size_t i = 0; i < len - k + 1; ++i) {
                    std::string kmer = seq_aa.substr(i, k);
                    index_[hash(kmer)].push_back({edge, 3 * i + shift});
                }
            }
        }
    }

    const std::vector<IdHolder> &candidates(const std::string &kmer) const {
        return index_[hash(kmer)];
    }

    size_t has_nt(const std::string &s, DebruijnGraphCursor::Context context) const {
        VERIFY(s.length() >= k_);
        std::vector<DebruijnGraphCursor> cursors;
        for (const auto &holder : candidates(s.substr(0, k_))) {
            cursors.push_back(DebruijnGraphCursor(holder.e(), holder.pos()));
        }

        return check(s, cursors, context);
    }

    size_t has_aa(const std::string &s, DebruijnGraphCursor::Context context) const {
        VERIFY(s.length() >= k_);
        std::vector<AAGraphCursor<DebruijnGraphCursor>> cursors;
        for (const auto &holder : candidates(s.substr(0, k_))) {
            DebruijnGraphCursor cursor(holder.e(), holder.pos());
            for (const auto &aa_cursor : make_aa_cursors(cursor, context)) {
                cursors.push_back(aa_cursor);
            }
        }

        return check(s, cursors, context);
    }

    template <typename Cursor>
    static size_t check(const std::string &s, std::vector<Cursor> cursors, typename Cursor::Context context) {
        for (size_t i = 0; i < s.length(); ++i) {
            // Filter
            auto it = std::remove_if(cursors.begin(), cursors.end(),
                                     [&](const auto &cursor) -> bool { return cursor.letter(context) != s[i]; });
            cursors.erase(it, cursors.end());
            // ->
            if (i != s.length() - 1) {
                std::vector<Cursor> nexts;
                for (const auto &cursor : cursors) {
                    for (const auto &next : cursor.next(context)) {
                        nexts.push_back(next);
                    }
                }
                cursors = nexts;
            }
        }
        return cursors.size();
    }

private:
    std::vector<std::vector<IdHolder>> index_;
    size_t k_;
};
