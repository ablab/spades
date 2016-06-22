//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "data_structures/sequence/sequence_tools.hpp"
#include "data_structures/sequence/runtime_k.hpp"
#include "utils/adt/kmer_vector.hpp"
#include "edge_index.hpp"

#include "kmer_map.hpp"

#include <set>
#include <cstdlib>

namespace debruijn_graph {
template<class Graph>
class KmerMapper : public omnigraph::GraphActionHandler<Graph> {
    typedef omnigraph::GraphActionHandler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef runtime_k::RtSeq Kmer;
    typedef runtime_k::RtSeq Seq;
    typedef typename Seq::DataType RawSeqData;

    unsigned k_;
    KMerMap mapping_;
    bool verification_on_;
    bool normalized_;

    bool CheckAllDifferent(const Sequence &old_s, const Sequence &new_s) const {
        std::set<Kmer> kmers;
        Kmer kmer = old_s.start<Kmer>(k_) >> 0;
        for (size_t i = k_ - 1; i < old_s.size(); ++i) {
            kmer <<= old_s[i];
            kmers.insert(kmer);
        }
        kmer = new_s.start<Kmer>(k_) >> 0;
        for (size_t i = k_ - 1; i < new_s.size(); ++i) {
            kmer <<= new_s[i];
            kmers.insert(kmer);
        }
        return kmers.size() == old_s.size() - k_ + 1 + new_s.size() - k_ + 1;
    }

public:
    KmerMapper(const Graph &g, bool verification_on = true) :
            base(g, "KmerMapper"), k_(unsigned(g.k() + 1)), mapping_(k_), verification_on_(verification_on), normalized_(false) {
    }

    virtual ~KmerMapper() {}

    unsigned get_k() const { return k_; }

    auto begin() const -> decltype(mapping_.begin()) {
        return mapping_.begin();
    }

    auto end() const -> decltype(mapping_.end()) {
        return mapping_.end();
    }

    void Normalize() {
        if (normalized_)
            return;

        KMerVector<Kmer> all(k_, size());
        for (auto it = begin(); it != end(); ++it)
            all.push_back(it->first);

        for (auto it = all.begin(); it != all.end(); ++it) {
            Seq val(k_, it.data());
            Normalize(val);
        }
        normalized_ = true;
    }

    void Revert(const Kmer &kmer) {
        Kmer old_value = Substitute(kmer);
        if (old_value != kmer) {
            mapping_.erase(kmer);
            mapping_.set(old_value, kmer);
            normalized_ = false;
        }
    }

    void Normalize(const Kmer &kmer) {
        mapping_.set(kmer, Substitute(kmer));
    }

    bool CheckCanRemap(const Sequence &old_s, const Sequence &new_s) const {
        if (!CheckAllDifferent(old_s, new_s))
            return false;

        size_t old_length = old_s.size() - k_ + 1;
        size_t new_length = new_s.size() - k_ + 1;
        UniformPositionAligner aligner(old_s.size() - k_ + 1,
                                       new_s.size() - k_ + 1);
        Kmer old_kmer = old_s.start<Kmer>(k_);
        old_kmer >>= 0;
        for (size_t i = k_ - 1; i < old_s.size(); ++i) {
            old_kmer <<= old_s[i];
            size_t old_kmer_offset = i - k_ + 1;
            size_t new_kmer_offest = aligner.GetPosition(old_kmer_offset);
            if (old_kmer_offset * 2 + 1 == old_length && new_length % 2 == 0) {
                Kmer middle(k_ - 1, new_s, new_length / 2);
                if (typename Kmer::less2()(middle, !middle)) {
                    new_kmer_offest = new_length - 1 - new_kmer_offest;
                }
            }
            Kmer new_kmer(k_, new_s, new_kmer_offest);
            if (mapping_.count(new_kmer)) {
                if (Substitute(new_kmer) != old_kmer) {
                    return false;
                }
            }
        }
        return true;
    }

    void RemapKmers(const Sequence &old_s, const Sequence &new_s) {
        VERIFY(this->IsAttached());
        size_t old_length = old_s.size() - k_ + 1;
        size_t new_length = new_s.size() - k_ + 1;
        UniformPositionAligner aligner(old_s.size() - k_ + 1,
                                       new_s.size() - k_ + 1);
        Kmer old_kmer = old_s.start<Kmer>(k_);

        for (size_t i = k_ - 1; i < old_s.size(); ++i) {
            // Instead of shifting right
            if (i != k_ - 1) {
                old_kmer <<= old_s[i];
            }

            size_t old_kmer_offset = i - k_ + 1;
            size_t new_kmer_offest = aligner.GetPosition(old_kmer_offset);
            if (old_kmer_offset * 2 + 1 == old_length && new_length % 2 == 0) {
                Kmer middle(k_-1, new_s, new_length / 2);
                if (typename Kmer::less2()(middle, !middle)) {
                    new_kmer_offest = new_length - 1 - new_kmer_offest;
                }
            }
            Kmer new_kmer(k_, new_s, new_kmer_offest);
            if (mapping_.count(new_kmer)) {
                if (verification_on_)
                    VERIFY(Substitute(new_kmer) == old_kmer);
                mapping_.erase(new_kmer);
            }
            if (old_kmer != new_kmer) {
                mapping_.set(old_kmer, new_kmer);
                normalized_ = false;
            }
        }
    }

    void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) override {
        VERIFY(this->g().EdgeNucls(new_edge) == this->g().EdgeNucls(edge2));
        RemapKmers(this->g().EdgeNucls(edge1), this->g().EdgeNucls(edge2));
    }

    Kmer Substitute(const Kmer &kmer) const {
        VERIFY(this->IsAttached());
        Kmer answer = kmer;
        const auto *rawval = mapping_.find(answer);
        while (rawval != nullptr) {
            Seq val(k_, rawval);
            if (verification_on_)
                VERIFY(answer != val);

            answer = val;
            rawval = mapping_.find(answer);
        }
        return answer;
    }

    bool CanSubstitute(const Kmer &kmer) const {
        const auto *rawval = mapping_.find(kmer);
        return rawval != nullptr;
    }

    void BinWrite(std::ostream &file) const {
        uint32_t sz = (uint32_t)size();
        file.write((const char *) &sz, sizeof(uint32_t));

        for (auto iter = begin(); iter != end(); ++iter) {
            Kmer::BinWrite(file, iter->first);
            Kmer::BinWrite(file, iter->second);
        }
    }

    void BinRead(std::istream &file) {
        clear();

        uint32_t size;
        file.read((char *) &size, sizeof(uint32_t));
        for (uint32_t i = 0; i < size; ++i) {
            Kmer key(k_);
            Seq value(k_);
            Kmer::BinRead(file, &key);
            Seq::BinRead(file, &value);
            mapping_.set(key, value);
        }
        normalized_ = false;
    }

    bool CompareTo(KmerMapper<Graph> const &m) {
        if (size() != m.size()) {
            INFO("Unequal sizes");
            return false;
        }

        for (auto iter = begin(); iter != end(); ++iter) {
            auto cmp = m.mapping_.find(iter.first());
            if (cmp == m.mapping_.end() || cmp.second() != iter.second()) {
                return false;
            }
        }
        return true;
    }

    void clear() {
        normalized_ = false;
        return mapping_.clear();
    }

    size_t size() const {
        return mapping_.size();
    }

    // "turn on = true" means turning of all verifies
    void SetUnsafeMode(bool turn_on) {
        verification_on_ = !turn_on;
    }
};

}
