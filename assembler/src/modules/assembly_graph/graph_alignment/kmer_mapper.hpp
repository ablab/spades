//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "data_structures/sequence/sequence_tools.hpp"
#include "data_structures/sequence/runtime_k.hpp"
#include "edge_index.hpp"

#include "htrie/hat-trie.h"

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
    size_t rawcnt_;
    bool verification_on_;

    hattrie_t *mapping_;

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

    value_t* internal_tryget(const Kmer &key) const {
        return hattrie_tryget(mapping_, (const char *)key.data(), rawcnt_ * sizeof(RawSeqData));
    }

    value_t* internal_get(const Kmer &key) const {
        return hattrie_get(mapping_, (const char *)key.data(), rawcnt_ * sizeof(RawSeqData));
    }

    int internal_erase(const Kmer &key) {
        return hattrie_del(mapping_, (const char *)key.data(), rawcnt_ * sizeof(RawSeqData));
    }

    void erase(const Kmer &key) {
        value_t *vp = internal_tryget(key);
        if (vp == nullptr)
            return;

        RawSeqData *value = reinterpret_cast<RawSeqData*>(*vp);
        delete[] value;
        int res = internal_erase(key);
        VERIFY_MSG(res == 0, "Failed to delete from kmer mapper");
    }

    void set(const Kmer &key, const Seq &value) {
        value_t *vp = internal_tryget(key);
        if (vp == nullptr) {
            vp = internal_get(key);
        } else {
            RawSeqData *value = reinterpret_cast<RawSeqData*>(*vp);
            delete[] value;
        }

        RawSeqData *rawvalue = new RawSeqData[rawcnt_];
        memcpy(rawvalue, value.data(), rawcnt_ * sizeof(RawSeqData));
        *vp = reinterpret_cast<uintptr_t>(rawvalue);
    }

    bool count(const Kmer &key) const {
        return internal_tryget(key) != nullptr;
    }

    const RawSeqData *find(const Kmer &key) const {
        value_t *vp = internal_tryget(key);
        if (vp == nullptr)
            return nullptr;

        return reinterpret_cast<const RawSeqData*>(*vp);
    }

    class KmerMapperIterator : public boost::iterator_facade<KmerMapperIterator,
                                                             const std::pair<Kmer, Seq>,
                                                             std::forward_iterator_tag,
                                                             const std::pair<Kmer, Seq>> {
      public:
        KmerMapperIterator(unsigned k, hattrie_iter_t *start = nullptr)
                : k_(k), iter_(start, [](hattrie_iter_t *p) { hattrie_iter_free(p); }) {}

      private:
        friend class boost::iterator_core_access;

        void increment() {
            hattrie_iter_next(iter_.get());
        }

        bool equal(const KmerMapperIterator &other) const {
            // Special case: NULL and finished are equal
            if (iter_.get() == nullptr || hattrie_iter_finished(iter_.get()))
                return other.iter_.get() == nullptr || hattrie_iter_finished(other.iter_.get());

            if (other.iter_.get() == nullptr)
                return false;

            return hattrie_iter_equal(iter_.get(), other.iter_.get());
        }

        const std::pair<Kmer, Seq> dereference() const {
            size_t len;
            Kmer k(k_, (const RawSeqData*)hattrie_iter_key(iter_.get(), &len));
            Seq s(k_, (const RawSeqData*)(*hattrie_iter_val(iter_.get())));
            return std::make_pair(k, s);
        }

        unsigned k_;
        std::shared_ptr<hattrie_iter_t> iter_;
    };

public:

    KmerMapper(const Graph &g, bool verification_on = true) :
            base(g, "KmerMapper"), k_(unsigned(g.k() + 1)), verification_on_(verification_on) {
        mapping_ = hattrie_create();
        rawcnt_ = Seq::GetDataSize(k_);
    }

    virtual ~KmerMapper() {
        hattrie_free(mapping_);
    }

    unsigned get_k() const { return k_; }

    KmerMapperIterator begin() const {
        return KmerMapperIterator(k_, hattrie_iter_begin(mapping_, false));
    }

    KmerMapperIterator end() const {
        return KmerMapperIterator(k_);
    }

    void Normalize() {
        std::vector<Kmer> all;
        for (auto it = begin(); it != end(); ++it) {
            all.push_back(it->first);
        }
        for (auto it = all.begin(); it != all.end(); ++it) {
            Normalize(*it);
        }
    }

    void Revert(const Kmer &kmer) {
        Kmer old_value = Substitute(kmer);
        if (old_value != kmer) {
            erase(kmer);
            set(old_value, kmer);
        }
    }

    void Normalize(const Kmer &kmer) {
        set(kmer, Substitute(kmer));
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
            if (count(new_kmer)) {
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
            if (count(new_kmer)) {
                if (verification_on_)
                    VERIFY(Substitute(new_kmer) == old_kmer);
                erase(new_kmer);
            }
            if (old_kmer != new_kmer)
                set(old_kmer, new_kmer);
        }
    }

    void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) override {
        VERIFY(this->g().EdgeNucls(new_edge) == this->g().EdgeNucls(edge2));
        RemapKmers(this->g().EdgeNucls(edge1), this->g().EdgeNucls(edge2));
    }

    Kmer Substitute(const Kmer &kmer) const {
        VERIFY(this->IsAttached());
        Kmer answer = kmer;
        const auto *rawval = find(answer);
        while (rawval != nullptr) {
            Seq val(k_, rawval);
            if (verification_on_)
                VERIFY(answer != val);

            answer = val;
            rawval = find(answer);
        }
        return answer;
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
            set(key, value);
        }
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
        // Delete all the values
        auto *iter = hattrie_iter_begin(mapping_, false);
        while (!hattrie_iter_finished(iter)) {
            RawSeqData *value = (RawSeqData*)(*hattrie_iter_val(iter));
            delete[] value;
            hattrie_iter_next(iter);
        }
        hattrie_iter_free(iter);
        // Delete the mapping and all the keys
        hattrie_clear(mapping_);
    }

    size_t size() const {
        return hattrie_size(mapping_);
    }

    // "turn on = true" means turning of all verifies
    void SetUnsafeMode(bool turn_on) {
        verification_on_ = !turn_on;
    }
};

}
