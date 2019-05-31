//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>
#include <bitset>

#include "common/utils/verify.hpp"

// Serialization
#include <cereal/types/common.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>

class CachedAACursorContext;

class CachedAACursor {
public:
    using Context = const CachedAACursorContext *;
    CachedAACursor(size_t index, unsigned char mask = 0b111) : index_{index}, mask_{mask} {}
    CachedAACursor() : index_{size_t(-1)}, mask_{0b000} {}
    bool is_empty() const { return mask_ == 0b000; }
    bool operator==(const CachedAACursor &other) const { return to_size_t() == other.to_size_t(); }
    const std::vector<CachedAACursor> &next(Context context) const;
    const std::vector<CachedAACursor> &prev(Context context) const;
    char letter(Context context) const;
    size_t index() const { return index_; }
    unsigned char mask() const { return mask_; }
    // uint64_t to_size_t() const { return *reinterpret_cast<const uint64_t *>(this); }
    uint64_t to_size_t() const { return (index_ << 3) + mask_; }
    const std::vector<CachedAACursor> &next_frame_shift(Context context) const;

    CachedAACursor triplet_form() const {
        CachedAACursor result = *this;
        result.mask_ = 0b111;
        return result;
    }

    std::vector<size_t> nucl_cursor_indices(Context context) const;

    bool operator<(const CachedAACursor &other) const { return to_size_t() < other.to_size_t(); }

    template <class Archive>
    void serialize(Archive &archive) {
        archive(cereal_as_pod(*this));
    }

private:
    size_t index_ : 61;
    unsigned char mask_ : 3;
};

namespace std {
template <>
struct hash<CachedAACursor> {
    std::size_t operator()(const CachedAACursor &c) const { return std::hash<size_t>()(c.to_size_t()); }
};

}  // namespace std

inline std::ostream &operator<<(std::ostream &os, const CachedAACursor &c) {
    if (c.is_empty()) {
        return os << "(@)";
    } else {
        return os << "(" << c.index() << ":" << std::bitset<3>(c.mask()) << ")";
    }
}

class CachedAACursorContext {
public:
    using Index = uint32_t;

    template <typename Cursor>
    auto UnpackCursor(const CachedAACursor &cursor, const std::vector<Cursor> &cursors) {
        const auto &triplet = triplets_[cursor.index()];
        return AAGraphCursor<Cursor>(cursors[triplet[0]], cursors[triplet[1]], cursors[triplet[2]], cursor.mask());
    }

    template <typename Cursor>
    std::vector<AAGraphCursor<Cursor>> UnpackPath(const std::vector<CachedAACursor> &path, const std::vector<Cursor> &cursors) {
        std::vector<AAGraphCursor<Cursor>> result;
        result.reserve(path.size());
        for (const auto &cursor : path) {
            result.push_back(UnpackCursor(cursor, cursors));
        }
        return result;
    }

    std::vector<CachedAACursor> Cursors() const {
        std::vector<CachedAACursor> result;
        size_t size = letters_.size();
        result.reserve(size);
        for (size_t i = 0; i < size; ++i) {
            result.push_back(CachedAACursor(i, 0b111));
        }
        return result;
    }

    template <typename Cursor>
    CachedAACursorContext(const std::vector<Cursor> &cursors, typename Cursor::Context context) {
        VERIFY(std::numeric_limits<Index>::max() > cursors.size());
        std::unordered_map<Cursor, Index> cursor2index;
        for (Index i = 0; i < cursors.size(); ++i) {
            cursor2index[cursors[i]] = i;
        }

        std::unordered_map<AAGraphCursor<Cursor>, size_t> aa_cursor2index;
        auto get_or_create = [&](const auto &aa_cursor) -> size_t {
            VERIFY(!aa_cursor.is_empty());
            auto triplet_form = aa_cursor.triplet_form();
            VERIFY(!triplet_form.is_empty());
            auto triplet = aa_cursor.triplet_cursors();
            if (!aa_cursor2index.count(triplet_form)) {
                aa_cursor2index[triplet_form] = triplets_.size();
                triplets_.push_back({cursor2index[triplet[0]], cursor2index[triplet[1]], cursor2index[triplet[2]]});
                letters_.push_back(triplet_form.letter(context));
                return triplets_.size() - 1;
            } else {
                return aa_cursor2index[triplet_form];
            }
        };

        auto get = [&](const auto &aa_cursor) -> size_t {
            auto triplet_form = aa_cursor.triplet_form();
            auto triplet = aa_cursor.triplet_cursors();
            auto it = aa_cursor2index.find(triplet_form);
            VERIFY(it != aa_cursor2index.cend());
            return it->second;
        };

        auto aa_cursors = make_aa_cursors(cursors, context);
        for (const auto &cursor : aa_cursors) {
            get_or_create(cursor);
            for (const auto &n : cursor.next(context)) {
                get_or_create(n);
            }
            for (const auto &p : cursor.prev(context)) {
                get_or_create(p);
            }
            for (const auto &n : cursor.next_frame_shift(context)) {
                get_or_create(n);
            }
        }

        nexts_.resize(triplets_.size());
        prevs_.resize(triplets_.size());
        nexts_frame_shift_.resize(triplets_.size());
        for (size_t i = 0; i < triplets_.size(); ++i) {
            auto cc = CachedAACursor(i, 0b111);
            auto cursor = UnpackCursor(cc, cursors);
            for (const auto &c : cursor.next(context)) {
                nexts_[i].emplace_back(get(c), c.mask());
            }
            for (const auto &c : cursor.prev(context)) {
                prevs_[i].emplace_back(get(c), c.mask());
            }
            for (const auto &c : cursor.next_frame_shift(context)) {
                nexts_frame_shift_[i].emplace_back(get(c), c.mask());
            }
        }
    }

    friend class CachedAACursor;

    template <class Archive>
    void serialize(Archive &archive) {
        archive(triplets_, letters_, nexts_, prevs_, nexts_frame_shift_);
    }

private:
    std::vector<std::array<Index, 3>> triplets_;
    std::vector<char> letters_;
    std::vector<std::vector<CachedAACursor>> nexts_;
    std::vector<std::vector<CachedAACursor>> prevs_;
    std::vector<std::vector<CachedAACursor>> nexts_frame_shift_;
};

// FIXME add cpp

inline char CachedAACursor::letter(CachedAACursor::Context context) const {
    switch (mask_) {
        case 0b111:
            return context->letters_[index_];
        case 0b110:
            return '=';
        case 0b100:
            return '-';
        default:
            return '*';
    }
}

inline std::vector<size_t> CachedAACursor::nucl_cursor_indices(Context context) const {
    if (is_empty()) return {};
    std::vector<size_t> result;
    const auto &triplet = context->triplets_[index_];
    if (mask_ & 1) result.push_back(triplet[0]);
    if ((mask_ >> 1) & 1) result.push_back(triplet[1]);
    if ((mask_ >> 2) & 1) result.push_back(triplet[2]);
    return result;
}

inline const std::vector<CachedAACursor> &CachedAACursor::next(CachedAACursor::Context context) const {
    VERIFY(!is_empty());
    return context->nexts_[index_];
}
inline const std::vector<CachedAACursor> &CachedAACursor::prev(CachedAACursor::Context context) const {
    VERIFY(!is_empty());
    return context->prevs_[index_];
}
inline const std::vector<CachedAACursor> &CachedAACursor::next_frame_shift(CachedAACursor::Context context) const {
    VERIFY(!is_empty());
    return context->nexts_frame_shift_[index_];
}

// vim: set ts=4 sw=4 et :
