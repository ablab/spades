//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <limits>
#include <unordered_map>
#include <vector>
#include <cstdint>

#include "common/utils/verify.hpp"

using Index = uint32_t;

class CachedCursorContext;

class CachedCursor {
public:
    using Context = CachedCursorContext;
    CachedCursor(Index index = Index(-1)) : index_{index} {}
    bool is_empty() const { return index_ == Index(-1); }
    bool operator==(const CachedCursor &other) const { return index_ == other.index_; }
    const std::vector<CachedCursor> &next(const void *context) const;
    const std::vector<CachedCursor> &prev(const void *context) const;
    char letter(const void *context) const;
    Index index() const { return index_; }

private:
    Index index_;
};

namespace std {
template <>
struct hash<CachedCursor> {
    std::size_t operator()(const CachedCursor &c) const { return std::hash<size_t>()(c.index()); }
};

}  // namespace std

inline std::ostream &operator<<(std::ostream &os, const CachedCursor &c) {
    if (c.is_empty()) {
        return os << "(@)";
    } else {
        return os << "(" << c.index() << ")";
    }
}

class CachedCursorContext {
public:
    template <typename Cursor>
    static std::vector<Cursor> UnpackPath(const std::vector<CachedCursor> &path, const std::vector<Cursor> &cursors) {
        std::vector<Cursor> result;
        result.reserve(path.size());
        for (const auto &cursor : path) {
            result.push_back(cursors[cursor.index()]);
        }
        return result;
    }

    std::vector<CachedCursor> Cursors() const {
        std::vector<CachedCursor> result;
        size_t size = letters_.size();
        VERIFY(std::numeric_limits<Index>::max() > size);
        result.reserve(size);
        for (Index i = 0; i < size; ++i) {
            result.emplace_back(i);
        }
        return result;
    }

    template <typename Cursor>
    CachedCursorContext(const std::vector<Cursor> &cursors, const void *context) {
        VERIFY(std::numeric_limits<Index>::max() > cursors.size());
        std::unordered_map<Cursor, Index> cursor2index;
        for (size_t i = 0; i < cursors.size(); ++i) {
            cursor2index[cursors[i]] = i;
        }

        letters_.resize(cursors.size());
        nexts_.resize(cursors.size());
        prevs_.resize(cursors.size());
        for (size_t i = 0; i < cursors.size(); ++i) {
            const auto &cursor = cursors[i];
            letters_[i] = cursor.letter(context);
            for (const auto &c : cursor.next(context)) {
                nexts_[i].push_back(cursor2index[c]);
                nexts_[i].shrink_to_fit();
            }
            // FIXME prev is not implemented for aa-cursor
            // for (const auto &c : cursor.prev(context)) {
            //     prevs_[i].push_back(cursor2index[c]);
            //     prevs_[i].shrink_to_fit();
            // }
        }
    }

    friend class CachedCursor;

private:
    std::vector<char> letters_;
    std::vector<std::vector<CachedCursor>> nexts_;
    std::vector<std::vector<CachedCursor>> prevs_;
};

// FIXME add cpp

inline char CachedCursor::letter(const void *context) const {
    return static_cast<const CachedCursorContext *>(context)->letters_[index_];
}


inline const std::vector<CachedCursor> &CachedCursor::next(const void *context) const {
    return static_cast<const CachedCursorContext *>(context)->nexts_[index_];
}

// It is not required for find_best_path
// inline const std::vector<CachedCursor> &CachedCursor::prev(const void *context) const {
//     return static_cast<const CachedCursorContext *>(context)->prevs_[index_];
// }
