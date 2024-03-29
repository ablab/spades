//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include <vector>

class StringCursor {
public:
    using Context = const std::string*;

    StringCursor(size_t pos = -1) : pos_{pos} {}
    char letter(Context context) const {
        return (*context)[pos_];
    }

    bool is_empty() const { return pos_ == size_t(-1); }
    size_t edge() const { return 0; }  // FIXME this method is required for AAGraphCursor wrapper only
    bool operator==(const StringCursor &other) const { return pos_ == other.pos_; }

    StringCursor(const StringCursor &) = default;
    StringCursor(StringCursor &&) = default;
    StringCursor &operator=(const StringCursor &) = default;
    StringCursor &operator=(StringCursor &&) = default;
    ~StringCursor() noexcept = default;

    std::vector<StringCursor> prev(Context) const {
        if (pos_ == 0) {
            return {};
        } else {
            return {StringCursor(pos_ - 1)};
        }
    }

    std::vector<StringCursor> next(Context context) const {
        const std::string &s = *context;
        if (pos_ == s.length() - 1) {
            return {};
        } else {
            return {StringCursor(pos_ + 1)};
        }
    }

    size_t position() const { return pos_; }

private:
    size_t pos_;
};

namespace std {
template <>
struct hash<StringCursor> {
    std::size_t operator()(const StringCursor &c) const { return std::hash<size_t>()(c.position()); }
};
}  // namespace std

inline bool operator<(const StringCursor &c1, const StringCursor &c2) { return c1.position() < c2.position(); }

inline std::ostream &operator<<(std::ostream &os, const StringCursor &c) {
    if (c.is_empty()) {
        return os << "(@)";
    } else {
        return os << "(" << c.position() << ")";
    }
}
