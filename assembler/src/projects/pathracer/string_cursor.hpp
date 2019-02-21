//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include <vector>

class StringCursor {
public:
    StringCursor(size_t pos = -1) : pos_{pos} {}
    char letter(const void *context) const {
        const std::string &s = *static_cast<const std::string *>(context);
        return s[pos_];
    }

    bool is_empty() const { return pos_ == size_t(-1); }
    bool operator==(const StringCursor &other) const { return pos_ == other.pos_; }

    StringCursor(const StringCursor &) = default;
    StringCursor(StringCursor &&) = default;
    StringCursor &operator=(const StringCursor &) = default;
    StringCursor &operator=(StringCursor &&) = default;
    ~StringCursor() noexcept = default;

    std::vector<StringCursor> prev(const void *) const {
        if (pos_ == 0) {
            return {};
        } else {
            return {StringCursor(pos_ - 1)};
        }
    }

    std::vector<StringCursor> next(const void *context) const {
        const std::string &s = *static_cast<const std::string *>(context);
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

inline std::ostream &operator<<(std::ostream &os, const StringCursor &c) {
    if (c.is_empty()) {
        return os << "(@)";
    } else {
        return os << "(" << c.position() << ")";
    }
}

