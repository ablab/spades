#pragma once
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <unordered_set>
#include <utility>
#include <vector>
#include <sys/stat.h>


#include "utils/stl_utils.hpp"
using namespace utils;

// template <typename T>
// std::ostream &operator<<(std::ostream &os, const std::unordered_set<T> &v);
//
// template <typename T>
// std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
//   os << "[ ";
//   for (const auto &e : v) {
//     os << e << " ";
//   }
//   os << "]";
//   return os;
// }
//
// template <typename T>
// std::ostream &operator<<(std::ostream &os, const std::vector<T> &v);

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::unordered_set<T> &v) {
    os << "{ ";
    for (const auto &e : v) {
        os << e << " ";
    }
    os << "}";
    return os;
}

template <typename T>
void remove_duplicates(std::vector<T> &v) {
  // Remove duplicated items
  std::sort(v.begin(), v.end());
  auto it = std::unique_copy(v.cbegin(), v.cend(), v.begin());
  v.resize(std::distance(v.begin(), it));
}

template <class InputIt, class OutputIt, class Key>
OutputIt unique_copy_by(InputIt first, InputIt last, OutputIt d_first, Key key) {
    return std::unique_copy(first, last, d_first,
                            [&key](const auto &e1, const auto &e2) -> bool { return key(e1) == key(e2); });
}

inline size_t hash_size_t_pair(size_t s0, size_t s1) {
  s1 ^= s1 << 23;  // a
  return (s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26)) + s0;
}

template <typename T>
auto hash_value(const T& v) {
    return std::hash<T>{}(v);
}

template <typename T>
std::string int_to_hex(const T &i) {
  std::stringstream ss;
  ss << std::setfill('0') << std::setw(sizeof(T) * 2) << std::hex << i;
  return ss.str();
}

namespace std {
template <typename T>
struct hash<std::vector<T>> {
    std::size_t operator()(const std::vector<T> &v) const {
        size_t result = 0xDEADBEEF;
        for (const auto &entry : v) result = hash_size_t_pair(result, std::hash<T>()(entry));
        return result;
    }
};
}  // namespace std

template <typename Iter, typename Key>
void sort_by(Iter b, Iter e, const Key &key) {
    std::sort(b, e, [&key](const auto &a, const auto &b) -> bool { return key(a) < key(b); });
}

template <typename Iter, typename Key>
void stable_sort_by(Iter b, Iter e, const Key &key) {
    std::stable_sort(b, e, [&key](const auto &a, const auto &b) -> bool { return key(a) < key(b); });
}

template <typename Range, typename Sep>
std::string join(const Range &range, const Sep &sep) {
    std::stringstream ss;
    size_t inserted = 0;
    for (const auto &e : range) {
        if (inserted > 0) {
            ss << sep;
        }
        ss << e;
        ++inserted;
    }

    return ss.str();
}

inline bool ends_with(const std::string &s, const std::string &p) {
    if (s.size() < p.size())
        return false;

    return (s.compare(s.size() - p.size(), p.size(), p) == 0);
}

inline bool is_dir(const std::string &path) {
    struct stat sb;
    return (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode));
}

inline std::string compress_alignment(const std::string &alignment, bool x_as_m = false) {
    size_t count = 0;
    char prev_c = '\0';
    std::string result;
    for (size_t i = 0; i <= alignment.size(); ++i) {
        char c = alignment[i];
        if (x_as_m && c == 'X') {
            c = 'M';
        }
        if (c == prev_c) {
            ++count;
        } else {
            if (prev_c != '\0') {
                result += std::to_string(count);
                if (prev_c == '-') {
                    prev_c = 'D';
                }
                result += prev_c;
            }
            count = 1;
        }
        prev_c = c;
    }

    return result;
}
