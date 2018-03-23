#pragma once
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <vector>

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v);
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::unordered_set<T> &v);

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  os << "[ ";
  for (const auto &e : v) {
    os << e << " ";
  }
  os << "]";
  return os;
}

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

inline size_t hash_size_t_pair(size_t s0, size_t s1) {
  s1 ^= s1 << 23;  // a
  return (s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26)) + s0;
}

namespace std {
template <typename T>
struct hash<std::vector<T>> {
    std::size_t operator()(const std::vector<T> &v) const {
        size_t result = 0xDEADBEEF;
        for (const auto &entry :v)
            result = hash_size_t_pair(result, std::hash<T>()(entry));
        return result;
    }
};
}
