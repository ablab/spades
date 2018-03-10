#pragma once
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <vector>

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &p);
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

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &p) {
  return os << "(" << p.first << "," << p.second << ")";
}

template <typename T>
void remove_duplicates(std::vector<T> &v) {
  // Remove duplicated items
  std::sort(v.begin(), v.end());
  auto it = std::unique_copy(v.cbegin(), v.cend(), v.begin());
  v.resize(std::distance(v.begin(), it));
}
