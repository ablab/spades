#pragma once
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>
#include <sys/stat.h>


#include "utils/logger/log_writers.hpp"
#include "utils/logger/log_writers_thread.hpp"
#include "utils/stl_utils.hpp"
#include <cereal/cereal.hpp>
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

namespace impl {
template <size_t N, typename... Ts>
struct _hash {
    size_t operator()(const std::tuple<Ts...> &t) const {
        return hash_size_t_pair(_hash<N - 1, Ts...>()(t), hash_value(std::get<N - 1>(t)));
    }
};

template <typename... Ts>
struct _hash<0, Ts...> {
    size_t operator()(const std::tuple<Ts...> &) const { return 0xDEADBEEF; }
};
}  // namespace impl

template <typename... Ts>
struct hash<std::tuple<Ts...>> {
    std::size_t operator()(const std::tuple<Ts...> &t) const { return impl::_hash<sizeof...(Ts), Ts...>{}(t); }
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

inline void to_upper_case(std::string &s) {
    for (char &c: s) c = static_cast<char>(std::toupper(c));
}

template <typename IntType, typename Pred>
IntType int_max_binsearch(const Pred &pred, IntType surely_true, IntType surely_false) {
    while (surely_true + 1 < surely_false) {
        IntType mid = (surely_true + surely_false) / 2;
        if (pred(mid)) {
            surely_true = mid;
        } else {
            surely_false = mid;
        }
    }

    return surely_true;
}

inline void create_console_logger(const std::string &filename = "") {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<console_writer_thread>()));
    if (filename != "") {
        lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<file_writer_thread>(filename)));
    }
    attach_logger(lg);
}

template <typename T>
auto cereal_as_pod(T &ref) {
    return cereal::binary_data(&ref, sizeof(ref));
}

template <typename T>
T saturated_inc(T v, const T max = std::numeric_limits<T>::max()) {
    return v == max ? v : v + 1;
}

template <typename T>
T saturated_sum(T v1, T v2) {
    static_assert(std::is_unsigned<T>::value, "T should be unsigned");
    T result = v1 + v2;
    if (result < v1) {  // overflow
        result = std::numeric_limits<T>::max();
    }
    return result;
}

template <typename T>
constexpr bool is_power_of_two_or_zero(T v) {
    return (v & (v - 1)) == 0;
}

inline bool double_equal(double a, double b) {
    return std::fabs(a - b) < 10 * std::numeric_limits<double>::epsilon() || (a == std::numeric_limits<double>::infinity() && b == std::numeric_limits<double>::infinity());
}
