//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/verify.hpp"

#include <algorithm>
#include <memory>
#include <map>
#include <string>
#include <set>
#include <vector>
#include <ostream>

namespace utils {

template<class Container>
std::string ContainerToString(const Container &c) {
    std::ostringstream ss;
    ss << "Size " << c.size() << ": [";
    for (const auto &el : c)
        ss << el << ", ";
    ss << "]";
    return ss.str();
}

template<typename T>
inline const std::pair<T, T> ReversePair(std::pair<T, T> ep) {
    return std::pair<T, T>(ep.second, ep.first);
}

template<class ContainerT1, class ContainerT2>
void push_back_all(ContainerT1 &target, const ContainerT2 &to_insert) {
    target.insert(target.end(), to_insert.begin(), to_insert.end());
}

template<class ContainerT1, class ContainerT2>
void insert_all(ContainerT1 &target, const ContainerT2 &to_insert) {
    target.insert(to_insert.begin(), to_insert.end());
}

template<class MapT>
std::set<typename MapT::key_type> key_set(const MapT &m) {
    std::set<typename MapT::key_type> answer;
    for (auto it = m.begin(); it != m.end(); ++it) {
        answer.insert(it->first);
    }
    return answer;
}

template<class MapT>
std::set<typename MapT::mapped_type> value_set(const MapT &m) {
    std::set<typename MapT::mapped_type> answer;
    for (auto it = m.begin(); it != m.end(); ++it) {
        answer.insert(it->second);
    }
    return answer;
}

template<class MapT>
const typename MapT::mapped_type &get(const MapT &from, const typename MapT::key_type &key) {
    auto it = from.find(key);
    VERIFY(it != from.end());
    return it->second;
}

template<class MapT>
typename MapT::mapped_type &get(MapT &from, const typename MapT::key_type &key) {
    auto it = from.find(key);
    VERIFY(it != from.end());
    return it->second;
}

template<class MMapT>
const std::vector<typename MMapT::mapped_type> get_all(const MMapT &from, const typename MMapT::key_type &key) {
    std::vector<typename MMapT::mapped_type> answer;
    for (auto it = from.lower_bound(key); it != from.upper_bound(key); ++it) {
        answer.push_back(it->second);
    }
    return answer;
}

template<class Container, class F>
std::string join(const Container &c,
                 const std::string &delim = ", ",
                 F str_f = [](typename Container::value_type t) { return std::to_string(t); }) {
    std::stringstream ss;
    std::string d = "";
    for (const auto &item : c) {
        ss << d << str_f(item);
        d = delim;
    }
    return ss.str();
}

template<typename Out>
void split(const std::string_view s, const std::string_view delims, Out result,
           bool compress = false) {
    size_t last_pos = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (delims.find(s[i]) == delims.npos)
            continue;
        auto item = s.substr(last_pos, i - last_pos);
        if (!compress || item.size() > 0)
            *result++ = item;
        last_pos = i + 1;
    }

    if (last_pos != s.size()) {
        auto item = s.substr(last_pos);
        if (!compress || item.size() > 0)
            *result++ = item;
    }
}

template<typename Out, class Mapper>
void split(const std::string_view s, const std::string_view delims,
           Out result, Mapper f) {
    std::string item;
    size_t last_pos = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (delims.find(s[i]) == delims.npos)
            continue;
        *result++ = f(s.substr(last_pos, i - last_pos));
        last_pos = i + 1;
    }

    if (last_pos != s.size()) {
        *result++ = f(s.substr(last_pos));
    }
}

static inline std::vector<std::string_view>
split(const std::string &s, const std::string_view delim,
      bool compress = false) {
    std::vector<std::string_view> elems;
    split(s, delim, std::back_inserter(elems), compress);
    return elems;
}

static inline bool starts_with(const std::string_view s, const std::string_view p) {
    if (s.size() < p.size())
        return false;

    return (s.compare(0, p.size(), p) == 0);
}

static inline bool ends_with(const std::string_view s, const std::string_view p) {
    if (s.size() < p.size())
        return false;

    return (s.compare(s.size() - p.size(), p.size(), p) == 0);
}

static inline std::string str_tolower(std::string_view s) {
    std::string result(s);

    std::transform(s.begin(), s.end(), result.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return result;
}

static inline std::string str_toupper(std::string_view s) {
    std::string result(s);

    std::transform(s.begin(), s.end(), result.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    return result;
}

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(),
            std::find_if(s.begin(), s.end(),
                         [](unsigned char ch) {
                             return !std::isspace(ch);
                         }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [](unsigned char ch) {
                             return !std::isspace(ch);
                         }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}

}

namespace std {
template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> const &pair) {
    return os << "(" << pair.first << ", " << pair.second << ")";
}

template<class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
    os << "[";
    std::string delim = "";
    for (auto it = v.begin(); it != v.end(); ++it) {
        os << delim << *it;
        delim = ", ";
    }
//     std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, ", "));
    os << "]";
    return os;
}

template<class T>
std::ostream &operator<<(std::ostream &os, const std::set<T> &set) {
    os << "{";
    bool delim = false;
    for (const auto &i : set) {
        if (delim) os << ", ";
        os << i;
        delim = true;
    }
    os << "}";
    return os;
}

template<typename K, typename V>
std::ostream &operator<<(std::ostream &os, const std::map<K, V> &map) {
    os << "{";
    bool delim = false;
    for (const auto &i : map) {
        if (delim) os << ", ";
        os << i.first << ": " << i.second;
        delim = true;
    }
    os << "}";
    return os;
}

}
