//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <type_traits>
#include <istream>
#include "llvm/Support/LEB128.h"

namespace io {

namespace binary {

enum class Encoding {
    Unknown,
    Raw,
    SLEB128,
    ULEB128
};

static const char LEB_BUF_SIZE = 128 / 7 + 1;

template<typename T>
static constexpr Encoding GetEncoding() {
    return !std::is_pod<T>::value ? Encoding::Unknown :
               std::is_unsigned<T>::value && sizeof(T) > 1 ? Encoding::ULEB128 :
               std::is_signed<T>::value && sizeof(T) > 1 ? Encoding::SLEB128 :
               Encoding::Raw;
}

template<typename T>
typename std::enable_if_t<GetEncoding<T>() == Encoding::Raw> BinWrite(std::ostream &str, const T &value) {
    str.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

template<typename T>
typename std::enable_if_t<GetEncoding<T>() == Encoding::Raw> BinRead(std::istream &str, T &value) {
    str.read(reinterpret_cast<char *>(&value), sizeof(T));
}

template<typename T>
typename std::enable_if_t<GetEncoding<T>() == Encoding::ULEB128> BinWrite(std::ostream &str, const T &value) {
    uint8_t buf[LEB_BUF_SIZE];
    auto count = llvm::encodeULEB128(value, buf);
    str.write(reinterpret_cast<const char *>(buf), count);
}

template<typename T>
typename std::enable_if_t<GetEncoding<T>() == Encoding::ULEB128> BinRead(std::istream &str, T &value) {
    uint8_t buf[LEB_BUF_SIZE];
    auto pos = reinterpret_cast<char *>(buf);
    char count = 0;
    do {
        str.read(pos, 1);
        VERIFY_MSG(++count < LEB_BUF_SIZE, "Malformed LEB128 sequence");
    } while (*(pos++) & 0x80);
    value = static_cast<T>(llvm::decodeULEB128(buf));
}

// Implementation for classes that have corresponding BinWrite/BinRead methods.
template<typename T>
typename std::enable_if_t<GetEncoding<T>() == Encoding::SLEB128> BinWrite(std::ostream &str, const T &value) {
    uint8_t buf[LEB_BUF_SIZE];
    auto count = llvm::encodeSLEB128(value, buf);
    str.write(reinterpret_cast<const char *>(buf), count);
}

template<typename T>
typename std::enable_if_t<GetEncoding<T>() == Encoding::SLEB128> BinRead(std::istream &str, T &value) {
    uint8_t buf[LEB_BUF_SIZE];
    auto pos = reinterpret_cast<char *>(buf);
    char count = 0;
    do {
        str.read(pos, 1);
        VERIFY_MSG(++count < LEB_BUF_SIZE, "Malformed LEB128 sequence");
    } while (*(pos++) & 0x80);
    value = static_cast<T>(llvm::decodeSLEB128(buf));
}

template<typename T>
auto BinWrite(std::ostream &str, const T &value) -> decltype(std::declval<const T>().BinWrite(str), void()) {
    value.BinWrite(str);
}

template<typename T>
auto BinRead(std::istream &str, T &value) -> decltype(std::declval<T>().BinRead(str), void()) {
    value.BinRead(str);
}

//---- Ad-hoc overloads ------------------------------------------------------------------------------------------------
template <typename T, typename... Ts>
void BinWrite(std::ostream &os, const T &v, const Ts &... vs);
template <typename T, typename... Ts>
void BinRead(std::istream &is, T &v, Ts &... vs);

inline void BinWrite(std::ostream &str, const std::string &value);
inline void BinRead(std::istream &str, std::string &value);

template <typename T>
std::enable_if_t<!std::is_same<T, bool>::value> BinRead(std::istream &is, std::vector<T> &v);
template <typename T>
std::enable_if_t<!std::is_same<T, bool>::value> BinWrite(std::ostream &os, const std::vector<T> &v);

template <typename T1, typename T2>
void BinWrite(std::ostream &os, const std::pair<T1, T2> &p);
template <typename T1, typename T2>
void BinRead(std::istream &is, std::pair<T1, T2> &p);

template <typename K, typename V, typename... Args>
void BinWrite(std::ostream &os, const std::map<K, V, Args...> &m);
template <typename K, typename V, typename... Args>
void BinRead(std::istream &is, std::map<K, V, Args...> &m);

template <typename K, typename V, typename... Args>
void BinWrite(std::ostream &os, const std::unordered_map<K, V, Args...> &m);
template <typename K, typename V, typename... Args>
void BinRead(std::istream &is, std::unordered_map<K, V, Args...> &m);

// Arrays

template<typename T, size_t N>
void BinWrite(std::ostream &str, const T (&value)[N]) {
    for (size_t i = 0; i < N; ++i)
        BinWrite(str, value[i]);
}

template<typename T, size_t N>
void BinRead(std::istream &str, T (&value)[N]) {
    for (size_t i = 0; i < N; ++i)
        BinRead(str, value[i]);
}

// Strings

inline void BinWrite(std::ostream &str, const std::string &value) {
    BinWrite(str, static_cast<size_t>(value.length()));
    str.write(value.data(), value.length());
}

inline void BinRead(std::istream &str, std::string &value) {
    size_t size;
    BinRead(str, size);
    value.resize(size);
    str.read(const_cast<char *>(value.data()), value.length());
}

// Vectors

template <typename T>
std::enable_if_t<!std::is_same<T, bool>::value> BinWrite(std::ostream &os, const std::vector<T> &v) {
    BinWrite(os, v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        BinWrite(os, v[i]);
    }
}

template <typename T>
std::enable_if_t<!std::is_same<T, bool>::value> BinRead(std::istream &is, std::vector<T> &v) {
    size_t size;
    BinRead(is, size);
    v.resize(size);
    for (size_t i = 0; i < size; ++i) {
        BinRead(is, v[i]);
    }
}

// Tuples

// TODO implement for std::tuple
template <typename T1, typename T2>
void BinWrite(std::ostream &os, const std::pair<T1, T2> &p) {
    BinWrite(os, p.first, p.second);
}

template <typename T1, typename T2>
void BinRead(std::istream &is, std::pair<T1, T2> &p) {
    BinRead(is, p.first, p.second);
}

// Maps

template <typename K, typename V, typename... Args>
void BinWrite(std::ostream &os, const std::map<K, V, Args...> &m) {
    size_t size = m.size();
    BinWrite(os, size);
    for (const auto &kv : m) {
        BinWrite(os, kv.first, kv.second);
    }
}

template <typename K, typename V, typename... Args>
void BinRead(std::istream &is, std::map<K, V, Args...> &m) {
    m.clear();
    size_t size;
    BinRead(is, size);
    for (size_t i = 0; i < size; ++i) {
        K k;
        V v;
        BinRead(is, k, v);
        m.insert({std::move(k), std::move(v)});
    }
}

template <typename K, typename V, typename... Args>
void BinWrite(std::ostream &os, const std::unordered_map<K, V, Args...> &m) {
    size_t size = m.size();
    BinWrite(os, size);
    for (const auto &kv : m) {
        BinWrite(os, kv.first, kv.second);
    }
}

template <typename K, typename V, typename... Args>
void BinRead(std::istream &is, std::unordered_map<K, V, Args...> &m) {
    m.clear();
    size_t size;
    BinRead(is, size);
    for (size_t i = 0; i < size; ++i) {
        K k;
        V v;
        BinRead(is, k, v);
        m.insert({std::move(k), std::move(v)});
    }
}

//---- Helpers ---------------------------------------------------------------------------------------------------------

/**
 * @brief A value-returning variant of BinRead to save some LOCs.
 */
template<typename T>
T BinRead(std::istream &str) {
    T result;
    BinRead(str, result);
    return result;
}

/**
 * @brief Multi-value wrappers to save even more LOCs.
 */
template <typename T, typename... Ts>
void BinWrite(std::ostream &os, const T &v, const Ts &... vs) {
    BinWrite(os, v);
    BinWrite(os, vs...);
}

template <typename T, typename... Ts>
void BinRead(std::istream &is, T &v, Ts &... vs) {
    BinRead(is, v);
    BinRead(is, vs...);
}

} // namespace binary

} // namespace io
