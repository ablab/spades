//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <type_traits>
#include <iostream>
#include "llvm/Support/LEB128.h"

namespace io {

namespace binary {

enum EncodingType {
    Unknown,
    Raw,
    LEB128
};

static const char LEB_BUF_SIZE = 128 / 7 + 1;

template<typename T>
static constexpr EncodingType GetEncodingType() {
    return std::is_unsigned<T>::value ? LEB128 :
           (std::is_pod<T>::value && !std::is_array<T>::value) ? Raw :
           Unknown;
}

template<typename T>
typename std::enable_if<GetEncodingType<T>() == Raw>::type BinWrite(std::ostream &str, const T &value) {
    str.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

template<typename T>
typename std::enable_if<GetEncodingType<T>() == Raw, T>::type BinRead(std::istream &str, T &value) {
    str.read(reinterpret_cast<char *>(&value), sizeof(T));
}

template<typename T>
typename std::enable_if<GetEncodingType<T>() == LEB128>::type BinWrite(std::ostream &str, const T &value) {
    uint8_t buf[LEB_BUF_SIZE];
    auto count = llvm::encodeULEB128(value, buf);
    str.write(reinterpret_cast<const char *>(buf), count);
}

template<typename T>
typename std::enable_if<GetEncodingType<T>() == LEB128>::type BinRead(std::istream &str, T &value) {
    uint8_t buf[LEB_BUF_SIZE];
    auto pos = reinterpret_cast<char *>(buf);
    char count = 0;
    do {
        str.read(pos, 1);
        VERIFY_MSG(++count < LEB_BUF_SIZE, "Malformed LEB128 sequence");
    } while (*(pos++) & 0x80);
    value = static_cast<T>(llvm::decodeULEB128(buf));
}

template<typename T>
auto BinWrite(std::ostream &str, const T &value) -> decltype(std::declval<const T>().BinWrite(str), void()) {
    value.BinWrite(str);
}

template<typename T>
auto BinRead(std::istream &str, T &value) -> decltype(std::declval<T>().BinRead(str), void()) {
    value.BinRead(str);
}

//Specialization for sized containers
/*
template<typename T>
decltype(std::declval<const T>().begin(), std::declval<const T>().end(), std::declval<const T>().size(), void())
BinWrite(std::ostream &str, const T &value) {
    BinWrite(str, (size_t)value.size());
    for (const auto &i : value)
        BinWrite(str, i);
}
 */

//Ad-hoc overloads
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

inline void BinWrite(std::ostream &str, const std::string &value) {
    BinWrite(str, (size_t)value.length());
    str.write(value.data(), value.length());
}

inline void BinRead(std::istream &str, std::string &value) {
    size_t size;
    BinRead(str, size);
    value.reserve(size);
    str.read(const_cast<char *>(value.data()), value.length());
}

//
template<typename T>
T BinRead(std::istream &str) {
    T result;
    BinRead(str, result);
    return result;
}

}

}
