//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

#include "utils/verify.hpp"
#include "io/binary/access.hpp"

namespace io {

namespace binary {

namespace impl {

enum class Encoding {
    Unknown,
    Raw,
    SLEB128,
    ULEB128
};

// TODO implement special wrapper for pointers and disable direct serialization
template<typename T>
static constexpr Encoding GetEncoding() {
    return !std::is_pod<T>::value ? Encoding::Unknown :
           std::is_floating_point<T>::value || sizeof(T) == 1 || std::is_pointer<T>::value ? Encoding::Raw : //Required because is_signed is true for floats
           std::is_unsigned<T>::value ? Encoding::ULEB128 :
           std::is_signed<T>::value ? Encoding::SLEB128 :
           Encoding::Unknown;
}

template <typename T, typename Enable = void>
class Serializer;

// Motivated by declval, the same stuff but returns lvalue reference
template <class T>
typename std::add_lvalue_reference<T>::type declref() noexcept;

template <class T>
typename std::add_lvalue_reference<const T>::type declcref() noexcept;

namespace impl {
template <typename T, typename Enable = void>
constexpr bool is_serializable = false;

template <typename T>
constexpr bool is_serializable<T, decltype(Serializer<T>::Write(declref<std::ostream>(), declcref<T>()),
                                           Serializer<T>::Read(declref<std::istream>(), declref<T>()), void())> = true;

template <typename... Ts>
constexpr std::enable_if_t<sizeof...(Ts) == 0, bool> is_all_serializable_f() {
    return true;
}

template <typename T, typename... Ts>
constexpr bool is_all_serializable_f() {
    return is_serializable<T> && is_all_serializable_f<Ts...>();
}

}  // namespace impl

template <typename... Ts>
constexpr bool is_serializable = impl::is_all_serializable_f<Ts...>();

template <typename T>
std::enable_if_t<is_serializable<T>> BinWrite(std::ostream &os, const T &v) {
    Serializer<T>::Write(os, v);
}

template <typename T>
std::enable_if_t<is_serializable<T>> BinRead(std::istream &is, T &v) {
    Serializer<T>::Read(is, v);
}

// Multi-value wrappers to save even more LOCs.
template <typename T, typename... Ts>
std::enable_if_t<is_serializable<T, Ts...> && sizeof...(Ts) != 0> BinWrite(std::ostream &os, const T &v, const Ts &... vs) {
    BinWrite(os, v);
    BinWrite(os, vs...);
}

template <typename T, typename... Ts>
std::enable_if_t<is_serializable<T, Ts...> && sizeof...(Ts) != 0> BinRead(std::istream &is, T &v, Ts &... vs) {
    BinRead(is, v);
    BinRead(is, vs...);
}

template <typename T>
class Serializer<T, std::enable_if_t<GetEncoding<T>() == Encoding::Raw>> {
public:
    static void Write(std::ostream &os, const T &value) {
        os.write(reinterpret_cast<const char *>(&value), sizeof(T));
    }

    static void Read(std::istream &is, T &value) {
        is.read(reinterpret_cast<char *>(&value), sizeof(T));
    }
};

inline static void encodeULEB128(uint64_t Value, std::ostream &os) {
    static const char LEB_BUF_SIZE = 64 / 7 + 1;
    uint8_t buf[LEB_BUF_SIZE];
    unsigned count = 0;
    do {
        uint8_t Byte = Value & 0x7f;
        Value >>= 7;
        if (Value != 0)
            Byte |= 0x80; // Mark this byte to show that more bytes will follow.
        buf[count] = Byte;
        ++count;
    } while (Value != 0);
    os.write(reinterpret_cast<const char *>(buf), count);
}

inline static uint64_t decodeULEB128(std::istream &is) {
    uint64_t Value = 0;
    unsigned Shift = 0;
    for (uint8_t Byte = static_cast<uint8_t>(is.get()); Shift < 70; Byte = static_cast<uint8_t>(is.get())) {
        Value |= uint64_t(Byte & 0x7f) << Shift;
        Shift += 7;
        if ((Byte & 0x80) == 0) break;
    }
    VERIFY(Shift < 70 && is);
    return Value;
}

template <typename T>
class Serializer<T, std::enable_if_t<GetEncoding<T>() == Encoding::ULEB128>> {
public:
    static void Write(std::ostream &os, const T &value) {
        encodeULEB128(value, os);
    }

    static void Read(std::istream &is, T &value) {
        value = static_cast<T>(decodeULEB128(is));
    }
};

inline static void encodeSLEB128(int64_t Value, std::ostream &os) {
    static const char LEB_BUF_SIZE = 64 / 7 + 1;
    uint8_t buf[LEB_BUF_SIZE];
    unsigned count = 0;
    bool More;
    do {
        uint8_t Byte = Value & 0x7f;
        // NOTE: this assumes that this signed shift is an arithmetic right shift.
        Value >>= 7;
        More = !((((Value == 0 ) && ((Byte & 0x40) == 0)) ||
                  ((Value == -1) && ((Byte & 0x40) != 0))));
        if (More)
            Byte |= 0x80; // Mark this byte to show that more bytes will follow.
        buf[count] = Byte;
        ++count;
    } while (More);
    os.write(reinterpret_cast<const char *>(buf), count);
}

inline static int64_t decodeSLEB128(std::istream &is) {
    int64_t Value = 0;
    unsigned Shift = 0;
    uint8_t Byte = static_cast<uint8_t>(is.get());
    for (;; Byte = static_cast<uint8_t>(is.get())) {
        Value |= (int64_t(Byte & 0x7f) << Shift);
        Shift += 7;
        if ((Byte & 0x80) == 0) break;
    }
    VERIFY(Shift < 70 && is);
    // Sign extend negative numbers.
    if (Byte & 0x40)
        Value |= (-1ULL) << Shift;
    return Value;
}

template <typename T>
class Serializer<T, std::enable_if_t<GetEncoding<T>() == Encoding::SLEB128>> {
public:
    static void Write(std::ostream &os, const T &value) {
        encodeSLEB128(value, os);
    }

    static void Read(std::istream &is, T &value) {
        value = static_cast<T>(decodeSLEB128(is));
    }
};

namespace impl {
template <typename T, typename Enable = void>
constexpr bool has_binread_binwrite = false;

template <typename T>
constexpr bool has_binread_binwrite<T, decltype(declcref<T>().BinWrite(declref<std::ostream>()), declref<T>().BinRead(declref<std::istream>()), void())> = true;
}

template <typename T>
constexpr bool has_binread_binwrite = impl::has_binread_binwrite<T>;

// Implementation for classes having corresponding BinWrite/BinRead methods:
template <typename T>
class Serializer<T, std::enable_if_t<has_binread_binwrite<T>>> {
public:
    static void Write(std::ostream &os, const T &value) {
        value.BinWrite(os);
    }

    static void Read(std::istream &is, T &value) {
        value.BinRead(is);
    }
};

// Implementation for classes having BinArchive method
class StreamArchiveWriter {
public:
    StreamArchiveWriter(std::ostream &os) : os_{os} {}

    template <typename T>
    void operator()(const T &value) {
        BinWrite(os_, value);
    }

    template <typename T, typename... Ts>
    void operator()(const T &v, const Ts &... vs) {
        operator()(v);
        operator()(vs...);
    }

    template <typename T>
    std::enable_if_t<std::is_pod<T>::value> raw(const T &value) {
        os_.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    template <typename T, typename... Ts>
    auto raw(const T &v, const Ts &... vs) -> decltype(raw(v), raw(vs...), void()) {
        raw(v);
        raw(vs...);
    }

    template <typename T>
    std::enable_if_t<std::is_pod<T>::value> raw_array(const T *p, size_t count) {
        os_.write(reinterpret_cast<const char*>(p), count * sizeof(T));
    }

    std::ostream &stream() {
        return os_;
    }

private:
    std::ostream &os_;
};

class StreamArchiveReader {
public:
    StreamArchiveReader(std::istream &is) : is_{is} {}

    template <typename T>
    void operator()(T &value) {
        BinRead(is_, value);
    }

    template <typename T, typename... Ts>
    void operator()(T &v, Ts &... vs) {
        operator()(v);
        operator()(vs...);
    }

    template <typename T>
    std::enable_if_t<std::is_pod<T>::value> raw(T &value) {
        is_.read(reinterpret_cast<char*>(&value), sizeof(T));
    }

    template <typename T, typename... Ts>
    auto raw(const T &v, const Ts &... vs) -> decltype(raw(v), raw(vs...), void()) {
        raw(v);
        raw(vs...);
    }

    template <typename T>
    std::enable_if_t<std::is_pod<T>::value> raw_array(T *p, size_t count) {
        is_.read(reinterpret_cast<char*>(p), count * sizeof(T));
    }

    std::istream &stream() {
        return is_;
    }

    // That would require 'template' keyword on call
    // template <typename T>
    // T get() {
    //     T v;
    //     operator()(v);
    //     return v;
    // }

    template <typename T>
    T get(T) {
        T v;
        operator()(v);
        return v;
    }

private:
    std::istream &is_;
};

namespace impl {
template <typename T, typename Enable = void>
constexpr bool has_binarchive = false;

template <typename T>
constexpr bool has_binarchive<T, decltype(access().BinArchive(declref<T>(), declref<StreamArchiveWriter>()), access().BinArchive(declref<T>(), declref<StreamArchiveReader>()), void())> = true;

template <typename T, typename Enable = void>
constexpr bool has_binarchive_save_load = false;

template <typename T>
constexpr bool has_binarchive_save_load<T, decltype(access().BinArchiveSave(declcref<T>(), declref<StreamArchiveWriter>()), access().BinArchiveLoad(declref<T>(), declref<StreamArchiveReader>()), void())> = true;
}

template <typename T>
constexpr bool has_binarchive = impl::has_binarchive<T>;

template <typename T>
constexpr bool has_binarchive_save_load = impl::has_binarchive_save_load<T>;

template <typename T>
constexpr bool has_inside_serialization = has_binread_binwrite<T> || has_binarchive<T> || has_binarchive_save_load<T>;

template <typename T>
class Serializer<T, std::enable_if_t<has_binarchive<T>>> {
public:
    static void Write(std::ostream &os, const T &value) {
        StreamArchiveWriter saw(os);
        access().BinArchive(const_cast<T&>(value), saw);
    }

    static void Read(std::istream &is, T &value) {
        StreamArchiveReader sar(is);
        access().BinArchive(value, sar);
    }
};

template <typename T>
class Serializer<T, std::enable_if_t<has_binarchive_save_load<T>>> {
public:
    static void Write(std::ostream &os, const T &value) {
        StreamArchiveWriter saw(os);
        access().BinArchiveSave(value, saw);
    }

    static void Read(std::istream &is, T &value) {
        StreamArchiveReader sar(is);
        access().BinArchiveLoad(value, sar);
    }
};

// Arrays
// TODO Disable direct array serialization
template <typename T, size_t N>
class Serializer<T[N], std::enable_if_t<is_serializable<T>>> {
public:
    static void Write(std::ostream &os, const T (&value)[N]) {
        for (size_t i = 0; i < N; ++i)
            BinWrite(os, value[i]);
    }

    static void Read(std::istream &is, T (&value)[N]) {
        for (size_t i = 0; i < N; ++i)
            BinRead(is, value[i]);
    }
};

// Strings
template <>
class Serializer<std::string> {
public:
    static void Write(std::ostream &os, const std::string &value) {
        BinWrite(os, static_cast<size_t>(value.length()));
        os.write(value.data(), value.length());
    }

    static void Read(std::istream &is, std::string &value) {
        size_t size;
        BinRead(is, size);
        value.resize(size);
        is.read(const_cast<char *>(value.data()), value.length());
    }
};


// std::array
template <typename T, size_t N>
class Serializer<std::array<T, N>, std::enable_if_t<is_serializable<T>>> {
public:
    static void Write(std::ostream &os, const std::array<T, N> &v) {
        for (size_t i = 0; i < v.size(); ++i) {
            BinWrite(os, v[i]);
        }
    }

    static void Read(std::istream &is, std::array<T, N> &v) {
        for (size_t i = 0; i < v.size(); ++i) {
            BinRead(is, v[i]);
        }
    }
};

// std::vector
template <typename T>
class Serializer<std::vector<T>, std::enable_if_t<!std::is_same<T, bool>::value && is_serializable<T>>> {
public:
    static void Write(std::ostream &os, const std::vector<T> &v) {
        BinWrite(os, v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            BinWrite(os, v[i]);
        }
    }

    static void Read(std::istream &is, std::vector<T> &v) {
        size_t size;
        BinRead(is, size);
        v.resize(size);
        for (size_t i = 0; i < size; ++i) {
            BinRead(is, v[i]);
        }
    }
};

// std::tuple
namespace detail {
template <class F, class Tuple, std::size_t... I>
constexpr decltype(auto) apply_impl(F &&f, Tuple &&t, std::index_sequence<I...>) {
    return std::forward<F>(f)(std::get<I>(std::forward<Tuple>(t))...);
}

template <class T>
constexpr std::size_t tuple_size_v = std::tuple_size<T>::value;
}  // namespace detail

// Just std::apply from C++17
template <class F, class Tuple>
constexpr decltype(auto) apply(F &&f, Tuple &&t) {
    return detail::apply_impl(std::forward<F>(f), std::forward<Tuple>(t),
                              std::make_index_sequence<detail::tuple_size_v<std::remove_reference_t<Tuple>>>{});
}

template <typename... Ts>
class Serializer<std::tuple<Ts...>, std::enable_if_t<is_serializable<Ts...>>> {
public:
    static void Write(std::ostream &os, const std::tuple<Ts...> &t) {
        auto binwriter = [&os](const auto& v) { BinWrite(os, v); };
        apply(binwriter, t);
    }

    static void BinRead(std::istream &is, const std::tuple<Ts...> &t) {
        auto binreader = [&is](auto& v) { BinRead(is, v); };
        apply(binreader, t);
    }
};

// std::pair
template <typename T1, typename T2>
class Serializer<std::pair<T1, T2>, std::enable_if_t<is_serializable<T1, T2>>> {
public:
    static void Write(std::ostream &os, const std::pair<T1, T2> &p) {
        BinWrite(os, p.first, p.second);
    }

    static void Read(std::istream &is, std::pair<T1, T2> &p) {
        BinRead(is, p.first, p.second);
    }
};

// std::map & std::unordered_map
template <typename M>
class MapSerializer {
public:
    static void Write(std::ostream &os, const M &m) {
        size_t size = m.size();
        BinWrite(os, size);
        for (const auto &kv : m) {
            BinWrite(os, kv.first, kv.second);
        }
    }

    static void Read(std::istream &is, M &m) {
        m.clear();
        size_t size;
        BinRead(is, size);
        for (size_t i = 0; i < size; ++i) {
            typename M::key_type k;
            typename M::mapped_type v;
            BinRead(is, k, v);
            m.insert({std::move(k), std::move(v)});
        }
    }
};

// std::unique_ptr
template <typename T>
class Serializer<std::unique_ptr<T>, std::enable_if_t<is_serializable<T>>> {
public:
    static void Write(std::ostream &os, const std::unique_ptr<T> &v) {
        if (!v) {
            BinWrite(os, false);
            return;
        }
        BinWrite(os, true);
        BinWrite(os, *v);
    }

    static void Read(std::istream &is, std::unique_ptr<T> &v) {
        bool hasValue;
        BinRead(is, hasValue);
        if (!hasValue) {
            v = nullptr;
            return;
        }

        if (!v)
            v.reset(new T());
        BinRead(is, *v);
    }
};

template <typename K, typename V, typename... Args>
class Serializer<std::map<K, V, Args...>, std::enable_if_t<is_serializable<K, V>>> : public MapSerializer<std::map<K, V, Args...>> {};

template <typename K, typename V, typename... Args>
class Serializer<std::unordered_map<K, V, Args...>, std::enable_if_t<is_serializable<K, V>>> : public MapSerializer<std::unordered_map<K, V, Args...>> {};

//---- Helpers ---------------------------------------------------------------------------------------------------------

/**
 * @brief A value-returning variant of BinRead to save some LOCs.
 */
template<typename T>
T BinRead(std::istream &is) {
    T result;
    BinRead(is, result);
    return result;
}

/**
 * @brief  A small wrapper for binary serialization into an STL ostream.
 */
class BinOStream {
public:
    BinOStream(std::ostream &str)
            : str_(str) {
    }

    template<typename T>
    BinOStream &operator<<(const T &value) {
        BinWrite(str_, value);
        return *this;
    }

    operator bool() const {
        return (bool) str_;
    }

    bool operator!() const {
        return !str_;
    }

private:
    std::ostream &str_;
};

/**
 * @brief  A small wrapper for binary deserialization from an STL istream.
 */
class BinIStream {
public:
    BinIStream(std::istream &str)
            : str_(str) {
    }

    template<typename T>
    BinIStream &operator>>(T &value) {
        BinRead(str_, value);
        //VERIFY(!str.fail());
        return *this;
    }

    template<typename T>
    T Read() {
        T result;
        (*this) >> result;
        return result;
    }

    operator bool() const {
        return (bool) str_;
    }

    bool operator!() const {
        return !str_;
    }

    std::istream &stream() {
        return str_;
    }

private:
    std::istream &str_;
};

}  // namespace impl

using impl::BinWrite;
using impl::BinRead;
using impl::BinIStream;
using impl::BinOStream;
using impl::Serializer;
using impl::is_serializable;
using impl::has_binread_binwrite;
using impl::has_binarchive;
using impl::has_binarchive_save_load;
using impl::has_inside_serialization;
using impl::MapSerializer;

}  // namespace binary

}  // namespace io
