/// TODO: Remove this awful file when C++17 comes!

#pragma once

#include <type_traits>
#include <cstddef>

namespace traits {

namespace details {

template<class T>
constexpr static bool ContainsHandler() {
    return false;
}

template<class T, class Head, class ... Tail>
constexpr static bool ContainsHandler() {
    return (std::is_same<T, Head>::value ? true : ContainsHandler<T, Tail ...>());
}


template<class T, T el>
constexpr bool ContainsHandler() {
    return false;
}

template<class T, T el, T head, T ... tail>
constexpr bool ContainsHandler() {
    return (el == head ? true : ContainsHandler<T, el, tail ...>());
}


template<size_t pos, class T, T el>
constexpr static size_t GetIndexHandler() {
    return -1ull;
}

template<size_t pos, class T, T el, T head, T ... tail>
constexpr static size_t GetIndexHandler() {
    return (el == head ? pos : GetIndexHandler<pos+1, T, el, tail ...>());
}

} // namespace details

template<class T, class ... Types>
constexpr bool Contains() {
    return details::ContainsHandler<T, Types ...>();
}

template<class T, T el, T ... elements>
constexpr bool Contains() {
    return details::ContainsHandler<T, el, elements ...>();
}

template<class T, T el, T ... elements>
constexpr static size_t GetIndex() {
    return details::GetIndexHandler<0, T, el, elements ...>();
}


} // namespace traits
