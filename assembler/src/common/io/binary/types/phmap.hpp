//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <parallel_hashmap/phmap.h>
#include "common/io/binary/binary.hpp"

namespace io {
namespace binary {
namespace impl {

// TODO Use nested template here
template <typename K, typename V, typename... Args>
class Serializer<phmap::flat_hash_map<K, V, Args...>, std::enable_if_t<io::binary::is_serializable<K, V>>>
    : public io::binary::MapSerializer<phmap::flat_hash_map<K, V, Args...>> {};

template <typename K, typename V, typename... Args>
class Serializer<phmap::node_hash_map<K, V, Args...>, std::enable_if_t<io::binary::is_serializable<K, V>>>
    : public io::binary::MapSerializer<phmap::node_hash_map<K, V, Args...>> {};
}  // namespace impl
}  // namespace binary
}  // namespace io
