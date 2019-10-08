//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/io/binary/binary.hpp"
#include <parallel_hashmap/phmap.h>

template <typename K, typename V, typename... Args>
class io::binary::Serializer<phmap::flat_hash_map<K, V, Args...>, std::enable_if_t<io::binary::is_serializable<K, V>>> : public io::binary::MapSerializer<phmap::flat_hash_map<K, V, Args...>> {};

template <typename K, typename V, typename... Args>
class io::binary::Serializer<phmap::node_hash_map<K, V, Args...>, std::enable_if_t<io::binary::is_serializable<K, V>>> : public io::binary::MapSerializer<phmap::node_hash_map<K, V, Args...>> {};
