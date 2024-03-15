//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/io/binary/binary.hpp"
#include <filesystem>

namespace io {
namespace binary {
namespace impl {

template<>
class Serializer<std::filesystem::path> {
public:
    static void Write(std::ostream &os, const std::filesystem::path &p) {
        io::binary::BinWrite(os, p.native());
    }

    static void Read(std::istream &is, std::filesystem::path &p) {
        p = io::binary::BinRead<std::string>(is);
    }
};
}  // namespace impl
}  // namespace binary
}  // namespace io
