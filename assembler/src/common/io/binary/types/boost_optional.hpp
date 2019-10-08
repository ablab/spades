//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/io/binary/binary.hpp"
#include <boost/optional.hpp>

template <typename T>
class io::binary::Serializer<boost::optional<T>, std::enable_if_t<io::binary::is_serializable<T>>> {
public:
    static void Write(std::ostream &os, const boost::optional<T> &v) {
        if (v) {
            io::binary::BinWrite(os, true, *v);
        } else {
            io::binary::BinWrite(os, false);
        }
    }

    static void Read(std::istream &is, boost::optional<T> &v) {
        auto present = io::binary::BinRead<bool>(is);
        if (present) {
            auto val = io::binary::BinRead<T>(is);
            v = boost::make_optional(val);
        } else {
            v = boost::none;
        }
    }
};
