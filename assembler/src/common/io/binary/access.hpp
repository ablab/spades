//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace io {
namespace binary {
class access {
public:
    template <typename T, typename Archive>
    auto BinArchive(T &v, Archive &ar) -> decltype(v.BinArchive(ar), void()) {
        v.BinArchive(ar);
    }

    template <typename T, typename Archive>
    auto BinArchiveLoad(T &v, Archive &ar) -> decltype(v.BinArchiveLoad(ar), void()) {
        v.BinArchiveLoad(ar);
    }

    template <typename T, typename Archive>
    auto BinArchiveSave(T &v, Archive &ar) -> decltype(v.BinArchiveSave(ar), void()) {
        v.BinArchiveSave(ar);
    }
};

}  // namespace binary
}  // namespace io
