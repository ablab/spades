//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "pipeline/genomic_info.hpp"

namespace io {

namespace binary {

class GenomicInfoIO : public IOBase<GenomicInfo> {
public:
    GenomicInfoIO()
            : IOBase<GenomicInfo>() {
    }

    void Save(const std::string &basename, const GenomicInfo &value) override {
        value.Save(basename + ext_);
    }

    bool Load(const std::string &basename, GenomicInfo &value) override {
        return value.Load(basename + ext_);
    }

    void Write(std::ostream &os, const GenomicInfo &value) {
        value.BinWrite(os);
    }

    bool Read(std::istream &is, GenomicInfo &value) {
        value.BinRead(is);
        return true;
    }

protected:
    static constexpr const char *ext_ = ".ginfo";
};

template<>
struct IOTraits<GenomicInfo> {
    typedef GenomicInfoIO Type;
};

} // namespace binary

} // namespace io
