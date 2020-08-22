//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "genomic_info.hpp"

#include "utils/filesystem/file_opener.hpp"

#include "llvm/Support/YAMLTraits.h"
#include "llvm/Support/Errc.h"
#include "llvm/Support/FileSystem.h"

#include <fstream>

using namespace llvm;

namespace llvm { namespace yaml {
template <>
struct MappingTraits<GenomicInfo> {
    static void mapping(yaml::IO &io, GenomicInfo &info) {
      info.yamlize(io);
    }
};


template <>
struct SequenceTraits<std::vector<std::size_t>> {
    static size_t size(IO &, std::vector<std::size_t> &seq) {
        return seq.size();
    }
    static size_t&
    element(IO &, std::vector<std::size_t> &seq, size_t index) {
        if (index >= seq.size())
            seq.resize(index+1);
        return seq[index];
    }
    static const bool flow = true;
};
}}

void GenomicInfo::yamlize(yaml::IO &io) {
    io.mapOptional("ec bound", ec_bound_, 0.0);
    io.mapOptional("estimated mean", estimated_mean_, 0.0);
    io.mapOptional("trusted bound", trusted_bound_, 0.0);
    io.mapOptional("genome size", genome_size_, size_t(0));
    io.mapOptional("coverage histogram", cov_histogram_);
}

bool GenomicInfo::Load(const std::string &filename) {
    auto ifs = fs::open_file(filename, std::ios::binary);
    BinRead(ifs);
    return true;
}

void GenomicInfo::Save(const std::string &filename) const {
    std::ofstream ofs(filename, std::ios::binary);
    BinWrite(ofs);
}

bool GenomicInfo::BinWrite(std::ostream &os) const {
    std::string buffer;
    llvm::raw_string_ostream raw_os(buffer);
    llvm::yaml::Output yout(raw_os);
    yout << const_cast<GenomicInfo&>(*this);
    raw_os.str();  // Flush content to the target string
    io::binary::BinWrite(os, buffer);
    return true;
}

void GenomicInfo::BinRead(std::istream &is) {
    std::string buffer;
    io::binary::BinRead(is, buffer);
    llvm::yaml::Input yin(buffer);
    yin >> *this;
    VERIFY(!yin.error());
}
