//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "pipeline/library.hpp"
#include "utils/filesystem/path_helper.hpp"

#include "llvm/Support/YAMLTraits.h"
#include "llvm/Support/Errc.h"
#include "llvm/Support/FileSystem.h"

#include <string>
#include <fstream>
#include <iostream>

using namespace llvm;
using namespace io;

namespace llvm { namespace yaml {
template <>
struct ScalarEnumerationTraits<LibraryOrientation> {
    static void enumeration(yaml::IO &io, LibraryOrientation &value) {
        io.enumCase(value, "fr", LibraryOrientation::FR);
        io.enumCase(value, "rf", LibraryOrientation::RF);
        io.enumCase(value, "ff", LibraryOrientation::FF);
        io.enumCase(value, "rr", LibraryOrientation::RR);
    }
};

template <>
struct ScalarEnumerationTraits<LibraryType> {
    static void enumeration(yaml::IO &io, LibraryType &value) {
        io.enumCase(value, "paired-end",          LibraryType::PairedEnd);
        io.enumCase(value, "mate-pairs",          LibraryType::MatePairs);
        io.enumCase(value, "hq-mate-pairs",       LibraryType::HQMatePairs);
        io.enumCase(value, "pacbio",              LibraryType::PacBioReads);
        io.enumCase(value, "single",              LibraryType::SingleReads);
        io.enumCase(value, "sanger",              LibraryType::SangerReads);
        io.enumCase(value, "nanopore",            LibraryType::NanoporeReads);
        io.enumCase(value, "tslr",                LibraryType::TSLReads);
        io.enumCase(value, "trusted-contigs",     LibraryType::TrustedContigs);
        io.enumCase(value, "untrusted-contigs",   LibraryType::UntrustedContigs);
        io.enumCase(value, "path-extend-contigs", LibraryType::PathExtendContigs);
    }
};

template <>
struct SequenceTraits<std::vector<std::string>> {
    static size_t size(IO &, std::vector<std::string> &seq) {
        return seq.size();
    }
    static std::string&
    element(IO &, std::vector<std::string> &seq, size_t index) {
        if (index >= seq.size())
            seq.resize(index+1);
        return seq[index];
    }
};
}}

namespace io {
template<>
void SequencingLibrary<io::NoData>::yamlize(llvm::yaml::IO &io) {
  SequencingLibraryBase::yamlize(io);
}
template<>
void SequencingLibrary<io::NoData>::validate(llvm::yaml::IO &io, llvm::StringRef &res) {
  SequencingLibraryBase::validate(io, res);
}
}

void SequencingLibraryBase::yamlize(llvm::yaml::IO &io) {
    io.mapRequired("type", type_);
    io.mapOptional("orientation", orientation_, LibraryOrientation::Undefined);
    io.mapOptional("left reads", left_paired_reads_);
    io.mapOptional("right reads", right_paired_reads_);
    io.mapOptional("merged reads", merged_reads_);
    io.mapOptional("single reads", single_reads_);
}

void SequencingLibraryBase::validate(llvm::yaml::IO &, llvm::StringRef &res) {
    switch (type_) {
    case LibraryType::PairedEnd:
    case LibraryType::MatePairs:
    case LibraryType::HQMatePairs:
        if (left_paired_reads_.size() != right_paired_reads_.size()) {
            res = "Left and right reads lists should have equal length";
            return;
        }

        if (orientation_ == LibraryOrientation::Undefined) {
            res = "Orientation for paired reads should be specified";
            return;
        }
        break;
    case LibraryType::SingleReads:
    case LibraryType::PacBioReads:
    case LibraryType::SangerReads:
    case LibraryType::NanoporeReads:
    case LibraryType::TSLReads:
    case LibraryType::TrustedContigs:
    case LibraryType::UntrustedContigs:
    case LibraryType::PathExtendContigs:
        if (left_paired_reads_.size() || right_paired_reads_.size()) {
            res = "Paired reads should not be set for this library type";
            return;
        }
      break;
    default:
        // Impossible
        res = "Unsupported library type";
        return;
  }
}

// FIXME: Lambda
struct update_relative_filename : public std::binary_function<std::string, std::string, std::string> {
    std::string operator() (const std::string &filename, const std::string &input_dir) const {
        if (filename[0] == '/')
            return filename;
        return input_dir + filename;
    }
};

void SequencingLibraryBase::update_relative_reads_filenames(const std::string &input_dir) {
    std::transform(left_paired_reads_.begin(), left_paired_reads_.end(), left_paired_reads_.begin(),
                   std::bind2nd(update_relative_filename(), input_dir));
    std::transform(right_paired_reads_.begin(), right_paired_reads_.end(), right_paired_reads_.begin(),
                   std::bind2nd(update_relative_filename(), input_dir));
    std::transform(merged_reads_.begin(), merged_reads_.end(), merged_reads_.begin(),
                   std::bind2nd(update_relative_filename(), input_dir));
    std::transform(single_reads_.begin(), single_reads_.end(), single_reads_.begin(),
                   std::bind2nd(update_relative_filename(), input_dir));
}

#include "pipeline/library.inl"

// Provide default implementation here (e.g. in case of Data == io::NoData)
template class io::DataSet<>;
