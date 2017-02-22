//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "library.hpp"

#include "library/library_fwd.hpp"
#include "llvm/Support/Errc.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/YAMLTraits.h"

#include <filesystem>
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
        io.enumCase(value, "tell-seq",            LibraryType::TellSeqReads);
        io.enumCase(value, "trusted-contigs",     LibraryType::TrustedContigs);
        io.enumCase(value, "untrusted-contigs",   LibraryType::UntrustedContigs);
        io.enumCase(value, "path-extend-contigs", LibraryType::PathExtendContigs);
        io.enumCase(value, "fl-rna",              LibraryType::FLRNAReads);
        io.enumCase(value, "assembly-graph",      LibraryType::AssemblyGraph);
        io.enumCase(value, "clouds10x",           LibraryType::Clouds10x);
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
    io.mapOptional("number", number_, -1);
    io.mapOptional("orientation", orientation_, LibraryOrientation::Undefined);
    io.mapOptional("left reads", left_paired_reads_);
    io.mapOptional("right reads", right_paired_reads_);
    io.mapOptional("interlaced reads", interlaced_reads_);
    io.mapOptional("merged reads", merged_reads_);
    io.mapOptional("single reads", single_reads_);
    io.mapOptional("aux", aux_reads_);
}

void SequencingLibraryBase::validate(llvm::yaml::IO &, llvm::StringRef &res) {
    switch (type_) {
    case LibraryType::PairedEnd:
    case LibraryType::MatePairs:
    case LibraryType::HQMatePairs:
    case LibraryType::TellSeqReads:
    case LibraryType::Clouds10x:
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
    case LibraryType::FLRNAReads:
        if (left_paired_reads_.size() || right_paired_reads_.size() || interlaced_reads_.size()) {
            res = "Paired reads should not be set for this library type";
            return;
        }
      break;
    case LibraryType::AssemblyGraph:
        if (left_paired_reads_.size() || right_paired_reads_.size() || interlaced_reads_.size()) {
            res = "Paired reads should not be set for this library type";
            return;
        }
        if (merged_reads_.size() || single_reads_.size() > 1) {
            res = "Only single assembly graph must be specified";
            return;
        }
      break;
    default:
        // Impossible
        res = "Unsupported library type";
        return;
  }

  if (type_  == LibraryType::TellSeqReads && aux_reads_.size() != left_paired_reads_.size()) {
    res = "Aux stream (with barcodes) must be specified for TellSeq reads";
    return;
  }
}

void SequencingLibraryBase::update_relative_reads_filenames(const std::filesystem::path &input_dir) {
    auto update_relative_filename = [&input_dir](const std::filesystem::path &filename) {
        return input_dir / filename;
    };

    std::transform(left_paired_reads_.begin(), left_paired_reads_.end(), left_paired_reads_.begin(),
                   update_relative_filename);
    std::transform(right_paired_reads_.begin(), right_paired_reads_.end(), right_paired_reads_.begin(),
                   update_relative_filename);
    std::transform(interlaced_reads_.begin(), interlaced_reads_.end(), interlaced_reads_.begin(),
                   update_relative_filename);
    std::transform(merged_reads_.begin(), merged_reads_.end(), merged_reads_.begin(),
                   update_relative_filename);
    std::transform(aux_reads_.begin(), aux_reads_.end(), aux_reads_.begin(),
                   update_relative_filename);
    std::transform(single_reads_.begin(), single_reads_.end(), single_reads_.begin(),
                   update_relative_filename);
}

#include "library.inl"

// Provide default implementation here (e.g. in case of Data == io::NoData)
template class io::DataSet<>;
