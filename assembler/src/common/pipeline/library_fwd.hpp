//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

// Forward decls for YAML API
namespace llvm { namespace yaml { class IO; template<typename T> struct MappingTraits; } }
namespace llvm { class StringRef;  }

namespace io {

enum class LibraryType {
    SingleReads,
    SangerReads,
    RNA10x,
    PacBioReads,
    NanoporeReads,
    PairedEnd,
    HQMatePairs,
    MatePairs,
    TrustedContigs,
    TSLReads,
    PathExtendContigs,
    UntrustedContigs,
    FLRNAReads
};

enum class LibraryOrientation {
  FR,
  FF,
  RF,
  RR,
  Undefined
};

class SequencingLibraryBase;

struct NoData {};

template<class Data = NoData>
class SequencingLibrary;

template<class Data = NoData>
class DataSet;

}

namespace llvm { namespace yaml {
template <>
struct MappingTraits<io::SequencingLibraryBase> {
    static void mapping(llvm::yaml::IO &io, io::SequencingLibraryBase &lib);
    static StringRef validate(llvm::yaml::IO &io, io::SequencingLibraryBase &lib);
};

template <class Data>
struct MappingTraits<io::SequencingLibrary<Data> > {
    static void mapping(llvm::yaml::IO &io, io::SequencingLibrary<Data> &lib);
    static StringRef validate(llvm::yaml::IO &io, io::SequencingLibrary<Data> &lib);
};

template <class Data>
struct MappingTraits<io::DataSet<Data> > {
    static void mapping(llvm::yaml::IO &io, io::DataSet<Data> &lib);
};

}}
