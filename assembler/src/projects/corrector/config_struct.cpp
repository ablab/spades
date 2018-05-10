//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "config_struct.hpp"

#include "utils/parallel/openmp_wrapper.h"

#include "llvm/Support/YAMLParser.h"
#include "llvm/Support/YAMLTraits.h"

#include <string>

using namespace llvm;

namespace llvm { namespace yaml {
template <>
struct ScalarEnumerationTraits<corrector::Strategy> {
    static void enumeration(yaml::IO &io, corrector::Strategy &value) {
        io.enumCase(value, "all_reads",      corrector::Strategy::AllReads);
        io.enumCase(value, "majority_only",  corrector::Strategy::MajorityOnly);
        io.enumCase(value, "not_started",    corrector::Strategy::AllExceptJustStarted);
        io.enumCase(value, "mapped_squared", corrector::Strategy::MappedSquared);
    }
};
}}

// FIXME: This is temporary
class DataSetReader {
  public:
    DataSetReader(yaml::IO&) {}
    DataSetReader(yaml::IO&, io::DataSet<>&) {}

    io::DataSet<> denormalize(yaml::IO &) {
        return io::DataSet<>(path);
    }

    std::string path;
};

namespace llvm { namespace yaml {
template <>
struct MappingTraits<corrector::corrector_config> {
    static void mapping(yaml::IO &io, corrector::corrector_config &cfg) {
        yaml::MappingNormalization<DataSetReader, io::DataSet<>> dataset(io, cfg.dataset);

        io.mapRequired("dataset", dataset->path);
        io.mapOptional("work_dir", cfg.work_dir, std::string("."));
        io.mapOptional("output_dir", cfg.output_dir, std::string("."));
        io.mapOptional("max_nthreads", cfg.max_nthreads, 1u);
        io.mapRequired("strategy", cfg.strat);
        io.mapOptional("bwa", cfg.bwa, std::string("."));
        io.mapOptional("log_filename", cfg.log_filename, std::string("."));
    }
};
}}


namespace corrector {
void load(corrector_config& cfg, const std::string &filename) {
    ErrorOr<std::unique_ptr<MemoryBuffer>> Buf = MemoryBuffer::getFile(filename);
    if (!Buf)
        throw(std::string("Failed to load config file ") + filename);

    yaml::Input yin(*Buf.get());
    yin >> cfg;

    if (yin.error())
        throw(std::string("Failed to load config file ") + filename);

    cfg.max_nthreads = spades_set_omp_threads(cfg.max_nthreads);
}
}
