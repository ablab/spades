//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "config_struct.hpp"

#include "utils/openmp_wrapper.h"

#include "llvm/Support/YAMLParser.h"
#include "llvm/Support/YAMLTraits.h"

#include <string>

using namespace llvm;

namespace llvm { namespace yaml {
template <>
struct ScalarEnumerationTraits<hammer_config::HammerStage> {
    static void enumeration(yaml::IO &io, hammer_config::HammerStage &value) {
        io.enumCase(value, "count",       hammer_config::HammerStage::KMerCounting);
        io.enumCase(value, "hamcluster",  hammer_config::HammerStage::HammingClustering);
        io.enumCase(value, "subcluster",  hammer_config::HammerStage::SubClustering);
        io.enumCase(value, "correct",     hammer_config::HammerStage::ReadCorrection);
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
struct MappingTraits<hammer_config::hammer_config> {
    static void mapping(yaml::IO &io, hammer_config::hammer_config &cfg) {
        yaml::MappingNormalization<DataSetReader, io::DataSet<>> dataset(io, cfg.dataset);

        io.mapRequired("dataset", dataset->path);
        io.mapOptional("working_dir", cfg.working_dir, std::string("."));
        io.mapOptional("output_dir", cfg.output_dir, std::string("."));
        io.mapRequired("hard_memory_limit", cfg.hard_memory_limit);
        io.mapOptional("count_split_buffer", cfg.count_split_buffer, 0ul);
        io.mapOptional("max_nthreads", cfg.max_nthreads, 1u);
        io.mapRequired("kmer_qual_threshold", cfg.kmer_qual_threshold);
        io.mapRequired("center_qual_threshold", cfg.center_qual_threshold);
        io.mapRequired("delta_score_threshold", cfg.delta_score_threshold);
        io.mapRequired("keep_uncorrected_ends", cfg.keep_uncorrected_ends);
        io.mapRequired("tau", cfg.tau);
        io.mapOptional("debug_mode", cfg.debug_mode, false);
        io.mapOptional("start_stage", cfg.start_stage, hammer_config::HammerStage::KMerCounting);
    }
};
}}

namespace hammer_config {
void load(hammer_config& cfg, const std::string &filename) {
    ErrorOr<std::unique_ptr<MemoryBuffer>> Buf = MemoryBuffer::getFile(filename);
    if (!Buf)
        throw(std::string("Failed to load config file ") + filename);

    yaml::Input yin(*Buf.get());
    yin >> cfg;

    if (yin.error())
        throw(std::string("Failed to load config file ") + filename);

    // Fix number of threads according to OMP capabilities.
    cfg.max_nthreads = std::min(cfg.max_nthreads, (unsigned)omp_get_max_threads());
    // Inform OpenMP runtime about this :)
    omp_set_num_threads(cfg.max_nthreads);
}
}
