//
// Created by andrey on 10.10.16.
//

#pragma once


#include "modules/path_extend/paired_library.hpp"
#include "configs/config_struct.hpp"
#include "configs/pe_config_struct.hpp"

namespace path_extend {

using namespace debruijn_graph;

inline size_t FindMaxISRightQuantile(const config::dataset& dataset_info, bool include_mate_pairs = true) {
    size_t res = 0;
    for (const auto& lib : dataset_info.reads) {
        if (lib.is_paired()) {
            if (lib.is_mate_pair() && !include_mate_pairs)
                continue;
            res = std::max(res, (size_t) lib.data().insert_size_right_quantile);
        }
    }
    return res;
}

inline bool HasLongReads(const config::dataset& dataset_info) {
    for (const auto& lib : dataset_info.reads) {
        if (lib.is_long_read_lib() || lib.is_contig_lib()) {
            return true;
        }
    }
    return false;
}

struct PathExtendParamsContainer {

    PathExtendParamsContainer(const config::dataset& dataset_info,
                              const pe_config::MainPEParamsT& pe_cfg_,
                              const config::debruijn_config::strand_specificity& strand_specificity,
                              const std::filesystem::path& output_dir_,
                              config::pipeline_type mode_,
                              bool uneven_depth_,
                              bool avoid_rc_connections_,
                              bool use_scaffolder_):
        pe_cfg(pe_cfg_),
        pset(pe_cfg_.param_set),
        ss(strand_specificity),
        output_dir(output_dir_),
        etc_dir(output_dir / pe_cfg_.etc_dir),
        mode(mode_),
        uneven_depth(uneven_depth_),
        avoid_rc_connections(avoid_rc_connections_),
        use_scaffolder(use_scaffolder_),
        traverse_loops(true)
    {
        if (!(use_scaffolder && pset.scaffolder_options.enabled)) {
            traverse_loops = false;
        }
        if (mode_ == config::pipeline_type::rna)
            traverse_loops = false;

        //Parameters are subject to change
        max_polisher_gap = FindMaxISRightQuantile(dataset_info);
        //TODO: params
        if (HasLongReads(dataset_info))
            max_polisher_gap = std::max(max_polisher_gap, size_t(10000));

        min_edge_len = 0;

        max_path_diff = FindMaxISRightQuantile(dataset_info);
        if (mode == config::pipeline_type::rna || config::PipelineHelper::IsMetagenomicPipeline(mode))
            max_path_diff = 0;
    }

    const pe_config::MainPEParamsT& pe_cfg;
    const pe_config::ParamSetT& pset;

    const config::debruijn_config::strand_specificity& ss;

    std::filesystem::path output_dir;
    std::filesystem::path etc_dir;

    config::pipeline_type mode;
    bool uneven_depth;

    bool avoid_rc_connections;
    bool use_scaffolder;
    bool traverse_loops;

    //todo move to config
    size_t min_edge_len;
    size_t max_path_diff;
    size_t max_polisher_gap;
    //TODO: move here size_t max_repeat_length;
};


class PELaunchSupport {
    const config::dataset& dataset_info_;
    const PathExtendParamsContainer& params_;

public:

    PELaunchSupport(const config::dataset& dataset_info,
                    const PathExtendParamsContainer& params):
        dataset_info_(dataset_info),
        params_(params) { }

    pe_config::ParamSetT::ExtensionOptionsT GetExtensionOpts(std::shared_ptr<PairedInfoLibrary> lib,
                                                             const pe_config::ParamSetT &pset) const;

    bool HasOnlyMPLibs() const;

    bool HasOnlySingleReads() const;

    bool IsForSingleReadExtender(const io::SequencingLibrary<config::LibraryData> &lib) const;

    bool IsForSingleReadScaffolder(const io::SequencingLibrary<config::LibraryData> &lib) const;

    bool IsForPEExtender(const io::SequencingLibrary<config::LibraryData> &lib) const;

    bool IsForShortLoopExtender(const io::SequencingLibrary<config::LibraryData> &lib) const;

    bool IsForScaffoldingExtender(const io::SequencingLibrary<config::LibraryData> &lib) const;

    bool UseCoverageResolverForSingleReads(const io::LibraryType& type) const;

    std::string LibStr(size_t count) const;

    pe_config::LongReads GetLongReadsConfig(const io::LibraryType &type) const;

    size_t FindMaxMPIS() const;

    bool HasLongReads() const;

    bool HasLongReadsScaffolding() const;

    bool HasMPReads() const;

    bool SingleReadsMapped() const;

    double EstimateLibCoverage(size_t lib_index) const;

    size_t TotalNuclsInGraph() const;

    bool NeedsUniqueEdgeStorage() const;

};

}
