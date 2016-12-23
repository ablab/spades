//
// Created by andrey on 10.10.16.
//

#pragma once

//#include "path_extend_launch.hpp"

#include <common/assembly_graph/paths/bidirectional_path.hpp>
#include "modules/path_extend/paired_library.hpp"
#include "pipeline/config_struct.hpp"
#include "modules/path_extend/pe_config_struct.hpp"

namespace path_extend {

using namespace debruijn_graph;

inline size_t FindMaxISRightQuantile(const config::dataset& dataset_info, bool include_mate_pairs = true) {
    size_t res = 0;
    for (const auto& lib : dataset_info.reads) {
        if (lib.is_paired()) {
            if (lib.is_mate_pair() && !include_mate_pairs)
                continue;
            res = max(res, (size_t) lib.data().insert_size_right_quantile);
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
                              const std::string& output_dir_,
                              const std::string& contigs_name_,
                              const std::string& scf_name_,
                              config::pipeline_type mode_,
                              bool uneven_depth_,
                              bool avoid_rc_connections_,
                              bool use_scaffolder_,
                              bool output_broken_scaffolds_ = true):
        pe_cfg(pe_cfg_),
        pset(pe_cfg_.param_set),
        output_dir(output_dir_),
        etc_dir(output_dir + pe_cfg_.etc_dir + "/"),
        contigs_name(scf_name_),
        broken_contigs(contigs_name_),
        mode(mode_),
        uneven_depth(uneven_depth_),
        avoid_rc_connections(avoid_rc_connections_),
        use_scaffolder(use_scaffolder_),
        traverse_loops(true),
        output_broken_scaffolds(output_broken_scaffolds_),
        detect_repeats_online(mode_ != config::pipeline_type::meta && mode_ != config::pipeline_type::rna)
    {
        if (!(use_scaffolder && pset.scaffolder_options.enabled)) {
            contigs_name = contigs_name_;
            traverse_loops = false;
            output_broken_scaffolds = false;
        }
        if (mode_ == config::pipeline_type::rna)
            traverse_loops = false;

        //Parameters are subject to change
        max_polisher_gap = FindMaxISRightQuantile(dataset_info);

        //TODO: params
        if (HasLongReads(dataset_info))
            max_polisher_gap = max(max_polisher_gap, size_t(10000));

        min_edge_len = 100;
        max_path_diff = mode == config::pipeline_type::rna ? 1 : max_polisher_gap;
    }

    const pe_config::MainPEParamsT& pe_cfg;
    const pe_config::ParamSetT& pset;

    std::string output_dir;
    std::string etc_dir;

    std::string contigs_name;
    std::string broken_contigs;

    config::pipeline_type mode;
    bool uneven_depth;

    bool avoid_rc_connections;
    bool use_scaffolder;
    bool traverse_loops;
    bool output_broken_scaffolds;
    bool detect_repeats_online;

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

    pe_config::ParamSetT::ExtensionOptionsT GetExtensionOpts(shared_ptr<PairedInfoLibrary> lib, const pe_config::ParamSetT& pset) const;

    void SetSingleThresholdForLib(shared_ptr<PairedInfoLibrary> lib, const pe_config::ParamSetT &pset, double threshold) const;

    bool HasOnlyMPLibs() const;

    bool IsForSingleReadExtender(const io::SequencingLibrary<config::DataSetData> &lib) const;

    bool IsForPEExtender(const io::SequencingLibrary<config::DataSetData> &lib) const;

    bool IsForShortLoopExtender(const io::SequencingLibrary<config::DataSetData> &lib) const;

    bool IsForScaffoldingExtender(const io::SequencingLibrary<config::DataSetData> &lib) const;

    bool UseCoverageResolverForSingleReads(const io::LibraryType& type) const;

    std::string LibStr(size_t count) const;

    pe_config::LongReads GetLongReadsConfig(const io::LibraryType &type) const;

    size_t FindMaxMPIS() const;

    bool HasLongReads() const;

    bool HasMPReads() const;

    bool SingleReadsMapped() const;

    double EstimateLibCoverage(size_t lib_index) const;

    size_t TotalNuclsInGraph() const;

    bool NeedsUniqueEdgeStorage() const;

};

}
