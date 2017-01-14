//
// Created by andrey on 10.10.16.
//

#include "launch_support.hpp"

namespace path_extend {

using namespace debruijn_graph;

bool PELaunchSupport::HasOnlyMPLibs() const {
    for (const auto &lib : dataset_info_.reads) {
        if (!(lib.is_mate_pair() && lib.data().mean_insert_size > 0.0)) {
            return false;
        }
    }
    return true;
}

pe_config::ParamSetT::ExtensionOptionsT PELaunchSupport::GetExtensionOpts(shared_ptr<PairedInfoLibrary> lib,
                                                                          const pe_config::ParamSetT &pset) const {
    return lib->IsMp() ? pset.mate_pair_options : pset.extension_options;
}

double PELaunchSupport::SingleThresholdForLib(const pe_config::ParamSetT &pset,
                                               double threshold) const {
    return pset.extension_options.use_default_single_threshold || math::le(threshold, 0.0) ?
               pset.extension_options.single_threshold : threshold;
}

bool PELaunchSupport::IsForSingleReadExtender(const io::SequencingLibrary<config::DataSetData> &lib) const {
    return (lib.data().single_reads_mapped || lib.is_long_read_lib() || lib.is_contig_lib());
}
bool PELaunchSupport::IsForSingleReadScaffolder(const io::SequencingLibrary<config::DataSetData> &lib) const {
    return (lib.is_long_read_lib() || (lib.is_contig_lib() && lib.type() != io::LibraryType::PathExtendContigs));
}

bool PELaunchSupport::IsForPEExtender(const io::SequencingLibrary<config::DataSetData> &lib) const {
    return (lib.type() == io::LibraryType::PairedEnd && lib.data().mean_insert_size > 0.0);
}

bool PELaunchSupport::IsForShortLoopExtender(const io::SequencingLibrary<config::DataSetData> &lib) const {
    return (lib.type() == io::LibraryType::PairedEnd && lib.data().mean_insert_size > 0.0);
}

bool PELaunchSupport::IsForScaffoldingExtender(const io::SequencingLibrary<config::DataSetData> &lib) const {
    return (lib.type() == io::LibraryType::PairedEnd && lib.data().mean_insert_size > 0.0);
}

//TODO: review usage
bool PELaunchSupport::UseCoverageResolverForSingleReads(const io::LibraryType &type) const {
    return HasOnlyMPLibs() && (type == io::LibraryType::HQMatePairs);
}

std::string PELaunchSupport::LibStr(size_t count) const {
    return count == 1 ? "library" : "libraries";
}

pe_config::LongReads PELaunchSupport::GetLongReadsConfig(const io::LibraryType &type) const {
    if (io::SequencingLibraryBase::is_long_read_lib(type)) {
        return params_.pe_cfg.long_reads.pacbio_reads;
    } else if (type == io::LibraryType::PathExtendContigs) {
        return params_.pe_cfg.long_reads.meta_contigs;
    } else if (io::SequencingLibraryBase::is_contig_lib(type)) {
        return params_.pe_cfg.long_reads.contigs;
    }
    return params_.pe_cfg.long_reads.single_reads;
}

size_t PELaunchSupport::FindMaxMPIS() const {
    size_t max_is = 0;
    for (size_t i = 0; i < dataset_info_.reads.lib_count(); ++i) {
        if (dataset_info_.reads[i].is_mate_pair()) {
            max_is = max(max_is, (size_t) dataset_info_.reads[i].data().mean_insert_size);
        }
    }
    return max_is;
}

bool PELaunchSupport::HasLongReads() const {
    return path_extend::HasLongReads(dataset_info_);
}

bool PELaunchSupport::HasLongReadsScaffolding() const {
    for (const auto &lib : dataset_info_.reads) {
        if (IsForSingleReadScaffolder(lib))
            return true;
    }
    return false;
}

bool PELaunchSupport::HasMPReads() const {
    for (const auto &lib : dataset_info_.reads) {
        if (lib.is_mate_pair()) {
            return true;
        }
    }
    return false;
}
bool PELaunchSupport::SingleReadsMapped() const {
    for (const auto &lib : dataset_info_.reads) {
        if (lib.data().single_reads_mapped) {
            return true;
        }
    }
    return false;
}

double PELaunchSupport::EstimateLibCoverage(size_t lib_index) const {
    double cov_fraction = double(dataset_info_.reads[lib_index].data().total_nucls) / double(TotalNuclsInGraph());
    return cov_fraction * dataset_info_.avg_coverage();
}

size_t PELaunchSupport::TotalNuclsInGraph() const {
    size_t total_nc_count = 0;
    for (const auto &lib: dataset_info_.reads) {
        if (lib.is_graph_contructable())
            total_nc_count += lib.data().total_nucls;
    }
    return total_nc_count;
}


bool PELaunchSupport::NeedsUniqueEdgeStorage() const {
    return !(params_.pset.sm == sm_old ||
        (params_.pset.sm == sm_old_pe_2015 && !HasLongReadsScaffolding() && !HasMPReads()));
}
}
