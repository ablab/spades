//
// Created by andrey on 10.10.16.
//

#include "launch_support.hpp"

namespace path_extend {

using namespace debruijn_graph;

bool PELaunchSupport::HasOnlyMPLibs() const {
    for (const auto& lib : dataset_info_.reads)  {
        if (!(lib.is_mate_pair() && lib.data().mean_insert_size > 0.0)) {
            return false;
        }
    }
    return true;
}

//maybe pass lib as const ref? But ok if always used as shared_ptr<> always
//ANSWER: it is always shared_ptr
pe_config::ParamSetT::ExtensionOptionsT PELaunchSupport::GetExtensionOpts(shared_ptr<PairedInfoLibrary> lib,
                                                                          const pe_config::ParamSetT& pset) const {
    return lib->IsMp() ? pset.mate_pair_options : pset.extension_options;
}

void PELaunchSupport::SetSingleThresholdForLib(shared_ptr<PairedInfoLibrary> lib, const pe_config::ParamSetT &pset, double threshold) const {
    double t = pset.extension_options.use_default_single_threshold || math::le(threshold, 0.0) ?
               pset.extension_options.single_threshold : threshold;
    lib->SetSingleThreshold(t);
}

bool PELaunchSupport::IsForSingleReadExtender(const io::SequencingLibrary<config::DataSetData> &lib) const {
    return (lib.data().single_reads_mapped || lib.is_long_read_lib() || lib.is_contig_lib());
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

//What is it by the way?
//ANSWER: used coverage loop resolver for nextera pipeline, will review it
bool PELaunchSupport::UseCoverageResolverForSingleReads(const io::LibraryType& type) const {
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
}
