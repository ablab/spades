//
// Created by andrey on 14.11.16.
//

#pragma once

#include "modules/path_extend/scaffolder2015/extension_chooser2015.hpp"
#include "pipeline/config_struct.hpp"
#include "modules/path_extend/pe_config_struct.hpp"
#include "modules/path_extend/path_extender.hpp"
#include "modules/path_extend/scaffolder2015/connection_condition2015.hpp"
#include "launch_support.hpp"

namespace path_extend {

using namespace debruijn_graph;

//is here right place?
//FIXME no right place for this stuff; let's try to design better :)
//ANSWER: to utils or smth like it? actually it was only for debug, can remove it now, but the function may be useful
template<typename Base, typename T>
inline bool instanceof(const T *ptr) {
    return dynamic_cast<const Base *>(ptr) != nullptr;
}

class ExtendersGenerator {
    const config::dataset &dataset_info_;
    const PathExtendParamsContainer &params_;
    const conj_graph_pack &gp_;

    const GraphCoverageMap &cover_map_;

    PELaunchSupport support_;

public:
    ExtendersGenerator(const config::dataset &dataset_info,
                       const PathExtendParamsContainer &params,
                       const conj_graph_pack &gp,
                       const GraphCoverageMap &cover_map) :
        dataset_info_(dataset_info),
        params_(params),
        gp_(gp),
        cover_map_(cover_map),
        support_(dataset_info, params) { }

    typedef vector<shared_ptr<PathExtender>> Extenders;

    Extenders MakeMPExtenders(const ScaffoldingUniqueEdgeStorage &storage) const;


    Extenders MakePBScaffoldingExtenders(ScaffoldingUniqueEdgeStorage &unique_storage_pb,
                                         vector<PathContainer> &long_reads_paths,
                                         vector<shared_ptr<GraphCoverageMap>> &long_reads_cov_map) const;

    Extenders MakeBasicExtenders(const ScaffoldingUniqueEdgeStorage &storage) const;

private:


    void AddPathsToContainer(const std::vector<PathInfo<Graph>> &paths,
                             size_t size_threshold,
                             PathContainer &result) const;

    shared_ptr<ExtensionChooser> MakeLongReadsExtensionChooser(const config::dataset::Library &lib,
                                                               size_t lib_index) const;

    shared_ptr<SimpleExtender> MakeLongReadsExtender(size_t lib_index) const;

    shared_ptr<SimpleExtender> MakeLongEdgePEExtender(size_t lib_index,
                                                      bool investigate_loops) const;

    shared_ptr<SimpleExtensionChooser> MakeMetaExtensionChooser(shared_ptr<PairedInfoLibrary> lib,
                                                                size_t read_length) const;

    shared_ptr<SimpleExtender> MakeMetaExtender(size_t lib_index, bool investigate_loops) const;


    shared_ptr<SimpleExtender> MakePEExtender(size_t lib_index, bool investigate_loops) const;


    shared_ptr<GapJoiner> MakeGapJoiners(shared_ptr<PairedInfoLibrary> paired_lib) const;


    shared_ptr<PathExtender> MakeScaffoldingExtender(size_t lib_index) const;


    shared_ptr<PathExtender> MakeRNAScaffoldingExtender(size_t lib_index) const;


    shared_ptr<PathExtender> MakeScaffolding2015Extender(size_t lib_index, const ScaffoldingUniqueEdgeStorage &storage) const;


    shared_ptr<SimpleExtender> MakeCoordCoverageExtender(size_t lib_index) const;


    shared_ptr<SimpleExtender> MakeRNAExtender(size_t lib_index, bool investigate_loops) const;


    void PrintExtenders(vector<shared_ptr<PathExtender>> &extenders) const;

};

}
