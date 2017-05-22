//
// Created by andrey on 14.11.16.
//

#ifndef PROJECT_LAUNCHER_H
#define PROJECT_LAUNCHER_H

#include "launch_support.hpp"
#include "extenders_logic.hpp"

#include "modules/path_extend/pe_resolver.hpp"
#include "modules/genome_consistance_checker.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

#include "modules/alignment/rna/ss_coverage.hpp"

namespace path_extend {

using namespace debruijn_graph;

class PathExtendLauncher {
    const config::dataset& dataset_info_;
    const PathExtendParamsContainer& params_;
    conj_graph_pack& gp_;
    PELaunchSupport support_;

    shared_ptr<ContigNameGenerator> contig_name_generator_;
    ContigWriter writer_;

    UniqueData unique_data_;

    vector<shared_ptr<ConnectionCondition>>
        ConstructPairedConnectionConditions(const ScaffoldingUniqueEdgeStorage& edge_storage) const;

    shared_ptr<scaffold_graph::ScaffoldGraph>
        ConstructScaffoldGraph(const ScaffoldingUniqueEdgeStorage& edge_storage) const;

    void PrintScaffoldGraph(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                            const set<EdgeId>& main_edge_set,
                            const debruijn_graph::GenomeConsistenceChecker& genome_checker,
                            const string& filename) const;

    void MakeAndOutputScaffoldGraph() const;

    void CountMisassembliesWithReference(const PathContainer& paths) const;

    void EstimateUniqueEdgesParams();

    void FillUniqueEdgeStorage();

    void FillPBUniqueEdgeStorages();

    void FillPathContainer(size_t lib_index, size_t size_threshold = 1);

    void FillLongReadsCoverageMaps();

    void DebugOutputPaths(const PathContainer& paths, const string& name) const;

    void RemoveOverlapsAndArtifacts(PathContainer &paths, GraphCoverageMap &cover_map, const PathExtendResolver &resolver) const;

    void CleanPaths(PathContainer &paths, const pe_config::ParamSetT::PathFiltrationT &path_filtration) const;

    size_t GetLengthCutoff(size_t abs_cutoff, double rel_cutoff) const;

    void TraverseLoops(PathContainer &paths, GraphCoverageMap &cover_map) const;

    void PolishPaths(const PathContainer &paths, PathContainer &result, const GraphCoverageMap &cover_map) const;

    Extenders ConstructExtenders(const GraphCoverageMap &cover_map, UsedUniqueStorage &used_unique_storage);

    Extenders ConstructMPExtenders(const ExtendersGenerator &generator);

    void AddScaffUniqueStorage(size_t uniqe_edge_len);

    Extenders ConstructPBExtenders(const ExtendersGenerator &generator);

    void FilterPaths();

public:

    PathExtendLauncher(const config::dataset& dataset_info,
                       const PathExtendParamsContainer& params,
                       conj_graph_pack& gp):
        dataset_info_(dataset_info),
        params_(params),
        gp_(gp),
        support_(dataset_info, params),
        contig_name_generator_(MakeContigNameGenerator(params_.mode, gp)),
        writer_(gp.g, contig_name_generator_),
        unique_data_()
    {
        unique_data_.min_unique_length_ = params.pset.scaffolding2015.unique_length_upper_bound;
        unique_data_.unique_variation_ = params.pset.uniqueness_analyser.unique_coverage_variation;
    }

    void Launch();

};

}

#endif //PROJECT_LAUNCHER_H
