//
// Created by andrey on 14.11.16.
//

#ifndef PROJECT_LAUNCHER_H
#define PROJECT_LAUNCHER_H

#include "extenders_logic.hpp"
#include "launch_support.hpp"

#include "modules/path_extend/pe_resolver.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "modules/genome_consistance_checker.hpp"

#include "alignment/rna/ss_coverage.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

namespace path_extend {

using namespace debruijn_graph;

class PathExtendLauncher {
    const config::dataset& dataset_info_;
    const PathExtendParamsContainer& params_;
    graph_pack::GraphPack& gp_;
    const Graph &graph_;
    PELaunchSupport support_;

    std::shared_ptr<ContigNameGenerator> contig_name_generator_;
    ContigWriter writer_;

    UniqueData unique_data_;

    std::vector<std::shared_ptr<ConnectionCondition>>
        ConstructPairedConnectionConditions(const ScaffoldingUniqueEdgeStorage &edge_storage) const;

    std::shared_ptr<scaffold_graph::ScaffoldGraph>
        ConstructScaffoldGraph(const ScaffoldingUniqueEdgeStorage &edge_storage) const;

    void PrintScaffoldGraph(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                            const std::set<EdgeId> &main_edge_set,
                            const debruijn_graph::GenomeConsistenceChecker &genome_checker,
                            const std::filesystem::path &filename) const;

    void MakeAndOutputScaffoldGraph() const;

    void CountMisassembliesWithReference(const PathContainer& paths) const;

    void EstimateUniqueEdgesParams();

    void CheckCoverageUniformity();

    void FillUniqueEdgeStorage();

    void FillPBUniqueEdgeStorages();

    void FillReadCloudUniqueEdgeStorage();

    void FillPathContainer(size_t lib_index, size_t size_threshold = 1);

    void FillLongReadsCoverageMaps();

    void DebugOutputPaths(const PathContainer &paths, const std::string &name) const;

    void RemoveOverlapsAndArtifacts(PathContainer &paths, GraphCoverageMap &cover_map, const PathExtendResolver &resolver) const;

    void CleanPaths(PathContainer &paths, const pe_config::ParamSetT::PathFiltrationT &path_filtration) const;

    size_t GetLengthCutoff(size_t abs_cutoff, double rel_cutoff) const;

    void TraverseLoops(PathContainer &paths, GraphCoverageMap &cover_map) const;

    void PolishPaths(const PathContainer &paths, PathContainer &result, const GraphCoverageMap &cover_map) const;

    Extenders ConstructExtenders(const GraphCoverageMap &cover_map, UsedUniqueStorage &used_unique_storage);

    Extenders ConstructMPExtenders(const ExtendersGenerator &generator);

    void AddScaffUniqueStorage(size_t uniqe_edge_len);

    Extenders ConstructPBExtenders(const ExtendersGenerator &generator);

    Extenders ConstructReadCloudExtenders(const ExtendersGenerator &generator);

    void FilterPaths(PathContainer& paths);

    void AddFLPaths(PathContainer& paths) const;

    void SelectStrandSpecificPaths(PathContainer &paths) const;

public:

    PathExtendLauncher(const config::dataset& dataset_info,
                       const PathExtendParamsContainer& params,
                       graph_pack::GraphPack& gp):
        dataset_info_(dataset_info),
        params_(params),
        gp_(gp),
        graph_(gp.get<Graph>()),
        support_(dataset_info, params),
        contig_name_generator_(MakeContigNameGenerator(params_.mode, gp)),
        writer_(graph_, contig_name_generator_),
        unique_data_() {
        unique_data_.min_unique_length_ = params.pset.scaffolding2015.unique_length_upper_bound;
        unique_data_.unique_variation_ = params.pset.uniqueness_analyser.unique_coverage_variation;
    }

    void Launch();

};

}

#endif //PROJECT_LAUNCHER_H
