//
// Created by andrey on 14.11.16.
//

#ifndef PROJECT_LAUNCHER_H
#define PROJECT_LAUNCHER_H

#include "launch_support.hpp"
#include "extenders_logic.hpp"

#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "modules/path_extend/pe_resolver.hpp"
#include "modules/path_extend/pe_io.hpp"
#include "modules/path_extend/path_visualizer.hpp"
#include "modules/path_extend/loop_traverser.hpp"
#include "modules/alignment/long_read_storage.hpp"
#include "modules/path_extend/scaffolder2015/extension_chooser2015.hpp"
#include "modules/genome_consistance_checker.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_visualizer.hpp"
#include "modules/path_extend/scaffolder2015/path_polisher.hpp"
#include "assembly_graph/graph_support/coverage_uniformity_analyzer.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"



namespace path_extend {

using namespace debruijn_graph;

class PathExtendLauncher {

private:
    const config::dataset& dataset_info_;
    const PathExtendParamsContainer& params_;
    const conj_graph_pack& gp_;

    ScaffoldingUniqueEdgeStorage main_unique_storage_;
    vector<ScaffoldingUniqueEdgeStorage*> unique_storages_;

    DefaultContigCorrector<ConjugateDeBruijnGraph> corrector_;
    DefaultContigConstructor<ConjugateDeBruijnGraph> constructor_;
    ContigWriter writer_;

    size_t min_unique_length_;
    double unique_variation_;

    void PrepareForRun() const;

    vector<shared_ptr<ConnectionCondition>>
        ConstructPairedConnectionConditions(const ScaffoldingUniqueEdgeStorage& edge_storage) const;

    shared_ptr<scaffold_graph::ScaffoldGraph>
        ConstructScaffoldGraph(const ScaffoldingUniqueEdgeStorage& edge_storage) const;

    void PrintScaffoldGraph(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                            const set<EdgeId>& main_edge_set,
                            const debruijn_graph::GenomeConsistenceChecker& genome_checker,
                            const string& filename) const;

    void ScaffoldGraphConstruction() const;

    void CountMisassembliesWithReference(const PathContainer& paths) const;

    void FillUniqueEdgeStorage();

    void DebugOutputPaths(const PathContainer& paths, const string& name) const;

    void OutputBrokenScaffolds(const PathContainer& paths, const std::string& filename) const;

    void FinalizePaths(PathContainer& paths, GraphCoverageMap &cover_map, const PathExtendResolver&resolver) const;

    void TraverseLoops(PathContainer& paths, GraphCoverageMap& cover_map) const;

    void PolishPaths(const PathContainer &paths, PathContainer &result) const;

    vector<shared_ptr<PathExtender>> ConstructExtenders(const GraphCoverageMap& cover_map);

public:

    PathExtendLauncher(const config::dataset& dataset_info,
                       const PathExtendParamsContainer& params,
                       const conj_graph_pack& gp):
        dataset_info_(dataset_info),
        params_(params),
        gp_(gp),
        main_unique_storage_(),
        unique_storages_(),
        corrector_(gp.g),
        constructor_(gp.g, corrector_),
        writer_(gp.g, constructor_, gp_.components, params_.mode == config::pipeline_type::plasmid),
        min_unique_length_(params.pset.scaffolding2015.unique_length_upper_bound),
        unique_variation_(params.pset.uniqueness_analyser.unique_coverage_variation)
    { }

    ~PathExtendLauncher() {
        for (auto storage : unique_storages_)
            delete storage;
    }

    void Launch();

};

}

#endif //PROJECT_LAUNCHER_H
