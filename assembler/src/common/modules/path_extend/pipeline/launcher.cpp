//
// Created by andrey on 14.11.16.
//

#include "launcher.hpp"

#include "modules/path_extend/path_visualizer.hpp"
#include "modules/path_extend/loop_traverser.hpp"
#include "modules/path_extend/path_extender.hpp"
#include "modules/alignment/long_read_storage.hpp"
#include "modules/path_extend/scaffolder2015/extension_chooser2015.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_visualizer.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "assembly_graph/graph_support/coverage_uniformity_analyzer.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "modules/path_extend/scaffolder2015/path_polisher.hpp"


namespace path_extend {

using namespace debruijn_graph;
using namespace std;

vector<shared_ptr<ConnectionCondition>>
    PathExtendLauncher::ConstructPairedConnectionConditions(const ScaffoldingUniqueEdgeStorage& edge_storage) const {

    vector<shared_ptr<ConnectionCondition>> conditions;
    const pe_config::ParamSetT::ScaffoldGraphParamsT &params = params_.pset.scaffold_graph_params;

    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info_.reads[lib_index];
        if (lib.is_paired()) {
            shared_ptr<PairedInfoLibrary> paired_lib;
            if (lib.is_mate_pair())
                paired_lib = MakeNewLib(gp_.g, lib, gp_.paired_indices[lib_index]);
            else if (lib.type() == io::LibraryType::PairedEnd)
                paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
            else {
                INFO("Unusable for scaffold graph paired lib #" << lib_index);
                continue;
            }
            conditions.push_back(make_shared<ScaffoldGraphPairedConnectionCondition>(gp_.g, edge_storage.unique_edges(),
                                                                                     paired_lib, lib_index,
                                                                                     params.always_add,
                                                                                     params.never_add,
                                                                                     params.relative_threshold));
        }
    }
    return conditions;
}

shared_ptr<scaffold_graph::ScaffoldGraph> PathExtendLauncher::ConstructScaffoldGraph(const ScaffoldingUniqueEdgeStorage &edge_storage) const {
    using namespace scaffold_graph;

    const pe_config::ParamSetT::ScaffoldGraphParamsT &params = params_.pset.scaffold_graph_params;

    INFO("Constructing connections");
    LengthLowerBound edge_condition(gp_.g, edge_storage.min_length());

    vector<shared_ptr<ConnectionCondition>> conditions =
        ConstructPairedConnectionConditions(edge_storage);

    if (params.use_graph_connectivity) {
        auto as_con = make_shared<AssemblyGraphConnectionCondition>(gp_.g, params.max_path_length, edge_storage);
        as_con->AddInterestingEdges(edge_condition);
        conditions.push_back(as_con);
    }

    INFO("Total conditions " << conditions.size());

    INFO("Constructing scaffold graph from set of size " << edge_storage.unique_edges().size());

    DefaultScaffoldGraphConstructor constructor(gp_.g, edge_storage.unique_edges(), conditions, edge_condition);
    auto scaffold_graph = constructor.Construct();

    INFO("Scaffold graph contains " << scaffold_graph->VertexCount() << " vertices and " << scaffold_graph->EdgeCount()
             << " edges");
    return scaffold_graph;
}

void PathExtendLauncher::PrintScaffoldGraph(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                                            const set<EdgeId> &main_edge_set,
                                            const debruijn_graph::GenomeConsistenceChecker &genome_checker,
                                            const string &filename) const {
    using namespace scaffold_graph;

    auto vertex_colorer = make_shared<ScaffoldVertexSetColorer>(main_edge_set);
    auto edge_colorer = make_shared<ScaffoldEdgeColorer>();
    graph_colorer::CompositeGraphColorer<ScaffoldGraph> colorer(vertex_colorer, edge_colorer);

    INFO("Visualizing scaffold graph");
    ScaffoldGraphVisualizer singleVisualizer(scaffold_graph, genome_checker.EdgeLabels());
    std::ofstream single_dot;
    single_dot.open((filename + "_single.dot").c_str());
    singleVisualizer.Visualize(single_dot, colorer);
    single_dot.close();

    INFO("Printing scaffold graph");
    std::ofstream data_stream;
    data_stream.open((filename + ".data").c_str());
    scaffold_graph.Print(data_stream);
    data_stream.close();
}


void PathExtendLauncher::MakeAndOutputScaffoldGraph() const {
    //Scaffold graph
    shared_ptr<scaffold_graph::ScaffoldGraph> scaffold_graph;
    if (params_.pset.scaffold_graph_params.construct) {
        debruijn_graph::GenomeConsistenceChecker genome_checker(gp_,
                                                                params_.pset.genome_consistency_checker.max_gap,
                                                                params_.pset.genome_consistency_checker.relative_max_gap,
                                                                unique_data_.main_unique_storage_.min_length(),
                                                                unique_data_.main_unique_storage_,
                                                                unique_data_.long_reads_cov_map_,
                                                                dataset_info_.reads);
        scaffold_graph = ConstructScaffoldGraph(unique_data_.main_unique_storage_);
        if (params_.pset.scaffold_graph_params.output) {
            PrintScaffoldGraph(*scaffold_graph,
                               unique_data_.main_unique_storage_.unique_edges(),
                               genome_checker,
                               params_.etc_dir + "scaffold_graph");
        }
    }
}

void PathExtendLauncher::CountMisassembliesWithReference(const PathContainer &paths) const {
    if (gp_.genome.size() == 0)
        return;
    bool use_main_storage = params_.pset.genome_consistency_checker.use_main_storage;
    size_t unresolvable_gap = unique_data_.main_unique_storage_.min_length();
    ScaffoldingUniqueEdgeStorage tmp_storage;
    if (!use_main_storage) {
        unresolvable_gap = params_.pset.genome_consistency_checker.unresolvable_jump;
        ScaffoldingUniqueEdgeAnalyzer tmp_analyzer(gp_, params_.pset.genome_consistency_checker.unique_length, unique_data_.unique_variation_);
        tmp_analyzer.FillUniqueEdgeStorage(tmp_storage);
    }
    debruijn_graph::GenomeConsistenceChecker genome_checker(gp_,
                                                            params_.pset.genome_consistency_checker.max_gap,
                                                            params_.pset.genome_consistency_checker.relative_max_gap,
                                                            unresolvable_gap,
                                                            use_main_storage ? unique_data_.main_unique_storage_ : tmp_storage,
                                                            unique_data_.long_reads_cov_map_,
                                                            dataset_info_.reads);

    size_t total_mis = 0, gap_mis = 0;
    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        BidirectionalPath *path = iter.get();
        auto map_res = genome_checker.CountMisassemblies(*path);
        if (map_res.misassemblies > 0) {
            INFO ("there are " << map_res.misassemblies << " misassemblies in path: ");
            path->PrintINFO();
            total_mis += map_res.misassemblies;
        }
        if (map_res.wrong_gap_size > 0) {
            INFO ("there are " << map_res.wrong_gap_size << " wrong gaps in path. ");
            path->PrintDEBUG();
            gap_mis += map_res.wrong_gap_size;
        }
        genome_checker.CheckPathEnd(*path);
        genome_checker.CheckPathEnd(path->Conjugate());
    }
    INFO ("In total found " << total_mis << " misassemblies " << " and " << gap_mis << " gaps.");
}


void PathExtendLauncher::EstimateUniqueEdgesParams() {
    bool uniform_coverage = false;
    if (params_.pset.uniqueness_analyser.enabled) {
        INFO("Autodetecting unique edge set parameters...");
        unique_data_.min_unique_length_ = max(unique_data_.min_unique_length_, support_.FindMaxMPIS());
        INFO("Minimal unique edge length set to the smallest MP library IS: " << unique_data_.min_unique_length_);

        CoverageUniformityAnalyzer coverage_analyzer(gp_.g, unique_data_.min_unique_length_);
        double median_coverage = coverage_analyzer.CountMedianCoverage();
        double uniformity_fraction = coverage_analyzer.UniformityFraction(unique_data_.unique_variation_, median_coverage);
        INFO ("median coverage for edges longer than " << unique_data_.min_unique_length_ << " is " << median_coverage <<
            " uniformity " << size_t(uniformity_fraction * 100) << "%");
        if (math::gr(uniformity_fraction, params_.pset.uniqueness_analyser.uniformity_fraction_threshold)) {
            uniform_coverage = true;
        }
        if (!uniform_coverage) {
            unique_data_.unique_variation_ = params_.pset.uniqueness_analyser.nonuniform_coverage_variation;
            INFO("Coverage is not uniform, we do not rely on coverage for long edge uniqueness");
        }

    } else {
        INFO("Unique edge set constructed with parameters from config : length " << unique_data_.min_unique_length_
                 << " variation " << unique_data_.unique_variation_);
    }
}


void PathExtendLauncher::FillUniqueEdgeStorage() {
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, unique_data_.min_unique_length_, unique_data_.unique_variation_);
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_data_.main_unique_storage_);
}

void PathExtendLauncher::DebugOutputPaths(const PathContainer &paths, const string &name) const {
    if (!params_.pe_cfg.debug_output) {
        return;
    }
    PathVisualizer visualizer;

    writer_.OutputPaths(paths, params_.etc_dir + name + ".fasta");
    if (params_.pe_cfg.output.write_paths) {
        std::ofstream oss(params_.etc_dir + name + ".dat");
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            iter.get()->Print(oss);
        }
        oss.close();
    }
    if (params_.pe_cfg.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp_, params_.etc_dir + name + ".dot", name, paths);
    }
}

void FilterInterstandBulges(PathContainer &paths) {
    DEBUG ("Try to delete paths with interstand bulges");
    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        if (EndsWithInterstrandBulge(*iter.get())) {
            iter.get()->PopBack();
        }
        if (EndsWithInterstrandBulge(*iter.getConjugate())) {
            iter.getConjugate()->PopBack();
        }
    }
    DEBUG("deleted paths with interstand bulges");
}

void PathExtendLauncher::RemoveOverlapsAndArtifacts(PathContainer &paths,
                                                    GraphCoverageMap &cover_map,
                                                    const PathExtendResolver &resolver) const {
    INFO("Finalizing paths");

    INFO("Deduplicating paths");
    Deduplicate(gp_.g, paths, cover_map, params_.min_edge_len,
                         params_.max_path_diff);

    INFO("Paths deduplicated");

    if (params_.pset.overlap_removal.enabled) {
        resolver.RemoveOverlaps(paths, cover_map, params_.min_edge_len, params_.max_path_diff,
                                params_.pset.overlap_removal.end_start_only,
                                params_.pset.overlap_removal.cut_all);
    } else {
        INFO("Overlaps will not be removed");
    }

    //TODO do we still need it?
    if (params_.avoid_rc_connections) {
        FilterInterstandBulges(paths);
    }
    resolver.AddUncoveredEdges(paths, cover_map);

    paths.SortByLength();
    INFO("Paths finalized");
}


void PathExtendLauncher::CleanPaths(PathContainer &paths, const pe_config::ParamSetT::PathFiltrationT &path_filtration) const {
    if (path_filtration.enabled) {
        paths.FilterPaths(LengthPathCondition(GetLengthCutoff(path_filtration.min_length, path_filtration.rel_cutoff)));
        paths.FilterPaths(func::And(CoveragePathCondition(gp_.g, path_filtration.min_coverage),
                                    LengthPathCondition(GetLengthCutoff(path_filtration.min_length_for_low_covered, path_filtration.rel_low_covered_cutoff))));
        paths.FilterPaths(func::And(IsolatedPathCondition(gp_.g),
                                    func::And(LengthPathCondition(GetLengthCutoff(path_filtration.isolated_min_length, path_filtration.rel_isolated_cutoff)),
                                              CoveragePathCondition(gp_.g, path_filtration.isolated_min_cov))));
    }

    paths.SortByLength();
}

size_t PathExtendLauncher::GetLengthCutoff(size_t abs_cutoff, double rel_cutoff) const {
    int rel_len = int(rel_cutoff * double(cfg::get().ds.RL)) - int(cfg::get().K);
    int abs_len = int(abs_cutoff) - int(cfg::get().K);
    size_t result = (size_t) max(0, max(rel_len, abs_len));

    INFO("Read length relative cutoff " << rel_cutoff << " converted to " << rel_len);
    INFO("Read length absolute cutoff " << abs_cutoff << " bp converted to " << result);
    INFO("Length cutoff: " << result);
    return result;
}

void PathExtendLauncher::TraverseLoops(PathContainer &paths, GraphCoverageMap &cover_map) const {
    INFO("Traversing tandem repeats");

    LoopTraverser
        loopTraverser(cover_map.graph(), cover_map,
                      params_.pset.loop_traversal.min_edge_length,
                      params_.pset.loop_traversal.max_component_size,
                      params_.pset.loop_traversal.max_path_length);
    size_t res = loopTraverser.TraverseAllLoops();
    paths.SortByLength();

    INFO("Traversed " << res << " loops");
}

void PathExtendLauncher::AddScaffUniqueStorage(size_t uniqe_edge_len) {
    ScaffoldingUniqueEdgeAnalyzer additional_edge_analyzer(gp_, (size_t) uniqe_edge_len,
                                                           unique_data_.unique_variation_);
    unique_data_.unique_storages_.push_back(ScaffoldingUniqueEdgeStorage());
    additional_edge_analyzer.FillUniqueEdgeStorage(unique_data_.unique_storages_.back());
}

Extenders PathExtendLauncher::ConstructMPExtenders(const ExtendersGenerator &generator) {
    const pe_config::ParamSetT &pset = params_.pset;

    size_t cur_length = unique_data_.min_unique_length_ - pset.scaffolding2015.unique_length_step;
    size_t lower_bound = max(pset.scaffolding2015.unique_length_lower_bound, pset.scaffolding2015.unique_length_step);

    while (cur_length > lower_bound) {
        INFO("Will add extenders for length " << cur_length);
        AddScaffUniqueStorage(cur_length);
        cur_length -= pset.scaffolding2015.unique_length_step;
    }
    if (unique_data_.min_unique_length_ > lower_bound) {
        INFO("Will add final extenders for length " << lower_bound);
        AddScaffUniqueStorage(lower_bound);
    }

    return generator.MakeMPExtenders();
}

void PathExtendLauncher::FillPathContainer(size_t lib_index, size_t size_threshold) {
    std::vector<PathInfo<Graph>> paths;
    gp_.single_long_reads[lib_index].SaveAllPaths(paths);
    for (const auto &path: paths) {
        const auto &edges = path.path();
        if (edges.size() <= size_threshold)
            continue;

        BidirectionalPath *new_path = new BidirectionalPath(gp_.g, edges);
        BidirectionalPath *conj_path = new BidirectionalPath(new_path->Conjugate());
        new_path->SetWeight((float) path.weight());
        conj_path->SetWeight((float) path.weight());
        unique_data_.long_reads_paths_[lib_index].AddPair(new_path, conj_path);
    }
    DEBUG("Long reads paths " << unique_data_.long_reads_paths_[lib_index].size());
    unique_data_.long_reads_cov_map_[lib_index].AddPaths(unique_data_.long_reads_paths_[lib_index]);
}


void PathExtendLauncher::FillLongReadsCoverageMaps() {
    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
        unique_data_.long_reads_paths_.push_back(PathContainer());
        unique_data_.long_reads_cov_map_.push_back(GraphCoverageMap(gp_.g));
        if (support_.IsForSingleReadExtender(dataset_info_.reads[lib_index])) {
            FillPathContainer(lib_index);
        }
    }
}

void  PathExtendLauncher::FillPBUniqueEdgeStorages() {
    //FIXME magic constants
    //FIXME need to change for correct usage of prelimnary contigs in loops
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer_pb(gp_, 500, 0.5);

    INFO("Filling backbone edges for long reads scaffolding...");
    if (params_.uneven_depth) {
        INFO(" with long reads paths");
        //TODO:: muiltiple libraries?
        for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
            if (support_.IsForSingleReadScaffolder(dataset_info_.reads[lib_index])) {
                unique_edge_analyzer_pb.FillUniqueEdgesWithLongReads(unique_data_.long_reads_cov_map_[lib_index],
                                                                     unique_data_.unique_pb_storage_,
                                                                     support_.GetLongReadsConfig(dataset_info_.reads[lib_index].type()));
            }
        }
        INFO("Removing fake unique with paired-end libs");
        for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
            if (dataset_info_.reads[lib_index].type() == io::LibraryType::PairedEnd) {
                unique_edge_analyzer_pb.ClearLongEdgesWithPairedLib(lib_index, unique_data_.unique_pb_storage_);
            }
        }

    } else {
        INFO(" with coverage")
        unique_edge_analyzer_pb.FillUniqueEdgeStorage(unique_data_.unique_pb_storage_);
    }
    INFO(unique_data_.unique_pb_storage_.size() << " unique edges");
}

Extenders PathExtendLauncher::ConstructPBExtenders(const ExtendersGenerator &generator) {
    FillPBUniqueEdgeStorages();
    return generator.MakePBScaffoldingExtenders();
}


Extenders PathExtendLauncher::ConstructExtenders(const GraphCoverageMap &cover_map,
                                                 UsedUniqueStorage &used_unique_storage) {
    INFO("Creating main extenders, unique edge length = " << unique_data_.min_unique_length_);
    if (support_.SingleReadsMapped() || support_.HasLongReads())
        FillLongReadsCoverageMaps();

    ExtendersGenerator generator(dataset_info_, params_, gp_, cover_map,
                                 unique_data_, used_unique_storage, support_);
    Extenders extenders = generator.MakeBasicExtenders();

    //long reads scaffolding extenders.
    if (support_.HasLongReads()) {
        if (params_.pset.sm == sm_old) {
            INFO("Will not use new long read scaffolding algorithm in this mode");
        } else {
            utils::push_back_all(extenders, ConstructPBExtenders(generator));
        }
    }

    if (support_.HasMPReads()) {
        if (params_.pset.sm == sm_old) {
            INFO("Will not use mate-pairs is this mode");
        } else {
            utils::push_back_all(extenders, ConstructMPExtenders(generator));
        }
    }

    if (params_.pset.use_coordinated_coverage)
        utils::push_back_all(extenders, generator.MakeCoverageExtenders());

    INFO("Total number of extenders is " << extenders.size());
    return extenders;
}

void PathExtendLauncher::PolishPaths(const PathContainer &paths, PathContainer &result,
                                     const GraphCoverageMap& /* cover_map */) const {
    //Fixes distances for paths gaps and tries to fill them in
    INFO("Closing gaps in paths");

    vector<shared_ptr<PathGapCloser>> gap_closers;

    gap_closers.push_back(make_shared<DijkstraGapCloser>(gp_.g, params_.max_polisher_gap));
    for (size_t i = 0; i < dataset_info_.reads.lib_count(); i++) {
        auto lib = dataset_info_.reads[i];
        if (lib.type() == io::LibraryType::HQMatePairs || lib.type() == io::LibraryType::MatePairs) {
            shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.paired_indices[i]);
            gap_closers.push_back(make_shared<MatePairGapCloser> (gp_.g, params_.max_polisher_gap, paired_lib,
                                                                   unique_data_.main_unique_storage_));
        }
    }

////TODO:: is it really empty?
//    UniqueData unique_data;
//    UsedUniqueStorage used_unique_storage(unique_data.main_unique_storage_);
//    ExtendersGenerator generator(dataset_info_, params_, gp_, cover_map,
//                                 unique_data, used_unique_storage, support_);
//    auto polisher_storage = ScaffoldingUniqueEdgeStorage();
//    for  (const auto& extender: generator.MakePEExtenders()) {
//        gap_closers.push_back(make_shared<PathExtenderGapCloser>(gp_.g, params_.max_polisher_gap, extender));
//    }
//FIXME: uncomment cover_map 

    PathPolisher polisher(gp_, gap_closers);
    result = polisher.PolishPaths(paths);
    result.SortByLength();
    INFO("Gap closing completed")
}

void PathExtendLauncher::FilterPaths() {
    PathContainer contig_paths_copy(gp_.contig_paths.begin(), gp_.contig_paths.end());
    for (const auto& it: params_.pset.path_filtration) {
        if (it.first == "default" && it.second.enabled) {
            INFO("Finalizing main paths");
            CleanPaths(gp_.contig_paths, it.second);
            DebugOutputPaths(gp_.contig_paths, "final_paths");
        }
        else if (it.second.enabled) {
            INFO("Finalizing paths - " + it.first);
            PathContainer to_clean(contig_paths_copy.begin(), contig_paths_copy.end());
            CleanPaths(to_clean, it.second);
            DebugOutputPaths(to_clean, it.first + "_final_paths");
            writer_.OutputPaths(to_clean, params_.output_dir + it.first + "_filtered_final_paths" + ".fasta");
        }
    }
}

void PathExtendLauncher::Launch() {
    INFO("ExSPAnder repeat resolving tool started");
    make_dir(params_.output_dir);
    make_dir(params_.etc_dir);

    if (support_.NeedsUniqueEdgeStorage()) {
        //Fill the storage to enable unique edge check
        EstimateUniqueEdgesParams();
        FillUniqueEdgeStorage();
    }

    MakeAndOutputScaffoldGraph();

    PathExtendResolver resolver(gp_.g);

    auto seeds = resolver.MakeSimpleSeeds();
    seeds.SortByLength();
    DebugOutputPaths(seeds, "init_paths");

    GraphCoverageMap cover_map(gp_.g);
    UsedUniqueStorage used_unique_storage(unique_data_.main_unique_storage_);
    Extenders extenders = ConstructExtenders(cover_map, used_unique_storage);
    CompositeExtender composite_extender(gp_.g, cover_map,
                                         used_unique_storage,
                                         extenders);

    auto paths = resolver.ExtendSeeds(seeds, composite_extender);
    DebugOutputPaths(paths, "raw_paths");

    RemoveOverlapsAndArtifacts(paths, cover_map, resolver);
    DebugOutputPaths(paths, "before_path_polishing");

    //TODO does path polishing correctly work with coverage map
    PolishPaths(paths, gp_.contig_paths, cover_map);
    //TODO use move assignment to original map here
    GraphCoverageMap polished_map(gp_.g, gp_.contig_paths, true);
    DebugOutputPaths(gp_.contig_paths, "polished_paths");

    TraverseLoops(gp_.contig_paths, polished_map);
    DebugOutputPaths(gp_.contig_paths, "loop_traveresed");

    RemoveOverlapsAndArtifacts(gp_.contig_paths, polished_map, resolver);
    DebugOutputPaths(gp_.contig_paths, "overlap_removed");

    if (params_.ss.ss_enabled) {
        PathContainerCoverageSwitcher switcher(gp_.g, gp_.ss_coverage.front(), params_.ss.antisense);
        switcher.Apply(gp_.contig_paths);
    }

    FilterPaths();

    CountMisassembliesWithReference(gp_.contig_paths);

    INFO("ExSPAnder repeat resolving tool finished");
}

}
