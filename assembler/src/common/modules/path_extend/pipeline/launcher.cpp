//
// Created by andrey on 14.11.16.
//

#include "launcher.hpp"

#include "modules/path_extend/path_visualizer.hpp"
#include "modules/path_extend/loop_traverser.hpp"
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
            conditions.push_back(make_shared<ScaffoldGraphPairedConnectionCondition>(gp_.g, edge_storage.GetSet(),
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
    LengthLowerBound edge_condition(gp_.g, edge_storage.GetMinLength());

    vector<shared_ptr<ConnectionCondition>> conditions =
        ConstructPairedConnectionConditions(edge_storage);

    if (params.use_graph_connectivity) {
        auto as_con = make_shared<AssemblyGraphConnectionCondition>(gp_.g, params.max_path_length, edge_storage);
        as_con->AddInterestingEdges(edge_condition);
        conditions.push_back(as_con);
    }

    INFO("Total conditions " << conditions.size());

    INFO("Constructing scaffold graph from set of size " << edge_storage.GetSet().size());

    DefaultScaffoldGraphConstructor constructor(gp_.g, edge_storage.GetSet(), conditions, edge_condition);
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

    INFO("Constructing reference labels");
    map<debruijn_graph::EdgeId, string> edge_labels;
    size_t count = 0;
    for (const auto &edge_coord_pair: genome_checker.ConstructEdgeOrder()) {
        if (edge_labels.find(edge_coord_pair.first) == edge_labels.end()) {
            edge_labels[edge_coord_pair.first] = "";
        }
        edge_labels[edge_coord_pair.first] += "order: " + ToString(count) +
            "\n mapped range: " + ToString(edge_coord_pair.second.mapped_range.start_pos) + " : "
            + ToString(edge_coord_pair.second.mapped_range.end_pos) +
            "\n init range: " + ToString(edge_coord_pair.second.initial_range.start_pos) + " : "
            + ToString(edge_coord_pair.second.initial_range.end_pos) + "\n";
        ++count;
    }

    auto vertex_colorer = make_shared<ScaffoldVertexSetColorer>(main_edge_set);
    auto edge_colorer = make_shared<ScaffoldEdgeColorer>();
    graph_colorer::CompositeGraphColorer<ScaffoldGraph> colorer(vertex_colorer, edge_colorer);

    INFO("Visualizing scaffold graph");
    ScaffoldGraphVisualizer singleVisualizer(scaffold_graph, edge_labels);
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
        debruijn_graph::GenomeConsistenceChecker genome_checker(gp_, unique_data_.main_unique_storage_,
                                                                params_.pset.genome_consistency_checker.max_gap,
                                                                params_.pset.genome_consistency_checker.relative_max_gap);
        scaffold_graph = ConstructScaffoldGraph(unique_data_.main_unique_storage_);
        if (params_.pset.scaffold_graph_params.output) {
            PrintScaffoldGraph(*scaffold_graph,
                               unique_data_.main_unique_storage_.GetSet(),
                               genome_checker,
                               params_.etc_dir + "scaffold_graph");
        }
    }
}

void PathExtendLauncher::CountMisassembliesWithReference(const PathContainer &paths) const {
    if (gp_.genome.size() == 0)
        return;

    debruijn_graph::GenomeConsistenceChecker genome_checker(gp_, unique_data_.main_unique_storage_,
                                                            params_.pset.genome_consistency_checker.max_gap,
                                                            params_.pset.genome_consistency_checker.relative_max_gap);

    size_t total_mis = 0, gap_mis = 0;
    genome_checker.SpellGenome();
    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        BidirectionalPath *path = iter.get();
        auto map_res = genome_checker.CountMisassemblies(*path);
        if (map_res.misassemblies > 0) {
            INFO ("there are " << map_res.misassemblies << " misassemblies in path: ");
            path->PrintInfo();
            total_mis += map_res.misassemblies;
        }
        if (map_res.wrong_gap_size > 0) {
            INFO ("there are " << map_res.wrong_gap_size << " wrong gaps in path: ");
            path->PrintInfo();
            gap_mis += map_res.wrong_gap_size;
        }
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
    PathInfoWriter path_writer;
    PathVisualizer visualizer;

    writer_.OutputPaths(paths, params_.etc_dir + name);
    if (params_.pe_cfg.output.write_paths) {
        path_writer.WritePaths(paths, params_.etc_dir + name + ".dat");
    }
    if (params_.pe_cfg.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp_, params_.etc_dir + name + ".dot", name, paths);
    }
}

void PathExtendLauncher::FinalizePaths(PathContainer &paths,
                                       GraphCoverageMap &cover_map,
                                       const PathExtendResolver &resolver) const {

    if (params_.pset.remove_overlaps) {
        resolver.RemoveOverlaps(paths, cover_map, params_.min_edge_len, params_.max_path_diff,
                                 params_.pset.cut_all_overlaps,
                                 (params_.mode == config::pipeline_type::moleculo));
    } else if (params_.mode == config::pipeline_type::rna) {
        resolver.RemoveRNAOverlaps(paths, cover_map, params_.min_edge_len, params_.max_path_diff);
    } else {
        resolver.RemoveEqualPaths(paths, cover_map, params_.min_edge_len);
    }

    if (params_.avoid_rc_connections) {
        paths.FilterInterstandBulges();
    }
    paths.FilterEmptyPaths();
    resolver.AddUncoveredEdges(paths, cover_map);

    if (params_.pset.path_filtration.enabled) {
        LengthPathFilter(gp_.g, params_.pset.path_filtration.min_length).filter(paths);;
        IsolatedPathFilter(gp_.g,
                           params_.pset.path_filtration.min_length_for_low_covered,
                           params_.pset.path_filtration.min_coverage).filter(paths);
        IsolatedPathFilter(gp_.g, params_.pset.path_filtration.isolated_min_length).filter(paths);
    }
    paths.SortByLength();
    for (auto &path : paths) {
        path.first->ResetOverlaps();
    }
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

Extenders PathExtendLauncher::ConstructMPExtender(const ExtendersGenerator &generator, size_t uniqe_edge_len) {
    ScaffoldingUniqueEdgeAnalyzer additional_edge_analyzer(gp_, (size_t) uniqe_edge_len, unique_data_.unique_variation_);
    unique_data_.unique_storages_.push_back(make_shared<ScaffoldingUniqueEdgeStorage>());
    additional_edge_analyzer.FillUniqueEdgeStorage(*unique_data_.unique_storages_.back());

    return generator.MakeMPExtenders(*unique_data_.unique_storages_.back());
}

Extenders PathExtendLauncher::ConstructMPExtenders(const ExtendersGenerator &generator) {
    const pe_config::ParamSetT &pset = params_.pset;

    Extenders extenders =  generator.MakeMPExtenders(unique_data_.main_unique_storage_);
    INFO("Using " << extenders.size() << " mate-pair " << support_.LibStr(extenders.size()));

    size_t cur_length = unique_data_.min_unique_length_ - pset.scaffolding2015.unique_length_step;
    size_t lower_bound = max(pset.scaffolding2015.unique_length_lower_bound, pset.scaffolding2015.unique_length_step);

    while (cur_length > lower_bound) {
        INFO("Adding extender with length " << cur_length);
        push_back_all(extenders, ConstructMPExtender(generator, cur_length));
        cur_length -= pset.scaffolding2015.unique_length_step;
    }
    if (unique_data_.min_unique_length_ > lower_bound) {
        INFO("Adding final extender with length " << lower_bound);
        push_back_all(extenders, ConstructMPExtender(generator, lower_bound));
    }

    return extenders;
}

void PathExtendLauncher::FillPathContainer(size_t lib_index, size_t size_threshold) {
    std::vector<PathInfo<Graph>> paths;
    gp_.single_long_reads[lib_index].SaveAllPaths(paths);
    for (const auto &path: paths) {
        const auto &edges = path.getPath();
        if (edges.size() <= size_threshold)
            continue;

        BidirectionalPath *new_path = new BidirectionalPath(gp_.g, edges);
        BidirectionalPath *conj_path = new BidirectionalPath(new_path->Conjugate());
        new_path->SetWeight((float) path.getWeight());
        conj_path->SetWeight((float) path.getWeight());
        unique_data_.long_reads_paths_[lib_index]->AddPair(new_path, conj_path);
    }
    DEBUG("Long reads paths " << unique_data_.long_reads_paths_[lib_index]->size());
    unique_data_.long_reads_cov_map_[lib_index]->AddPaths(*unique_data_.long_reads_paths_[lib_index]);
}


void PathExtendLauncher::FillLongReadsCoverageMaps() {
    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
        unique_data_.long_reads_paths_.push_back(make_shared<PathContainer>());
        unique_data_.long_reads_cov_map_.push_back(make_shared<GraphCoverageMap>(gp_.g));
        if (support_.IsForSingleReadExtender(dataset_info_.reads[lib_index])) {
            FillPathContainer(lib_index);
        }
    }
}

void  PathExtendLauncher::FillPBUniqueEdgeStorages() {
    //FIXME magic constants
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
    return generator.MakePBScaffoldingExtenders(unique_data_.unique_pb_storage_,
                                                unique_data_.long_reads_cov_map_);
}


Extenders PathExtendLauncher::ConstructExtenders(const GraphCoverageMap& cover_map) {
    INFO("Creating main extenders, unique edge length = " << unique_data_.min_unique_length_);
    if (support_.SingleReadsMapped() || support_.HasLongReads())
        FillLongReadsCoverageMaps();

    ExtendersGenerator generator(dataset_info_, params_, gp_, cover_map, support_);
    Extenders extenders = generator.MakeBasicExtenders(unique_data_.main_unique_storage_,
                                                       unique_data_.long_reads_cov_map_);

    //long reads scaffolding extenders.
    if (support_.HasLongReads()) {
        if (params_.pset.sm == sm_old) {
            INFO("Will not use new long read scaffolding algorithm in this mode");
        } else {
            push_back_all(extenders, ConstructPBExtenders(generator));
        }
    }

    if (support_.HasMPReads()) {
        if (params_.pset.sm == sm_old) {
            INFO("Will not use mate-pairs is this mode");
        } else {
            push_back_all(extenders, ConstructMPExtenders(generator));
        }
    }

    if (params_.pset.use_coordinated_coverage)
        push_back_all(extenders, generator.MakeCoverageExtenders());

    INFO("Total number of extenders is " << extenders.size());
    return extenders;
}

void PathExtendLauncher::PolishPaths(const PathContainer &paths, PathContainer &result) const {
    //Fixes distances for paths gaps and tries to fill them in
    INFO("Closing gaps in paths");
    PathPolisher polisher(gp_, dataset_info_, unique_data_.main_unique_storage_, params_.max_polisher_gap);
    polisher.PolishPaths(paths, result);
    result.SortByLength();
    INFO("Gap closing completed")
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
    Extenders extenders = ConstructExtenders(cover_map);
    shared_ptr<CompositeExtender> composite_extender = make_shared<CompositeExtender>(gp_.g, cover_map, extenders,
                                                                                      unique_data_.main_unique_storage_,
                                                                                      params_.max_path_diff,
                                                                                      params_.pset.extension_options.max_repeat_length,
                                                                                      params_.detect_repeats_online);

    auto paths = resolver.ExtendSeeds(seeds, *composite_extender);
    paths.FilterEmptyPaths();
    paths.SortByLength();
    DebugOutputPaths(paths, "raw_paths");

    FinalizePaths(paths, cover_map, resolver);
    DebugOutputPaths(paths, "before_loop_traversal");

    TraverseLoops(paths, cover_map);
    DebugOutputPaths(paths, "loop_traveresed");

    PolishPaths(paths, gp_.contig_paths);
    DebugOutputPaths(gp_.contig_paths, "polished_paths");
    
    GraphCoverageMap polished_map(gp_.g, gp_.contig_paths, true);
    FinalizePaths(gp_.contig_paths, polished_map, resolver);
    DebugOutputPaths(gp_.contig_paths, "final_paths");

    CountMisassembliesWithReference(gp_.contig_paths);

    INFO("ExSPAnder repeat resolving tool finished");
}

}
