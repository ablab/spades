//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2020-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "launcher.hpp"

#include "alignment/long_read_storage.hpp"
#include "alignment/rna/ss_coverage.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "assembly_graph/graph_support/coverage_uniformity_analyzer.hpp"

#include "modules/path_extend/loop_traverser.hpp"
#include "modules/path_extend/path_extender.hpp"
#include "modules/path_extend/path_visualizer.hpp"
#include "modules/path_extend/read_cloud_path_extend/extension_chooser_checker.hpp"
#include "modules/path_extend/read_cloud_path_extend/path_scaffolder.hpp"
#include "modules/path_extend/scaff_supplementary.hpp"
#include "modules/path_extend/scaffolder2015/path_polisher.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_visualizer.hpp"

#include <unordered_set>
namespace path_extend {

using namespace debruijn_graph;
using namespace omnigraph::de;

std::vector<std::shared_ptr<ConnectionCondition>>
PathExtendLauncher::ConstructPairedConnectionConditions(const ScaffoldingUniqueEdgeStorage& edge_storage) const {
    TIME_TRACE_SCOPE("ConstructPairedConnectionConditions");

    std::vector<std::shared_ptr<ConnectionCondition>> conditions;
    const pe_config::ParamSetT::ScaffoldGraphParamsT &params = params_.pset.scaffold_graph_params;

    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info_.reads[lib_index];
        if (lib.is_paired()) {
            std::shared_ptr<PairedInfoLibrary> paired_lib;
            if (lib.is_mate_pair())
                paired_lib = MakeNewLib(graph_, lib, gp_.get<UnclusteredPairedInfoIndicesT<Graph>>()[lib_index]);
            else if (lib.is_paired())
                paired_lib = MakeNewLib(graph_, lib, gp_.get<PairedInfoIndicesT<Graph>>("clustered_indices")[lib_index]);
            else {
                INFO("Unusable for scaffold graph paired lib #" << lib_index);
                continue;
            }
            conditions.push_back(std::make_shared<ScaffoldGraphPairedConnectionCondition>(graph_, edge_storage.unique_edges(),
                                                                                          paired_lib, lib_index,
                                                                                          params.always_add,
                                                                                          params.never_add,
                                                                                          params.relative_threshold));
        }
    }
    return conditions;
}

std::shared_ptr<scaffold_graph::ScaffoldGraph> PathExtendLauncher::ConstructScaffoldGraph(
        const ScaffoldingUniqueEdgeStorage &edge_storage) const {
    using namespace scaffold_graph;

    const pe_config::ParamSetT::ScaffoldGraphParamsT &params = params_.pset.scaffold_graph_params;

    INFO("Constructing connections");
    LengthLowerBound edge_condition(graph_, edge_storage.min_length());

    auto conditions = ConstructPairedConnectionConditions(edge_storage);

    if (params.use_graph_connectivity) {
        auto as_con = std::make_shared<AssemblyGraphConnectionCondition>(graph_, params.max_path_length, edge_storage);
        as_con->AddInterestingEdges(edge_condition);
        conditions.push_back(as_con);
    }

    INFO("Total conditions " << conditions.size());

    INFO("Constructing scaffold graph from set of size " << edge_storage.unique_edges().size());

    DefaultScaffoldGraphConstructor constructor(graph_, edge_storage.unique_edges(), conditions, edge_condition);
    auto scaffold_graph = constructor.Construct();

    INFO("Scaffold graph contains " << scaffold_graph->VertexCount() << " vertices and " << scaffold_graph->EdgeCount()
             << " edges");
    return scaffold_graph;
}

void PathExtendLauncher::PrintScaffoldGraph(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                                            const std::set<EdgeId> &main_edge_set,
                                            const debruijn_graph::GenomeConsistenceChecker &genome_checker,
                                            const std::filesystem::path &filename) const {
    using namespace scaffold_graph;

    set<ScaffoldVertex> scaff_vertex_set;
    for (const auto& edge: main_edge_set) {
        EdgeId copy = edge;
        scaff_vertex_set.insert(copy);
    }
    auto vertex_colorer = std::make_shared<ScaffoldVertexSetColorer>(scaff_vertex_set);
    auto edge_colorer = std::make_shared<ScaffoldEdgeColorer>();
    graph_colorer::CompositeGraphColorer<ScaffoldGraph> colorer(vertex_colorer, edge_colorer);

    INFO("Visualizing scaffold graph");
    map<ScaffoldVertex, string> scaff_vertex_labels;
    for (const auto& entry: genome_checker.EdgeLabels()) {
        scaff_vertex_labels.insert({entry.first, entry.second});
    }
    ScaffoldGraphVisualizer singleVisualizer(scaffold_graph, scaff_vertex_labels);
    std::ofstream single_dot;
    single_dot.open(filename.native() + "_single.dot");
    singleVisualizer.Visualize(single_dot, colorer);
    single_dot.close();

    INFO("Printing scaffold graph");
    std::ofstream data_stream;
    data_stream.open(filename.native() + ".data");
    scaffold_graph.Print(data_stream);
    data_stream.close();
}


void PathExtendLauncher::MakeAndOutputScaffoldGraph() const {
    TIME_TRACE_SCOPE("MakeAndOutputScaffoldGraph");

    auto scaffold_graph = ConstructScaffoldGraph(unique_data_.main_unique_storage_);
    if (params_.pset.scaffold_graph_params.output) {
        debruijn_graph::GenomeConsistenceChecker genome_checker(gp_,
                                                                params_.pset.genome_consistency_checker.max_gap,
                                                                params_.pset.genome_consistency_checker.relative_max_gap,
                                                                unique_data_.main_unique_storage_.min_length(),
                                                                unique_data_.main_unique_storage_,
                                                                unique_data_.long_reads_cov_map_,
                                                                dataset_info_.reads);
        PrintScaffoldGraph(*scaffold_graph,
                           unique_data_.main_unique_storage_.unique_edges(),
                           genome_checker,
                           params_.etc_dir / "scaffold_graph");
    }
}

void PathExtendLauncher::CountMisassembliesWithReference(const PathContainer &paths) const {
    if (!gp_.get<GenomeStorage>().size())
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
        const BidirectionalPath &path = iter.get();
        VERIFY(path.GetConjPath());
        auto map_res = genome_checker.CountMisassemblies(path);
        if (map_res.misassemblies > 0) {
            INFO ("there are " << map_res.misassemblies << " misassemblies in path: ");
            path.PrintINFO();
            total_mis += map_res.misassemblies;
        }
        if (map_res.wrong_gap_size > 0) {
            INFO ("there are " << map_res.wrong_gap_size << " wrong gaps in path. ");
            path.PrintDEBUG();
            gap_mis += map_res.wrong_gap_size;
        }
        genome_checker.CheckPathEnd(path);
        genome_checker.CheckPathEnd(*path.GetConjPath());
    }
    INFO ("In total found " << total_mis << " misassemblies " << " and " << gap_mis << " gaps.");
}

void PathExtendLauncher::CheckCoverageUniformity() {
    TIME_TRACE_SCOPE("CheckCoverageUniformity");
    if (params_.mode != config::pipeline_type::base)
        return;

    CoverageUniformityAnalyzer coverage_analyzer(graph_, std::min(size_t(1000), stats::Nx(graph_, 50) - 1));
    double median_coverage = coverage_analyzer.CountMedianCoverage();
    double uniformity_fraction = coverage_analyzer.UniformityFraction(unique_data_.unique_variation_, median_coverage);
    if (math::ge(uniformity_fraction, 0.8) and math::ge(median_coverage, 50.0)) {
        WARN("Your data seems to have high uniform coverage depth. It is strongly recommended to use --isolate option.");
    }
}

void PathExtendLauncher::EstimateUniqueEdgesParams() {
    TIME_TRACE_SCOPE("EstimateUniqueEdgesParams");
    bool uniform_coverage = false;
    if (params_.pset.uniqueness_analyser.enabled) {
        INFO("Autodetecting unique edge set parameters...");
        unique_data_.min_unique_length_ = std::max(unique_data_.min_unique_length_, support_.FindMaxMPIS());
        INFO("Minimal unique edge length set to the smallest MP library IS: " << unique_data_.min_unique_length_);

        CoverageUniformityAnalyzer coverage_analyzer(graph_, unique_data_.min_unique_length_);
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
    TIME_TRACE_SCOPE("FillUniqueEdgeStorage");

    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, unique_data_.min_unique_length_, unique_data_.unique_variation_);
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_data_.main_unique_storage_);
}

void PathExtendLauncher::DebugOutputPaths(const PathContainer &paths, const std::string &name) const {
    if (!params_.pe_cfg.debug_output)
        return;

    PathVisualizer visualizer;

    writer_.OutputPaths(paths, params_.etc_dir / (name + ".fasta"));
    if (params_.pe_cfg.output.write_paths) {
        std::ofstream oss(params_.etc_dir / (name + ".dat"));
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            iter.get().Print(oss);
        }
        oss.close();
    }

    if (params_.pe_cfg.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp_, params_.etc_dir / (name + ".dot"), name, paths);
    }
}

void FilterInterstandBulges(PathContainer &paths) {
    TIME_TRACE_SCOPE("FilterInterstandBulges");

    DEBUG ("Try to delete paths with interstand bulges");
    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        if (EndsWithInterstrandBulge(iter.get()))
            iter.get().PopBack();

        if (EndsWithInterstrandBulge(iter.getConjugate())) {
            iter.getConjugate().PopBack();
        }
    }
    DEBUG("deleted paths with interstand bulges");
}

void PathExtendLauncher::RemoveOverlapsAndArtifacts(PathContainer &paths,
                                                    GraphCoverageMap &cover_map,
                                                    const PathExtendResolver &resolver) const {
    TIME_TRACE_SCOPE("PathExtend::RemoveOverlapsAndArtifacts");

    INFO("Finalizing paths");

    INFO("Deduplicating paths");
    Deduplicate(graph_, paths, cover_map, params_.min_edge_len,
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
    TIME_TRACE_SCOPE("PathExtend::CleanPaths");

    if (path_filtration.enabled) {
        paths.FilterPaths(LengthPathCondition(GetLengthCutoff(path_filtration.min_length, path_filtration.rel_cutoff)));
        paths.FilterPaths(func::And(CoveragePathCondition(graph_, path_filtration.min_coverage),
                                    LengthPathCondition(GetLengthCutoff(path_filtration.min_length_for_low_covered, path_filtration.rel_low_covered_cutoff))));
        paths.FilterPaths(func::And(IsolatedPathCondition(graph_),
                                    func::And(LengthPathCondition(GetLengthCutoff(path_filtration.isolated_min_length, path_filtration.rel_isolated_cutoff)),
                                              CoveragePathCondition(graph_, path_filtration.isolated_min_cov))));
    }

    paths.SortByLength();
}

size_t PathExtendLauncher::GetLengthCutoff(size_t abs_cutoff, double rel_cutoff) const {
    int rel_len = int(rel_cutoff * double(cfg::get().ds.RL)) - int(cfg::get().K);
    int abs_len = int(abs_cutoff) - int(cfg::get().K);
    size_t result = (size_t) std::max(0, std::max(rel_len, abs_len));

    DEBUG("Read length relative cutoff " << rel_cutoff << " converted to " << rel_len);
    DEBUG("Read length absolute cutoff " << abs_cutoff << " bp converted to " << result);
    DEBUG("Length cutoff: " << result);
    return result;
}

void PathExtendLauncher::TraverseLoops(PathContainer &paths, GraphCoverageMap &cover_map) const {
    TIME_TRACE_SCOPE("PathExtend::TraverseLoops");

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
    size_t lower_bound = std::max(pset.scaffolding2015.unique_length_lower_bound, pset.scaffolding2015.unique_length_step);

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

Extenders PathExtendLauncher::ConstructReadCloudExtender(const ExtendersGenerator& generator) {
    return generator.MakeReadCloudExtenders(unique_data_.unique_read_cloud_storage_);
}

void PathExtendLauncher::FillPathContainer(size_t lib_index, size_t size_threshold) {
    TIME_TRACE_SCOPE("PathExtend::FillPathContainer");

    INFO("filling path container");
    if (dataset_info_.reads[lib_index].type() == io::LibraryType::TrustedContigs) {
        auto& trusted_paths = gp_.get_mutable<path_extend::TrustedPathsContainer>()[lib_index];
        for (auto & path : trusted_paths)
            unique_data_.long_reads_paths_[lib_index].Create(graph_, std::move(path));

        DebugOutputPaths(unique_data_.long_reads_paths_[lib_index], "trusted_contigs");
        trusted_paths.clear();
    } else {
        std::vector<PathInfo<Graph>> paths;
        gp_.get_mutable<LongReadContainer<Graph>>()[lib_index].SaveAllPaths(paths);
        for (const auto &path: paths) {
            const auto &edges = path.path();
            if (edges.size() <= size_threshold)
                continue;

            auto new_paths = unique_data_.long_reads_paths_[lib_index].CreatePair(graph_, std::move(edges));
            new_paths.first.SetWeight((float) path.weight());
            new_paths.second.SetWeight((float) path.weight());
        }
    }
    DEBUG("Long reads paths " << unique_data_.long_reads_paths_[lib_index].size());
    unique_data_.long_reads_cov_map_[lib_index].AddPaths(unique_data_.long_reads_paths_[lib_index]);
}


void PathExtendLauncher::FillLongReadsCoverageMaps() {
    TIME_TRACE_SCOPE("PathExtend::FillLongReadsCoverageMaps");

    DEBUG("long reads start ")
    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
        DEBUG("lib_index" << lib_index);
        unique_data_.long_reads_paths_.push_back(PathContainer());
        unique_data_.long_reads_cov_map_.push_back(GraphCoverageMap(graph_));
        if (support_.IsForSingleReadExtender(dataset_info_.reads[lib_index])) {
            FillPathContainer(lib_index);
        }
    }
}

void PathExtendLauncher::FillPBUniqueEdgeStorages() {
    TIME_TRACE_SCOPE("PathExtend::FillPBUniqueEdgeStorages");

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
            auto lib = cfg::get().ds.reads[lib_index];
            if (lib.is_paired() and not lib.is_mate_pair()) {
                unique_edge_analyzer_pb.ClearLongEdgesWithPairedLib(lib_index, unique_data_.unique_pb_storage_);
            }
        }

    } else {
        INFO(" with coverage")
        unique_edge_analyzer_pb.FillUniqueEdgeStorage(unique_data_.unique_pb_storage_);
    }
    INFO(unique_data_.unique_pb_storage_.size() << " unique edges");
}

void PathExtendLauncher::FillReadCloudUniqueEdgeStorage() {
    const size_t unique_edge_length = cfg::get().ts_res.edge_length_threshold;
    INFO("Unique edge length: " << unique_edge_length);
    ScaffoldingUniqueEdgeAnalyzer read_cloud_unique_edge_analyzer(gp_, unique_edge_length, unique_data_.unique_variation_);
    read_cloud_unique_edge_analyzer.FillUniqueEdgeStorage(unique_data_.unique_read_cloud_storage_);
}

Extenders PathExtendLauncher::ConstructPBExtenders(const ExtendersGenerator &generator) {
    FillPBUniqueEdgeStorages();
    return generator.MakePBScaffoldingExtenders();
}


Extenders PathExtendLauncher::ConstructExtenders(const GraphCoverageMap &cover_map,
                                                 UsedUniqueStorage &used_unique_storage) {
    TIME_TRACE_SCOPE("PathExtend::ConstructExtenders");

    INFO("Creating main extenders, unique edge length = " << unique_data_.min_unique_length_);
    if (!config::PipelineHelper::IsPlasmidPipeline(params_.mode) &&  (support_.SingleReadsMapped() || support_.HasLongReads()))
        FillLongReadsCoverageMaps();
    ExtendersGenerator generator(dataset_info_, params_, gp_, cover_map,
                                 unique_data_, used_unique_storage, support_);
    Extenders extenders = generator.MakeBasicExtenders();
    DEBUG("Total number of basic extenders is " << extenders.size());

    //long reads scaffolding extenders.

    if (!config::PipelineHelper::IsPlasmidPipeline(params_.mode) && support_.HasLongReads()) {
        if (params_.pset.sm == scaffolding_mode::sm_old) {
            INFO("Will not use new long read scaffolding algorithm in this mode");
        } else {
            utils::push_back_all(extenders, ConstructPBExtenders(generator));
        }
    }

    if (support_.HasReadClouds() and cfg::get().ts_res.read_cloud_resolution_on) {
        if (params_.pset.sm == sm_old) {
            INFO("Will not use read cloud path extend in this mode");
        } else {
            //Creating read cloud unique storage
            FillReadCloudUniqueEdgeStorage();
            INFO("Creating read cloud extenders");
            utils::push_back_all(extenders, ConstructReadCloudExtender(generator));
        }
    }

    if (support_.HasMPReads()) {
        if (params_.pset.sm == scaffolding_mode::sm_old) {
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
                                     const GraphCoverageMap& cover_map) const {
    TIME_TRACE_SCOPE("PathExtend::PolishPaths");
    //Fixes distances for paths gaps and tries to fill them in
    INFO("Closing gaps in paths");

    std::vector<std::shared_ptr<PathGapCloser>> gap_closers;

    gap_closers.push_back(std::make_shared<DijkstraGapCloser>(graph_, params_.max_polisher_gap));

    UniqueData unique_data;
    UsedUniqueStorage used_unique_storage(unique_data.main_unique_storage_);
    ExtendersGenerator generator(dataset_info_, params_, gp_, cover_map, unique_data_, used_unique_storage, support_);

    const auto &paired_indices = gp_.get<UnclusteredPairedInfoIndicesT<Graph>>();
    for (size_t i = 0; i < dataset_info_.reads.lib_count(); i++) {
        auto lib = dataset_info_.reads[i];
        if (lib.type() == io::LibraryType::HQMatePairs || lib.type() == io::LibraryType::MatePairs) {
            auto paired_lib = MakeNewLib(graph_, lib, paired_indices[i]);
            gap_closers.push_back(std::make_shared<MatePairGapCloser>(graph_, params_.max_polisher_gap, paired_lib,
                                                                      unique_data_.main_unique_storage_));
        }
        if (lib.type() == io::LibraryType::Clouds10x
                and cfg::get().ts_res.read_cloud_resolution_on
                and params_.pset.sm != sm_old) {
            INFO("Using read cloud path polisher");
            auto barcode_extractor_ptr =
                std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
            path_extend::ScaffolderParamsConstructor params_constructor;
            auto read_cloud_gap_closer_params = params_constructor.ConstructGapCloserParamsFromCfg(true);
            read_cloud_gap_closer_params.raw_score_threshold_ = cfg::get().ts_res.gap_closer_connection_score_threshold;
            read_cloud_gap_closer_params.relative_coverage_threshold_ = cfg::get().ts_res.gap_closer_relative_coverage_threshold;
            read_cloud_gap_closer_params.edge_length_threshold_ = cfg::get().ts_res.gap_closer_connection_length_threshold;
            size_t tail_threshold = cfg::get().ts_res.scaff_con.path_scaffolder_tail_threshold;
            size_t count_threshold = cfg::get().ts_res.scaff_con.path_scaffolder_count_threshold;
            size_t length_threshold = cfg::get().ts_res.scaff_con.min_edge_length_for_barcode_collection;
            //fixme move to configs
            const size_t scan_bound = 150;
            INFO("Tail threshold: " << tail_threshold);
            INFO("Count threshold: " << count_threshold);
            INFO("Length threshold: " << length_threshold);
            INFO("Connection score threshold: " << read_cloud_gap_closer_params.raw_score_threshold_);
            INFO("Connection relative coverage threshold: " << read_cloud_gap_closer_params.relative_coverage_threshold_);
            INFO("Connection length threshold: " << read_cloud_gap_closer_params.edge_length_threshold_);
            auto cloud_chooser_factory =
                std::make_shared<ReadCloudGapExtensionChooserFactory>(gp_.g,
                                                                   unique_data_.main_unique_storage_,
                                                                   barcode_extractor_ptr,
                                                                   tail_threshold, count_threshold,
                                                                   length_threshold,
                                                                   read_cloud_gap_closer_params,
                                                                   scan_bound);
            auto simple_chooser = generator.MakeSimpleExtensionChooser(i);
            auto simple_chooser_factory = std::make_shared<SameChooserFactory>(gp_.g, simple_chooser);
            auto composite_chooser_factory = std::make_shared<CompositeChooserFactory>(gp_.g, simple_chooser_factory,
                                                                                       cloud_chooser_factory);
            auto cloud_extender_factory = std::make_shared<SimpleExtenderFactory>(gp_, cover_map, used_unique_storage,
                                                                                  composite_chooser_factory);
            gap_closers.push_back(make_shared<PathExtenderGapCloser>(gp_.g, params_.max_polisher_gap, cloud_extender_factory));
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

    PathPolisher polisher(gp_.get<Graph>(), gap_closers);
    result = polisher.PolishPaths(paths);
    result.SortByLength();
    INFO("Gap closing completed")
}


void PathExtendLauncher::FilterPaths(PathContainer &contig_paths) {
    TIME_TRACE_SCOPE("PathExtend::FilterPaths");

    auto default_filtration = params_.pset.path_filtration.end();
    for (auto it = params_.pset.path_filtration.begin(); it != params_.pset.path_filtration.end(); ++it) {
        if (!it->second.enabled)
            continue;

        const auto& filtration_name = it->first;
        if (filtration_name == "default") {
            default_filtration = it;
        } else {
            auto file_name = filtration_name + "_filtered_final_paths.fasta";
            INFO("Finalizing paths - " << filtration_name << ", will be saved to " << file_name);
            PathContainer to_clean(contig_paths.begin(), contig_paths.end());
            CleanPaths(to_clean, it->second);
            DebugOutputPaths(to_clean, filtration_name + "_final_paths");
            writer_.OutputPaths(to_clean, params_.output_dir / file_name);
            contig_name_generator_->PrintStats();
        }
    }
    if (default_filtration != params_.pset.path_filtration.end()) {
        INFO("Finalizing main paths");
        CleanPaths(contig_paths, default_filtration->second);
        DebugOutputPaths(contig_paths, "final_paths");
        contig_name_generator_->PrintStats();
    }
}

void PathExtendLauncher::AddFLPaths(PathContainer &paths) const {
    TIME_TRACE_SCOPE("PathExtend::AddFLPaths");

    bool fl_paths_added = false;
    const auto &single_long_reads = gp_.get<LongReadContainer<Graph>>();
    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
        if (dataset_info_.reads[lib_index].is_full_length_rna_lib()) {
            INFO("Extracting FL paths for lib #" << lib_index);
            std::vector <PathInfo<Graph>> raw_paths;
            single_long_reads[lib_index].SaveAllPaths(raw_paths);
            for (const auto& path: raw_paths) {
                const auto& edges = path.path();

                auto new_paths = paths.CreatePair(graph_, edges);
                new_paths.first.SetWeight((float) path.weight());
                new_paths.second.SetWeight((float) path.weight());
            }
            fl_paths_added = true;
            INFO("Total " << paths.size() << " FL paths were extracted for lib #" << lib_index);
        }
    }

    if (fl_paths_added) {
        INFO("Deduplicating FL paths");
        GraphCoverageMap cover_map(graph_, paths);
        Deduplicate(graph_, paths, cover_map, params_.min_edge_len,  params_.max_path_diff);
        paths.SortByLength();
        INFO("FL Paths deduplicated");
    }
}

class PathContainerCoverageSwitcher {
    const Graph& g_;

    const debruijn_graph::SSCoverageStorage& coverage_storage_;

    bool antisense_;

    double CalculateCoverage(const BidirectionalPath& p, bool reverse) const {
        double res = 0.0;
        double len = 0;
        for (auto e : p) {
            res += coverage_storage_.GetCoverage(e, reverse) * double(g_.length(e));
            len += (double) g_.length(e);
        }
        return res / len;
    }

public:
    PathContainerCoverageSwitcher(const Graph& g, const debruijn_graph::SSCoverageStorage& coverage_storage, bool antisense):
        g_(g), coverage_storage_(coverage_storage), antisense_(antisense) {}


    void Apply(PathContainer& paths) const {
        for (size_t i = 0; i < paths.size(); ++i) {
            if (math::ls(CalculateCoverage(paths.Get(i), antisense_),
                         CalculateCoverage(paths.GetConjugate(i), antisense_))) {
                paths.Swap(i);
            }
        }
    }
};


void PathExtendLauncher::SelectStrandSpecificPaths(PathContainer &paths) const {
    if (!params_.ss.ss_enabled)
        return;

    TIME_TRACE_SCOPE("PathExtend::SelectStrandSpecificPaths");
    INFO("Paths will be printed according to strand-specific coverage");
    size_t lib_index = 0;
    while (lib_index < dataset_info_.reads.lib_count() && !dataset_info_.reads[lib_index].is_graph_constructable()) {
        ++lib_index;
    }
    PathContainerCoverageSwitcher switcher(graph_, gp_.get<SSCoverageContainer>()[lib_index], params_.ss.antisense);
    switcher.Apply(paths);
}

void MakeConjugateEdgePairsDump(ConjugateDeBruijnGraph const & graph) {
    std::ofstream out(cfg::get().output_dir / "conjugate_edge_pairs_dump.info");
    if (!out.is_open()) {
        FATAL_ERROR("Cannot open conjugate_edge_pairs_dump.info for writing");
        return;
    }

    for (EdgeId e : graph.canonical_edges())
        out << e << ' ' << graph.conjugate(e) << '\n';
}

void PathExtendLauncher::Launch() {
    TIME_TRACE_SCOPE("PathExtend");

    INFO("ExSPAnder repeat resolving tool started");
    create_directory(params_.output_dir);
    create_directory(params_.etc_dir);

    CheckCoverageUniformity();

    if (!config::PipelineHelper::IsPlasmidPipeline(params_.mode) && support_.NeedsUniqueEdgeStorage()) {
        TIME_TRACE_SCOPE("FillUniqueEdges");
        //Fill the storage to enable unique edge check
        EstimateUniqueEdgesParams();
        FillUniqueEdgeStorage();
    }

    if (params_.pset.scaffold_graph_params.construct)
        MakeAndOutputScaffoldGraph();

    PathContainer fl_paths;
    AddFLPaths(fl_paths);
    if (fl_paths.size() > 0)
        writer_.OutputPaths(fl_paths, params_.output_dir / "fl_transcripts.fasta");

    auto &contig_paths = gp_.get_mutable<PathContainer>("exSPAnder paths");
    PathExtendResolver resolver(graph_);
    {
        auto seeds = resolver.MakeSimpleSeeds();
        seeds.SortByLength();
        DebugOutputPaths(seeds, "init_paths");

        if (params_.pe_cfg.debug_output)
            MakeConjugateEdgePairsDump(graph_);

        GraphCoverageMap cover_map(graph_);
        UsedUniqueStorage used_unique_storage(unique_data_.main_unique_storage_, graph_);
        Extenders extenders = ConstructExtenders(cover_map, used_unique_storage);
        CompositeExtender composite_extender(graph_, cover_map,
                                             used_unique_storage,
                                             extenders);

        auto paths = resolver.ExtendSeeds(seeds, composite_extender);
        seeds.clear();
        DebugOutputPaths(paths, "raw_paths");

        RemoveOverlapsAndArtifacts(paths, cover_map, resolver);
        DebugOutputPaths(paths, "before_path_polishing");

    PathContainer polished_paths;

    //TODO does path polishing correctly work with coverage map
    PolishPaths(paths, polished_paths, cover_map);
    DebugOutputPaths(polished_paths, "polished_paths");

    //TODO use move assignment to original map here
    GraphCoverageMap polished_map(gp_.g, polished_paths, true);
    RemoveOverlapsAndArtifacts(polished_paths, polished_map, resolver);
    DebugOutputPaths(polished_paths, "overlap_removed");

    //todo discuss
    if (cfg::get().ts_res.path_scaffolding_on and params_.pset.sm != sm_old) {
        const size_t small_path_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
        const size_t large_path_length_threshold = cfg::get().ts_res.long_edge_length_upper_bound;
        PathScaffolder path_scaffolder(gp_, unique_data_.main_unique_storage_,
                                       small_path_length_threshold,
                                       large_path_length_threshold);
        path_scaffolder.MergePaths(polished_paths);

        PolishPaths(polished_paths, contig_paths, polished_map);
        DebugOutputPaths(contig_paths, "merged_polished_paths");

        GraphCoverageMap merged_polished_map(graph_, contig_paths, true);
    DebugOutputPaths(contig_paths, "polished_paths");
    TraverseLoops(contig_paths, polished_map);
        DebugOutputPaths(contig_paths, "loop_traveresed");

    RemoveOverlapsAndArtifacts(contig_paths, polished_map, resolver);
    DebugOutputPaths(contig_paths, "overlap_removed");

    RemoveOverlapsAndArtifacts(contig_paths, polished_map, resolver);
    DebugOutputPaths(contig_paths, "overlap_removed");

    AddFLPaths(contig_paths);

    SelectStrandSpecificPaths(contig_paths);

    CountMisassembliesWithReference(contig_paths);
    INFO("ExSPAnder repeat resolving tool finished");
}

}
