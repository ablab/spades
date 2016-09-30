//
// Created by andrey on 14.11.16.
//

#include "launcher.hpp"

namespace path_extend {

using namespace debruijn_graph;
using namespace std;


vector<shared_ptr<ConnectionCondition>>
    PathExtendLauncher::ConstructPairedConnectionConditions(const ScaffoldingUniqueEdgeStorage& edge_storage) const {

    vector<shared_ptr<ConnectionCondition>> conditions;
    const pe_config::ParamSetT::ScaffoldGraphParamsT &params = params_.pset.scaffold_graph_params;

    //extract method!
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
    LengthEdgeCondition edge_condition(gp_.g, edge_storage.GetMinLength());

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

//I would suggest to move it out of launcher
//ANSWER: now about no =)
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

//what does this method do? Does it even change the state?
//ANSWER: Construcnts and prints scaffold graph. I think it's a part of pipeline
void PathExtendLauncher::ScaffoldGraphConstruction() const {
    //Scaffold graph
    shared_ptr<scaffold_graph::ScaffoldGraph> scaffold_graph;
    if (params_.pset.scaffold_graph_params.construct) {
        debruijn_graph::GenomeConsistenceChecker genome_checker(gp_, main_unique_storage_,
                                                                params_.pset.genome_consistency_checker.max_gap,
                                                                params_.pset.genome_consistency_checker.relative_max_gap);
        scaffold_graph = ConstructScaffoldGraph(main_unique_storage_);
        if (params_.pset.scaffold_graph_params.output) {
            PrintScaffoldGraph(*scaffold_graph,
                               main_unique_storage_.GetSet(),
                               genome_checker,
                               params_.etc_dir + "scaffold_graph");
        }
    }
}

//Out of launcher! Suggest creating special class
//ANSWER:  I think no, it's part of pipeline, GenomeConsistenceChecker is a class and this is just printing things
void PathExtendLauncher::CountMisassembliesWithReference(const PathContainer &paths) const {
    if (gp_.genome.size() == 0)
        return;

    debruijn_graph::GenomeConsistenceChecker genome_checker(gp_, main_unique_storage_,
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

void PathExtendLauncher::FillUniqueEdgeStorage() {
    //do we still have this mode?
    //ANSWER: so far yes
    //Since new pipeline is experimental -- use it for multicell with long reads or mate-pairs only
    if (params_.pset.sm == sm_old ||
       (params_.pset.sm == sm_old_pe_2015 && !HasLongReads(dataset_info_) && !HasMPReads(dataset_info_)))
        return;

    bool uniform_coverage = false;
    if (params_.pset.uniqueness_analyser.enabled) {
        INFO("Autodetecting unique edge set parameters...");
        size_t max_is = min_unique_length_;
        //extract function
        for (size_t i = 0; i < dataset_info_.reads.lib_count(); ++i) {
            //only for mate pairs?
            //ANSWER: yes, if no MPs use default since PE insert size is usually too small for unique edge length threshold
            if (dataset_info_.reads[i].is_mate_pair()) {
                max_is = max(max_is, (size_t) dataset_info_.reads[i].data().mean_insert_size);
            }
        }
        //unpleasant side effect (not predictable from function name), cant we set it earlier in a separate method?
        //ANSWER: agreed, noticed long time ago -- DimaA heritage, will talk
        min_unique_length_ = max_is;
        INFO("Minimal unique edge length set to the smallest MP library IS: " << min_unique_length_);

        CoverageUniformityAnalyzer coverage_analyzer(gp_.g, min_unique_length_);
        double median_coverage = coverage_analyzer.CountMedianCoverage();
        double uniformity_fraction = coverage_analyzer.UniformityFraction(unique_variation_, median_coverage);
        INFO ("median coverage for edges longer than " << min_unique_length_ << " is " << median_coverage <<
            " uniformity " << size_t(uniformity_fraction * 100) << "%");
        if (math::gr(uniformity_fraction, params_.pset.uniqueness_analyser.uniformity_fraction_threshold)) {
            uniform_coverage = true;
        }
        if (!uniform_coverage) {
            unique_variation_ = params_.pset.uniqueness_analyser.nonuniform_coverage_variation;
            INFO("Coverage is not uniform, we do not rely on coverage for long edge uniqueness");
        }

    } else {
        INFO("Unique edge set constructed with parameters from config : length " << min_unique_length_
                 << " variation " << unique_variation_);
    }
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, min_unique_length_, unique_variation_);
    unique_edge_analyzer.FillUniqueEdgeStorage(main_unique_storage_);
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

void PathExtendLauncher::OutputBrokenScaffolds(const PathContainer &paths, const std::string &filename) const {
    if (!params_.pset.scaffolder_options.enabled ||
        !params_.use_scaffolder || params_.pe_cfg.obs == obs_none) {
        return;
    }

    int min_gap = int(params_.pe_cfg.obs == obs_break_all ? gp_.g.k() / 2 : gp_.g.k());

    ScaffoldBreaker breaker(min_gap, paths);
    breaker.container().SortByLength();
    writer_.OutputPaths(breaker.container(), filename);
}

void PathExtendLauncher::FinalizePaths(PathContainer &paths,
                                       GraphCoverageMap &cover_map,
                                       const PathExtendResolver&resolver) const {

    if (params_.pset.remove_overlaps) {
        resolver.RemoveOverlaps(paths, cover_map, params_.min_edge_len, params_.max_path_diff,
                                 params_.pset.cut_all_overlaps,
                                 (params_.mode == config::pipeline_type::moleculo));
    } else if (params_.mode == config::pipeline_type::rna) {
        //is it really different, or just decided not to bother with configuration? :)
        //ANSWER: Well, it is, Can make a a flag and add some ifs inside, but looks like we better not
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

void PathExtendLauncher::TraverseLoops(PathContainer &paths,
                                         GraphCoverageMap &cover_map) const {
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

vector<shared_ptr<PathExtender>>  PathExtendLauncher::ConstructExtenders(const GraphCoverageMap& cover_map) {
    const pe_config::ParamSetT &pset = params_.pset;

    INFO("Creating main exteders, unique edge length = " << min_unique_length_);
    ExtendersGenerator generator(dataset_info_, params_, gp_, cover_map);
    vector<shared_ptr<PathExtender>> extenders = generator.MakeBasicExtenders(main_unique_storage_);

    //long reads scaffolding extenders.
    if (HasLongReads(dataset_info_)) {
        //extract method OK
        vector<PathContainer> long_reads_paths(dataset_info_.reads.lib_count());
        ScaffoldingUniqueEdgeStorage unique_storage_pb;
        vector<shared_ptr<GraphCoverageMap>> long_reads_cov_map;
        push_back_all(extenders,
                      generator.MakePBScaffoldingExtenders(unique_storage_pb,
                                                          long_reads_paths,
                                                          long_reads_cov_map));
    }

    if (HasMPReads(dataset_info_)) {
        //extract method OK
        push_back_all(extenders,
                      generator.MakeMPExtenders(main_unique_storage_));

        size_t cur_length = min_unique_length_ - pset.scaffolding2015.unique_length_step;
        size_t lower_bound = max(pset.scaffolding2015.unique_length_lower_bound, pset.scaffolding2015.unique_length_step);

        while (cur_length > lower_bound) {
            INFO("Adding extender with length " << cur_length);
            ScaffoldingUniqueEdgeAnalyzer additional_edge_analyzer(gp_, (size_t) cur_length, unique_variation_);
            unique_storages_.push_back(new ScaffoldingUniqueEdgeStorage());
            additional_edge_analyzer.FillUniqueEdgeStorage(*unique_storages_.back());

            auto additional_extenders = generator.MakeMPExtenders(*unique_storages_.back());
            extenders.insert(extenders.end(), additional_extenders.begin(), additional_extenders.end());
            cur_length -= pset.scaffolding2015.unique_length_step;
        }
        if (min_unique_length_ > lower_bound) {
            INFO("Adding final extender with length " << lower_bound);
            ScaffoldingUniqueEdgeAnalyzer additional_edge_analyzer(gp_, lower_bound, unique_variation_);
            unique_storages_.push_back(new ScaffoldingUniqueEdgeStorage());
            additional_edge_analyzer.FillUniqueEdgeStorage(*unique_storages_.back());
            auto additional_extenders = generator.MakeMPExtenders(*unique_storages_.back());
            extenders.insert(extenders.end(), additional_extenders.begin(), additional_extenders.end());
        }
    }
    INFO("Total number of extenders is " << extenders.size());
    return extenders;
}

//PrepareDirs?
//inline in launch!
//ANSWER: Where?
void PathExtendLauncher::PrepareForRun() const {
    make_dir(params_.output_dir);
    make_dir(params_.etc_dir);
}

void PathExtendLauncher::PolishPaths(const PathContainer &paths, PathContainer &result) const {
    //Fixes distances for paths gaps and tries to fill them in
    INFO("Closing gaps in paths");
    PathPolisher polisher(gp_, dataset_info_, main_unique_storage_, params_.max_polisher_gap);
    polisher.PolishPaths(paths, result);
    result.SortByLength();
    INFO("Gap closing completed")
}

void PathExtendLauncher::Launch() {
    INFO("ExSPAnder repeat resolving tool started");
    PrepareForRun();

    //Fill the storage to enable unique edge check
    FillUniqueEdgeStorage();

    ScaffoldGraphConstruction();

    PathExtendResolver resolver(gp_.g);

    auto seeds = resolver.MakeSimpleSeeds();
    seeds.SortByLength();
    DebugOutputPaths(seeds, "init_paths");

    GraphCoverageMap cover_map(gp_.g);
    vector<shared_ptr<PathExtender>> extenders = ConstructExtenders(cover_map);
    shared_ptr<CompositeExtender> composite_extender = make_shared<CompositeExtender>(gp_.g, cover_map, extenders,
                                                                                      main_unique_storage_,
                                                                                      params_.max_path_diff,
                                                                                      params_.pset.extension_options.max_repeat_length,
                                                                                      params_.detect_repeats_online);

    auto paths = resolver.ExtendSeeds(seeds, *composite_extender);
    paths.FilterEmptyPaths();
    paths.SortByLength();
    DebugOutputPaths(paths, "raw_paths");

    TraverseLoops(paths, cover_map);
    DebugOutputPaths(paths, "loop_traveresed");

    PathContainer polished_paths;
    PolishPaths(paths, polished_paths);

    DebugOutputPaths(polished_paths, "polished_paths");
    GraphCoverageMap polished_map(gp_.g, polished_paths, true);
    FinalizePaths(polished_paths, polished_map, resolver);

    if (params_.output_broken_scaffolds) {
        OutputBrokenScaffolds(polished_paths, params_.output_dir + params_.broken_contigs);
    }

    DebugOutputPaths(polished_paths, "final_paths");
    writer_.OutputPaths(polished_paths, params_.output_dir + params_.contigs_name);

    CountMisassembliesWithReference(polished_paths);

    INFO("ExSPAnder repeat resolving tool finished");
}


}
