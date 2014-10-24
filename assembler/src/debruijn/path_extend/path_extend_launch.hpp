//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * lc_launch.hpp
 *
 *  Created on: Dec 1, 2011
 *      Author: andrey
 */

#ifndef PATH_EXTEND_LAUNCH_HPP_
#define PATH_EXTEND_LAUNCH_HPP_

#include "pe_config_struct.hpp"
#include "pe_resolver.hpp"
#include "path_extender.hpp"
#include "pe_io.hpp"
#include "path_visualizer.hpp"
#include "loop_traverser.hpp"
#include "long_read_storage.hpp"
#include "next_path_searcher.hpp"


namespace path_extend {

using namespace debruijn_graph;
typedef omnigraph::de::PairedInfoIndicesT<Graph> PairedInfoIndicesT;

inline size_t FindMaxOverlapedLen(const vector<PairedInfoLibrary*>& libes) {
    size_t max = 0;
    for (size_t i = 0; i < libes.size(); ++i) {
        max = std::max(libes[i]->GetISMax(), max);
    }
    return max;
}

inline string GetEtcDir(const std::string& output_dir) {
    return output_dir + cfg::get().pe_params.etc_dir + "/";
}

inline void DebugOutputPaths(const ContigWriter& writer, const conj_graph_pack& gp,
                      const std::string& output_dir, PathContainer& paths,
                      const string& name) {
    PathInfoWriter path_writer;
    PathVisualizer visualizer;
    string etcDir = GetEtcDir(output_dir);
    if (!cfg::get().pe_params.debug_output) {
        return;
    }
    if (cfg::get().pe_params.output.write_paths) {
        writer.writePathEdges(paths, etcDir + name + ".dat");
    }
    if (cfg::get().pe_params.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + name + ".dot", name,
                                             paths);
        path_writer.writePaths(paths, etcDir + name + ".data");
    }
}

inline double GetWeightThreshold(const PairedInfoLibrary* lib, const pe_config::ParamSetT& pset) {
    return lib->IsMp() ? pset.mate_pair_options.weight_threshold : pset.extension_options.weight_threshold;
}

inline double GetSingleThreshold(const PairedInfoLibrary* lib, const pe_config::ParamSetT& pset) {
    return lib->IsMp() ? pset.mate_pair_options.single_threshold : pset.extension_options.single_threshold;
}

inline double GetPriorityCoeff(const PairedInfoLibrary* lib, const pe_config::ParamSetT& pset) {
    return lib->IsMp() ? pset.mate_pair_options.priority_coeff : pset.extension_options.priority_coeff;
}

inline string MakeNewName(const std::string& contigs_name, const std::string& subname) {
    return contigs_name.substr(0, contigs_name.rfind(".fasta")) + "_" + subname + ".fasta";
}

inline void OutputBrokenScaffolds(PathContainer& paths, int k,
                           const ContigWriter& writer,
                           const std::string& filename) {
    if (!cfg::get().pe_params.param_set.scaffolder_options.on
            or !cfg::get().use_scaffolder
            or cfg::get().pe_params.obs == obs_none) {
        return;
    }

    int min_gap = cfg::get().pe_params.obs == obs_break_all ? k / 2 : k;

    ScaffoldBreaker breaker(min_gap);
    breaker.Split(paths);
    breaker.container().SortByLength();
    //writer.writePaths(breaker.container(), filename + ".fasta");
    writer.WritePathsToFASTG(breaker.container(), filename + ".fastg", filename + ".fasta");
}

inline void AddPathsToContainer(const conj_graph_pack& gp,
                         const std::vector<PathInfo<Graph> > paths,
                         size_t size_threshold, PathContainer& result) {
    for (size_t i = 0; i < paths.size(); ++i) {
        auto path = paths.at(i);
        vector<EdgeId> edges = path.getPath();
        if (edges.size() <= size_threshold) {
            continue;
        }
        BidirectionalPath* new_path = new BidirectionalPath(gp.g, edges);
        BidirectionalPath* conj_path = new BidirectionalPath(new_path->Conjugate());
        new_path->SetWeight((float) path.getWeight());
        conj_path->SetWeight((float) path.getWeight());
        result.AddPair(new_path, conj_path);
    }
    DEBUG("Long reads paths " << result.size() << " == ");
}

double GetSingleReadsFilteringThreshold(const io::LibraryType& type) {
    if (type == io::LibraryType::PacBioReads || type == io::LibraryType::SangerReads || type == io::LibraryType::NanoporeReads) {
        return cfg::get().pe_params.long_reads.pacbio_reads.filtering;
    } else if (type == io::LibraryType::TrustedContigs || type == io::LibraryType::UntrustedContigs) {
        return cfg::get().pe_params.long_reads.contigs.filtering;
    }
    return cfg::get().pe_params.long_reads.single_reads.filtering;
}

double GetSingleReadsWeightPriorityThreshold(const io::LibraryType& type) {
    if (type == io::LibraryType::PacBioReads || type == io::LibraryType::SangerReads || type == io::LibraryType::NanoporeReads) {
        return cfg::get().pe_params.long_reads.pacbio_reads.weight_priority;
    } else if (type == io::LibraryType::TrustedContigs || type == io::LibraryType::UntrustedContigs) {
        return cfg::get().pe_params.long_reads.contigs.weight_priority;
    }
    return cfg::get().pe_params.long_reads.single_reads.weight_priority;
}

double GetSingleReadsUniqueEdgePriorityThreshold(const io::LibraryType& type) {
    if (cfg::get().ds.single_cell &&
            (type == io::LibraryType::PacBioReads || type == io::LibraryType::SangerReads || type == io::LibraryType::NanoporeReads)) {
        return 10000.0;
    }
    if (type == io::LibraryType::PacBioReads || type == io::LibraryType::SangerReads || type == io::LibraryType::NanoporeReads) {
        return cfg::get().pe_params.long_reads.pacbio_reads.unique_edge_priority;
    } else if (type == io::LibraryType::TrustedContigs || type == io::LibraryType::UntrustedContigs) {
        return cfg::get().pe_params.long_reads.contigs.unique_edge_priority;
    }
    return cfg::get().pe_params.long_reads.single_reads.unique_edge_priority;
}

bool HasOnlyMPLibs() {
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (!((cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs || cfg::get().ds.reads[i].type() == io::LibraryType::HQMatePairs)
                && cfg::get().ds.reads[i].data().mean_insert_size > 0.0)) {
            return false;
        }
    }
    return true;
}

bool UseCoverageResolverForSingleReads(const io::LibraryType& type) {
    return HasOnlyMPLibs() && (type == io::LibraryType::HQMatePairs);

}

inline vector<SimpleExtender*> MakeLongReadsExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, size_t max_loops, const std::string& output_dir) {
    vector<SimpleExtender*> extends;
    for (size_t i = 0; i < gp.single_long_reads.size(); ++i) {
        PathContainer paths;
        AddPathsToContainer(gp, gp.single_long_reads[i].GetAllPaths(), 1, paths);
        if (paths.size() == 0) {
            continue;
        }
        ContigWriter writer(gp.g);
        DebugOutputPaths(writer, gp, output_dir, paths, "long_reads");
        ExtensionChooser * longReadEC = new LongReadsExtensionChooser(gp.g, paths, GetSingleReadsFilteringThreshold(cfg::get().ds.reads[i].type()),
                                                                      GetSingleReadsWeightPriorityThreshold(cfg::get().ds.reads[i].type()),
                                                                      GetSingleReadsUniqueEdgePriorityThreshold(cfg::get().ds.reads[i].type()));
        SimpleExtender * longReadExtender = new SimpleExtender(gp, cov_map, longReadEC, 10000,  //FIXME
                                                                       max_loops, true, UseCoverageResolverForSingleReads(cfg::get().ds.reads[i].type()));
        extends.push_back(longReadExtender);
    }
    return extends;
}

inline vector<SimpleExtender *> MakePEExtenders(const conj_graph_pack& gp, const pe_config::ParamSetT& pset, vector<PairedInfoLibrary*>& libs,
                                                const GraphCoverageMap& cov_map, bool investigate_loops) {
    vector<SimpleExtender *> extends;
    for (size_t i = 0; i < libs.size(); ++i) {
        WeightCounter* wc = new PathCoverWeightCounter(gp.g, libs[i], GetWeightThreshold(libs[i], pset), GetSingleThreshold(libs[i], pset));
        wc->setNormalizeWeight(pset.normalize_weight);
        SimpleExtensionChooser * extension = new SimpleExtensionChooser(gp.g, wc, GetPriorityCoeff(libs[i], pset));
        extends.push_back(new SimpleExtender(gp, cov_map, extension, libs[i]->GetISMax(), pset.loop_removal.max_loops, investigate_loops, false));
    }
    return extends;
}

inline vector<ScaffoldingPathExtender*> MakeScaffoldingExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, const pe_config::ParamSetT& pset,
                                                                vector<PairedInfoLibrary *>& libs) {
    GapJoiner * gapJoiner = new HammingGapJoiner(gp.g, pset.scaffolder_options.min_gap_score,
                                                 (int) (pset.scaffolder_options.max_must_overlap * (double) gp.g.k()),
                                                 (int) (pset.scaffolder_options.max_can_overlap * (double) gp.g.k()), pset.scaffolder_options.short_overlap);
    vector<ScaffoldingPathExtender*> scafPEs;
    for (size_t i = 0; i < libs.size(); ++i) {
        WeightCounter* counter = new ReadCountWeightCounter(gp.g, libs[i]);
        double prior_coef = GetPriorityCoeff(libs[i], pset);
        ScaffoldingExtensionChooser * scaff_chooser = new ScaffoldingExtensionChooser(gp.g, counter, prior_coef);
        scafPEs.push_back(new ScaffoldingPathExtender(gp, cov_map, scaff_chooser, gapJoiner, libs[i]->GetISMax(), pset.loop_removal.max_loops, false));
    }
    return scafPEs;
}
inline size_t CountEdgesInGraph(const Graph& g) {
    size_t count = 0;
    for (auto iter = g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        count++;
    }
    return count;
}

inline size_t GetNumberMPPaths(const Graph& g) {
    size_t count_edge = CountEdgesInGraph(g);
    if (count_edge < 1000) {
        return 1000;
    }
    if (count_edge < 10000) {
        return 100;
    }
    return 50;
}
inline vector<SimpleExtender *> MakeMPExtenders(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, const pe_config::ParamSetT& pset,
                                                const PathContainer& paths, vector<PairedInfoLibrary *>& libs) {
    size_t max_number_of_paths_to_search = GetNumberMPPaths(gp.g);
    DEBUG("max number of mp paths " << max_number_of_paths_to_search);
    vector<SimpleExtender *> mpPEs;
    for (size_t i = 0; i < libs.size(); ++i) {
        MatePairExtensionChooser* chooser = new MatePairExtensionChooser(gp.g, *libs[i], paths, max_number_of_paths_to_search);
        SimpleExtender* mp_extender = new SimpleExtender(gp, cov_map, chooser, libs[i]->GetISMax(), pset.loop_removal.mp_max_loops, true, false);
        mpPEs.push_back(mp_extender);
    }
    return mpPEs;
}

inline string LibStr(size_t count) {
    return count == 1 ? "library" : "libraries";
}

inline void ClonePathContainer(PathContainer& spaths, PathContainer& tpaths, GraphCoverageMap& tmap) {
    tpaths.clear();
    tmap.Clear();

    for (auto iter = spaths.begin(); iter != spaths.end(); ++iter) {
        BidirectionalPath& path = *iter.get();
        BidirectionalPath* new_path = new BidirectionalPath(path.graph());
        new_path->Subscribe(&tmap);
        new_path->PushBack(path);

        BidirectionalPath& cpath = *iter.getConjugate();
        BidirectionalPath* new_cpath = new BidirectionalPath(cpath.graph());
        new_cpath->Subscribe(&tmap);
        new_cpath->PushBack(cpath);

        tpaths.AddPair(new_path, new_cpath);
    }
}

inline void FinalizePaths(PathContainer& paths, GraphCoverageMap& cover_map, size_t max_overlap, bool mate_pairs = false) {
    ContigWriter writer(cover_map.graph());
    PathExtendResolver resolver(cover_map.graph());

    resolver.removeOverlaps(paths, cover_map, max_overlap, cfg::get().pe_params.param_set.remove_overlaps, cfg::get().pe_params.cut_all_overlaps);
    if (mate_pairs) {
        resolver.RemoveMatePairEnds(paths, max_overlap);
    }
    if (cfg::get().avoid_rc_connections) {
        paths.FilterInterstandBulges();
    }
    paths.FilterEmptyPaths();
    if (!mate_pairs) {
        resolver.addUncoveredEdges(paths, cover_map);
    }
    paths.SortByLength();
    paths.ResetPathsId();
}

inline void TraverseLoops(PathContainer& paths, GraphCoverageMap& cover_map, ContigsMaker* extender) {
    INFO("Traversing tandem repeats");
    LoopTraverser loopTraverser(cover_map.graph(), cover_map, extender);
    loopTraverser.TraverseAllLoops();
    paths.SortByLength();
    paths.ResetPathsId();
}


inline void ResolveRepeatsManyLibs(conj_graph_pack& gp,
		vector<PairedInfoLibrary *>& libs,
		vector<PairedInfoLibrary *>& scaff_libs,
		vector<PairedInfoLibrary *>& mp_libs,
		const std::string& output_dir,
		const std::string& contigs_name,
		bool traversLoops,
		boost::optional<std::string> broken_contigs) {

	INFO("Path-Extend repeat resolving tool started");

	make_dir(output_dir);
    make_dir(GetEtcDir(output_dir));
	const pe_config::ParamSetT& pset = cfg::get().pe_params.param_set;

	ContigWriter writer(gp.g);

//make pe + long reads extenders
    GraphCoverageMap cover_map(gp.g);
    vector<SimpleExtender *> usualPEs = MakePEExtenders(gp, pset, libs, cover_map, false);
    vector<SimpleExtender*> long_reads_extenders = MakeLongReadsExtender(gp, cover_map, pset.loop_removal.max_loops, output_dir);
    vector<SimpleExtender *> shortLoopPEs = MakePEExtenders(gp, pset, libs, cover_map, true);
    vector<ScaffoldingPathExtender*> scafPEs = MakeScaffoldingExtender(gp, cover_map, pset, scaff_libs);
    vector<PathExtender *> all_libs(long_reads_extenders.begin(), long_reads_extenders.end());
    all_libs.insert(all_libs.end(), usualPEs.begin(), usualPEs.end());
    all_libs.insert(all_libs.end(), shortLoopPEs.begin(), shortLoopPEs.end());
    all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
    INFO("Using " << libs.size() << " paired-end " << LibStr(libs.size()));
    INFO("Using " << scaff_libs.size() << " paired-end scaffolding " << LibStr(scaff_libs.size()));
    INFO("Using " << mp_libs.size() << " mate-pair " << LibStr(mp_libs.size()));
    INFO("Using " << long_reads_extenders.size() << " single read " << LibStr(long_reads_extenders.size()));
    INFO("Scaffolder is " << (pset.scaffolder_options.on ? "on" : "off"));

    size_t max_over = max(FindMaxOverlapedLen(libs), gp.g.k() + 100);

    CompositeExtender * mainPE = new CompositeExtender(gp.g, cover_map, all_libs, max_over);

//extend pe + long reads
    PathExtendResolver resolver(gp.g);
    auto seeds = resolver.makeSimpleSeeds();
	DebugOutputPaths(writer, gp, output_dir, seeds, "init_paths");
    seeds.SortByLength();
	seeds.ResetPathsId();
	INFO("Growing paths using paired-end and long single reads");
	auto paths = resolver.extendSeeds(seeds, *mainPE);
	paths.SortByLength();
	//paths.ResetPathsId();
	DebugOutputPaths(writer, gp, output_dir, paths, "pe_overlaped_paths");

    PathContainer clone_paths;
    GraphCoverageMap clone_map(gp.g);
    bool mp_exist = mp_libs.size() > 0;
    if (mp_exist) {
        ClonePathContainer(paths, clone_paths, clone_map);
    }

    FinalizePaths(paths, cover_map, max_over);
    if (broken_contigs.is_initialized()) {
        OutputBrokenScaffolds(paths, (int) gp.g.k(), writer,
                              output_dir + (mp_exist ? "pe_contigs" : broken_contigs.get()));
    }
    writer.WritePathsToFASTG(paths, GetEtcDir(output_dir) + "pe_before_traversal.fastg", GetEtcDir(output_dir) + "pe_before_traversal.fasta");
    DebugOutputPaths(writer, gp, output_dir, paths, "before_traverse_pe");
    if (traversLoops) {
        TraverseLoops(paths, cover_map, mainPE);
    }
    DebugOutputPaths(writer, gp, output_dir, paths, (mp_exist ? "final_pe_paths" : "final_paths"));
    writer.WritePathsToFASTG(paths,
                             output_dir + (mp_exist ? "pe_scaffolds" : contigs_name) + ".fastg",
                             output_dir + (mp_exist ? "pe_scaffolds" : contigs_name) + ".fasta" );

    cover_map.Clear();
    paths.DeleteAllPaths();
    if (!mp_exist) {
        return;
    }

//MP
    INFO("Adding mate-pairs");
    DebugOutputPaths(writer, gp, output_dir, clone_paths, "before_mp_paths");
    vector<SimpleExtender*> mpPEs = MakeMPExtenders(gp, clone_map, pset, clone_paths, mp_libs);
    max_over = std::max(FindMaxOverlapedLen(mp_libs), FindMaxOverlapedLen(libs));
	all_libs.clear();
	all_libs.insert(all_libs.end(), long_reads_extenders.begin(),
			long_reads_extenders.end());
	all_libs.insert(all_libs.end(), shortLoopPEs.begin(), shortLoopPEs.end());
	all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
	all_libs.insert(all_libs.end(), mpPEs.begin(), mpPEs.end());
	CompositeExtender* mp_main_pe = new CompositeExtender(gp.g, clone_map, all_libs, max_over);

	INFO("Growing paths using mate-pairs");
	auto mp_paths = resolver.extendSeeds(clone_paths, *mp_main_pe);
	FinalizePaths(mp_paths, clone_map, max_over, true);
	DebugOutputPaths(writer, gp, output_dir, mp_paths, "mp_final_paths");
	writer.WritePathsToFASTG(mp_paths, GetEtcDir(output_dir) + "mp_prefinal.fastg", GetEtcDir(output_dir) + "mp_prefinal.fasta");

    DEBUG("Paths are grown with mate-pairs");
    if (cfg::get().pe_params.debug_output) {
        writer.writePaths(mp_paths, output_dir + "mp_paths.fasta");
    }
//MP end

//pe again
    all_libs.clear();
    all_libs.insert(all_libs.end(), long_reads_extenders.begin(),
                    long_reads_extenders.end());
    all_libs.insert(all_libs.end(), shortLoopPEs.begin(), shortLoopPEs.end());
    all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
    max_over = FindMaxOverlapedLen(libs);
    CompositeExtender* last_extender = new CompositeExtender(gp.g, clone_map, all_libs, max_over);
    auto last_paths = resolver.extendSeeds(mp_paths, *last_extender);
    FinalizePaths(last_paths, clone_map, max_over);

    writer.WritePathsToFASTG(last_paths, GetEtcDir(output_dir) + "mp_before_traversal.fastg", GetEtcDir(output_dir) + "mp_before_traversal.fasta");
    DebugOutputPaths(writer, gp, output_dir, last_paths, "before_traverse_mp");
    TraverseLoops(last_paths, clone_map, last_extender);

//result
    if (broken_contigs.is_initialized()) {
        OutputBrokenScaffolds(last_paths, (int) gp.g.k(), writer, output_dir + broken_contigs.get());
    }
    DebugOutputPaths(writer, gp, output_dir, last_paths, "last_paths");
    writer.WritePathsToFASTG(last_paths, output_dir + contigs_name + ".fastg", output_dir + contigs_name + ".fasta");

    last_paths.DeleteAllPaths();
    seeds.DeleteAllPaths();
    mp_paths.DeleteAllPaths();
    clone_paths.DeleteAllPaths();

    INFO("Path-Extend repeat resolving tool finished");
}

template<class Index>
inline PairedInfoLibrary* MakeNewLib(conj_graph_pack::graph_t& g,
                                     const Index& paired_index,
                                     size_t index) {
    size_t read_length = cfg::get().ds.reads[index].data().read_length;
    size_t is = (size_t) cfg::get().ds.reads[index].data().mean_insert_size;
    int is_min = (int) cfg::get().ds.reads[index].data().insert_size_left_quantile;
    int is_max = (int) cfg::get().ds.reads[index].data().insert_size_right_quantile;
    int var = (int) cfg::get().ds.reads[index].data().insert_size_deviation;
    bool is_mp = cfg::get().ds.reads[index].type() == io::LibraryType::MatePairs;
    PairedInfoLibrary* lib = new PairedInfoLibraryWithIndex<decltype(paired_index[index])>(cfg::get().K, g, read_length,
                                                                   is, is_min > 0.0 ? size_t(is_min) : 0, is_max > 0.0 ? size_t(is_max) : 0,
                                                                   size_t(var),
                                                                   paired_index[index], is_mp,
                                                                   cfg::get().ds.reads[index].data().insert_size_distribution);
    return lib;
}

inline void DeleteLibs(vector<PairedInfoLibrary *>& libs) {
    for (size_t j = 0; j < libs.size(); ++j)
        delete libs[j];
}

inline bool InsertSizeCompare(const PairedInfoLibrary* lib1,
                              const PairedInfoLibrary* lib2) {
    return lib1->GetISMax() < lib2->GetISMax();
}

inline void ResolveRepeatsPe(conj_graph_pack& gp, const std::string& output_dir, const std::string& contigs_name, bool traverseLoops,
                             boost::optional<std::string> broken_contigs, bool use_auto_threshold = true) {
    vector<PairedInfoLibrary*> rr_libs;
    vector<PairedInfoLibrary*> mp_libs;
    vector<PairedInfoLibrary*> scaff_libs;
    for (size_t i = 0; i < gp.clustered_indices.size(); ++i) {
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd && cfg::get().ds.reads[i].data().mean_insert_size > 0.0) {
            PairedInfoLibrary* lib = MakeNewLib(gp.g, gp.clustered_indices, i);
            if (use_auto_threshold) {
                lib->SetSingleThreshold(cfg::get().ds.reads[i].data().pi_threshold);
                INFO("Threshold for library #" << i << " is " << cfg::get().ds.reads[i].data().pi_threshold);
            }
            rr_libs.push_back(lib);
        }
    }
    for (size_t i = 0; i < gp.paired_indices.size(); ++i) {
        if ((cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs || cfg::get().ds.reads[i].type() == io::LibraryType::HQMatePairs)
                && cfg::get().ds.reads[i].data().mean_insert_size > 0.0) {
            PairedInfoLibrary* lib = MakeNewLib(gp.g, gp.paired_indices, i);
            mp_libs.push_back(lib);
        }
    }
    std::sort(rr_libs.begin(), rr_libs.end(), InsertSizeCompare);
    if (cfg::get().use_scaffolder && cfg::get().pe_params.param_set.scaffolder_options.on) {
        for (size_t i = 0; i < gp.scaffolding_indices.size(); ++i) {
            if ((cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd) && cfg::get().ds.reads[i].data().mean_insert_size > 0.0) {
                PairedInfoLibrary* lib = MakeNewLib(gp.g, gp.scaffolding_indices, i);
                scaff_libs.push_back(lib);
            }
        }
    }
    std::sort(scaff_libs.begin(), scaff_libs.end(), InsertSizeCompare);
    ResolveRepeatsManyLibs(gp, rr_libs, scaff_libs, mp_libs, output_dir, contigs_name, traverseLoops, broken_contigs);
    DeleteLibs(rr_libs);
    DeleteLibs(mp_libs);
    DeleteLibs(scaff_libs);
}


} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
