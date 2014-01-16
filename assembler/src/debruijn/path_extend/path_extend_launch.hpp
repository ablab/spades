//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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
        writer.writePathEdges(paths, output_dir + name + ".dat");
    }
    if (cfg::get().pe_params.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + name + ".dot", name,
                                             paths);
        path_writer.writePaths(paths, etcDir + name + ".data");
    }
}

inline double GetWeightThreshold(const PairedInfoLibrary* lib,
                          const pe_config::ParamSetT& pset) {
    return lib->IsMp() ?
            pset.mate_pair_options.weight_threshold :
            pset.extension_options.weight_threshold;
}

inline double GetSingleThreshold(const PairedInfoLibrary* lib,
                          const pe_config::ParamSetT& pset) {
    return lib->IsMp() ?
            pset.mate_pair_options.single_threshold :
            pset.extension_options.single_threshold;
}

inline double GetPriorityCoeff(const PairedInfoLibrary* lib,
                        const pe_config::ParamSetT& pset) {
    return lib->IsMp() ?
            pset.mate_pair_options.priority_coeff :
            pset.extension_options.priority_coeff;
}

inline string MakeNewName(const std::string& contigs_name,
		const std::string& subname){
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
                         const std::vector<PathInfo<Graph> >& paths,
                         size_t size_threshold, PathContainer& result) {
    for (size_t i = 0; i < paths.size(); ++i) {
        PathInfo<Graph> path = paths[i];
        if (path.getPath().size() <= size_threshold) {
            continue;
        }
        vector<EdgeId> edges = path.getPath();
        BidirectionalPath* new_path = new BidirectionalPath(gp.g, edges);
        BidirectionalPath* conj_path = new BidirectionalPath(new_path->Conjugate());
        new_path->SetWeight((double) path.getWeight());
        conj_path->SetWeight((double) path.getWeight());
        result.AddPair(new_path, conj_path);
    }
    DEBUG("Long reads paths " << result.size() << " == ");
}

inline vector<SimpleExtender *> MakePEExtenders(const conj_graph_pack& gp,
                                         const pe_config::ParamSetT& pset,
                                         vector<PairedInfoLibrary*>& libs,
                                         bool investigate_loops) {
    vector<SimpleExtender *> extends;
    for (size_t i = 0; i < libs.size(); ++i) {
        WeightCounter* wc = new PathCoverWeightCounter(
                gp.g, libs[i], GetWeightThreshold(libs[i], pset),
                GetSingleThreshold(libs[i], pset));
        wc->setNormalizeWeight(pset.normalize_weight);
        SimpleExtensionChooser * extension = new SimpleExtensionChooser(
                gp.g, wc, GetPriorityCoeff(libs[i], pset));
        extends.push_back(new SimpleExtender(gp.g, extension, libs[i]->GetISMax(), pset.loop_removal.max_loops, investigate_loops));
    }
    return extends;
}

inline vector<SimpleExtender*> MakeLongReadsExtender(
        const conj_graph_pack& gp, const vector<PathStorageInfo<Graph> >& reads,
        size_t max_loops, const std::string& output_dir) {
    vector<SimpleExtender*> extends;
    for (size_t i = 0; i < reads.size(); ++i) {
        if (reads[i].GetPaths().size() == 0) {
            continue;
        }
        PathContainer paths;
        AddPathsToContainer(gp, reads[i].GetPaths(), 1, paths);
        ContigWriter writer(gp.g);
        DebugOutputPaths(writer, gp, output_dir, paths, "long_reads");
        ExtensionChooser * longReadEC = new LongReadsExtensionChooser(
                gp.g, paths, reads[i].GetFilteringThreshold(),
                reads[i].GetWeightPriorityThreshold(),
                reads[i].GetUniqueEdgePriorityThreshold());
        SimpleExtender * longReadExtender = new SimpleExtender(gp.g,
                                                               longReadEC,
                                                               10000, //FIXME
                                                               max_loops,
                                                               true);
        extends.push_back(longReadExtender);
    }
    return extends;
}
inline vector<ScaffoldingPathExtender*> MakeScaffoldingExtender(
        const conj_graph_pack& gp, const pe_config::ParamSetT& pset,
        vector<PairedInfoLibrary *>& libs) {
    GapJoiner * gapJoiner =
            new HammingGapJoiner(
                    gp.g,
                    pset.scaffolder_options.min_gap_score,
                    (int) (pset.scaffolder_options.max_must_overlap * (double) gp.g.k()),
                    (int) (pset.scaffolder_options.max_can_overlap * (double) gp.g.k()),
                    pset.scaffolder_options.short_overlap);
    vector<ScaffoldingPathExtender*> scafPEs;
    for (size_t i = 0; i < libs.size(); ++i) {
        WeightCounter* counter = new ReadCountWeightCounter(gp.g, libs[i]);
        double prior_coef = GetPriorityCoeff(libs[i], pset);
        ScaffoldingExtensionChooser * scaff_chooser =
                new ScaffoldingExtensionChooser(gp.g, counter, prior_coef);
        scafPEs.push_back(
                new ScaffoldingPathExtender(gp.g,
                                            scaff_chooser, gapJoiner,
                                            libs[i]->GetISMax(), pset.loop_removal.max_loops, false));
    }
    return scafPEs;
}

inline vector<SimpleExtender *> MakeMPExtenders(const conj_graph_pack& gp,
                                         const pe_config::ParamSetT& pset,
                                         const GraphCoverageMap& cover_map,
                                         vector<PairedInfoLibrary *>& libs) {
    vector<SimpleExtender *> mpPEs;
    for (size_t i = 0; i < libs.size(); ++i) {
        MatePairExtensionChooser* chooser = new MatePairExtensionChooser(
                gp.g, *libs[i], cover_map);
        SimpleExtender* mp_extender = new SimpleExtender(
                gp.g, chooser, libs[i]->GetISMax(),
                pset.loop_removal.mp_max_loops, false);
        mpPEs.push_back(mp_extender);
    }
    return mpPEs;
}

inline string LibStr(size_t count) {
    return count == 1 ? "library" : "libraries";
}


inline void ResolveRepeatsManyLibs(conj_graph_pack& gp,
		vector<PairedInfoLibrary *>& libs,
		vector<PairedInfoLibrary *>& scaff_libs,
		vector<PairedInfoLibrary *>& mp_libs,
		const vector<PathStorageInfo<Graph> >& long_reads,
		const std::string& output_dir,
		const std::string& contigs_name,
		bool traversLoops,
		boost::optional<std::string> broken_contigs) {

	INFO("Path-Extend repeat resolving tool started");

	make_dir(output_dir);
	if (cfg::get().developer_mode) {
	    make_dir(GetEtcDir(output_dir));
	}
	const pe_config::ParamSetT& pset = cfg::get().pe_params.param_set;

	INFO("Using " << libs.size() << " paired-end " << LibStr(libs.size()));
	INFO("Using " << scaff_libs.size() << " paired-end scaffolding " << LibStr(scaff_libs.size()));
	INFO("Using " << mp_libs.size() << " mate-pair " << LibStr(mp_libs.size()));
	INFO("Using " << long_reads.size() << " single read " << LibStr(long_reads.size()));
	INFO("Scaffolder is " << (pset.scaffolder_options.on ? "on" : "off"));

	ContigWriter writer(gp.g);

//make pe + long reads extenders
    vector<SimpleExtender *> usualPEs = MakePEExtenders(gp, pset, libs, false);
    vector<SimpleExtender*> long_reads_extenders = MakeLongReadsExtender(
            gp, long_reads, pset.loop_removal.max_loops, output_dir);
    vector<SimpleExtender *> shortLoopPEs = MakePEExtenders(gp, pset, libs, true);
    vector<ScaffoldingPathExtender*> scafPEs = MakeScaffoldingExtender(gp, pset, scaff_libs);
    vector<PathExtender *> all_libs(long_reads_extenders.begin(),
                                    long_reads_extenders.end());
    all_libs.insert(all_libs.end(), usualPEs.begin(), usualPEs.end() );
    all_libs.insert(all_libs.end(), shortLoopPEs.begin(), shortLoopPEs.end());
    all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
    size_t max_over = max(FindMaxOverlapedLen(libs), gp.g.k() + 100);
    CompositeExtender * mainPE = new CompositeExtender(
            gp.g, all_libs, max_over);

//extend pe + long reads
    PathExtendResolver resolver(gp.g);
    auto seeds = resolver.makeSimpleSeeds();
	seeds.SortByLength();
	seeds.ResetPathsId();
	INFO("Growing paths using paired-end and long single reads");
	auto paths = resolver.extendSeeds(seeds, *mainPE);
	DebugOutputPaths(writer, gp, output_dir, paths, "pe_overlaped_paths");
	paths.SortByLength();

    if (mp_libs.size() == 0) {
        resolver.removeOverlaps(paths, mainPE->GetCoverageMap(), max_over,
                                writer, output_dir);
        if (cfg::get().avoid_rc_connections) {
            paths.FilterInterstandBulges();
        }
        paths.FilterEmptyPaths();
        resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());
        paths.SortByLength();
        if (broken_contigs.is_initialized()) {
            OutputBrokenScaffolds(paths, (int) gp.g.k(), writer,
                                  output_dir + broken_contigs.get());
        }
        if (traversLoops) {
            INFO("Traversing tandem repeats");
            LoopTraverser loopTraverser(gp.g, mainPE->GetCoverageMap(), mainPE);
            loopTraverser.TraverseAllLoops();
            paths.SortByLength();
        }
        //writer.writePaths(paths, output_dir + contigs_name + ".fasta");
        DebugOutputPaths(writer, gp, output_dir, paths, "final_paths");
        writer.WritePathsToFASTG(paths, output_dir + contigs_name + ".fastg", output_dir + contigs_name + ".fasta");
        return;
    }
    //writer.writePaths(paths, output_dir + "pe_paths.fasta");
    if (cfg::get().pe_params.debug_output) {
        writer.WritePathsToFASTG(paths, output_dir + "pe_paths.fastg", output_dir + "pe_paths.fasta");
    }

//MP
    INFO("Adding mate-pairs");
    vector<SimpleExtender*> mpPEs = MakeMPExtenders(gp, pset,
                                                    mainPE->GetCoverageMap(),
                                                    mp_libs);
    max_over = std::max(FindMaxOverlapedLen(mp_libs), FindMaxOverlapedLen(libs));
	all_libs.clear();
	all_libs.insert(all_libs.end(), long_reads_extenders.begin(),
			long_reads_extenders.end());
	all_libs.insert(all_libs.end(), shortLoopPEs.begin(), shortLoopPEs.end());
	all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
	all_libs.insert(all_libs.end(), mpPEs.begin(), mpPEs.end());
	CompositeExtender* mp_main_pe = new CompositeExtender(gp.g, all_libs, max_over);
	INFO("Growing paths using mate-pairs");
	auto mp_paths = resolver.extendSeeds(paths, *mp_main_pe);
    resolver.removeOverlaps(mp_paths, mp_main_pe->GetCoverageMap(), max_over,
                            writer, output_dir);
    resolver.RemoveMatePairEnds(mp_paths, max_over);
	if (cfg::get().avoid_rc_connections) {
	    mp_paths.FilterInterstandBulges();
	}
	mp_paths.FilterEmptyPaths();
	//resolver.addUncoveredEdges(mp_paths, mp_main_pe->GetCoverageMap());
	mp_paths.SortByLength();
	DebugOutputPaths(writer, gp, output_dir, mp_paths, "mp_final_paths");
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
    CompositeExtender* last_extender = new CompositeExtender(gp.g, all_libs, max_over);
    auto last_paths = resolver.extendSeeds(mp_paths, *last_extender);
    resolver.removeOverlaps(last_paths, last_extender->GetCoverageMap(), max_over,
                    writer, output_dir);

    if (cfg::get().avoid_rc_connections) {
        last_paths.FilterInterstandBulges();
    }
    last_paths.FilterEmptyPaths();
    resolver.addUncoveredEdges(last_paths, last_extender->GetCoverageMap());
    last_paths.SortByLength();

//Traverse loops
    if (traversLoops) {
        INFO("Traversing tandem repeats");
        LoopTraverser loopTraverser(gp.g, last_extender->GetCoverageMap(), last_extender);
        loopTraverser.TraverseAllLoops();
        last_paths.SortByLength();
    }

//result
    if (broken_contigs.is_initialized()) {
        OutputBrokenScaffolds(last_paths, (int) gp.g.k(), writer,
                              output_dir + broken_contigs.get());
    }
    DebugOutputPaths(writer, gp, output_dir, last_paths, "last_paths");

    writer.WritePathsToFASTG(last_paths, output_dir + contigs_name + ".fastg", output_dir + contigs_name + ".fasta");
    //writer.writePaths(last_paths, output_dir + contigs_name+ "_temp" + ".fasta");
    last_paths.DeleteAllPaths();
    seeds.DeleteAllPaths();
    mp_paths.DeleteAllPaths();
    paths.DeleteAllPaths();

    INFO("Path-Extend repeat resolving tool finished");
    //TODO:DELETE ALL!!!!
}

inline PairedInfoLibrary* MakeNewLib(conj_graph_pack::graph_t& g,
                              const PairedInfoIndicesT& paired_index,
                              size_t index) {
    size_t read_length = cfg::get().ds.reads[index].data().read_length;
    size_t is = (size_t) cfg::get().ds.reads[index].data().mean_insert_size;
    int is_min = (int) cfg::get().ds.reads[index].data().insert_size_left_quantile;
    int is_max = (int) cfg::get().ds.reads[index].data().insert_size_right_quantile;
    int var = (int) cfg::get().ds.reads[index].data().insert_size_deviation;
    bool is_mp = cfg::get().ds.reads[index].type() == io::LibraryType::MatePairs;
    PairedInfoLibrary* lib = new PairedInfoLibrary(
            cfg::get().K, g, read_length, is, is_min > 0.0 ? size_t(is_min) : 0,
            is_max > 0.0 ? size_t(is_max) : 0, size_t(var), paired_index[index], is_mp,
            cfg::get().ds.reads[index].data().insert_size_distribution);
    return lib;
}

inline void DeleteLibs(vector<PairedInfoLibrary *>& libs) {
    for (size_t j = 0; j < libs.size(); ++j) {
        delete libs[j];
    }
}

inline bool InsertSizeCompare(const PairedInfoLibrary* lib1,
                       const PairedInfoLibrary* lib2) {
    return lib1->GetISMax() < lib2->GetISMax();
}

inline void ResolveRepeatsPe(conj_graph_pack& gp,
                      const vector<PathStorageInfo<Graph> >& long_reads,
                      const std::string& output_dir,
                      const std::string& contigs_name, bool traverseLoops,
                      boost::optional<std::string> broken_contigs,
                      bool use_auto_threshold = true) {
    vector<PairedInfoLibrary*> rr_libs;
    vector<PairedInfoLibrary*> mp_libs;
    vector<PairedInfoLibrary*> scaff_libs;
    for (size_t i = 0; i < gp.clustered_indices.size(); ++i) {
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd
                /*|| cfg::get().ds.reads[indexs[i]].type()
                        == io::LibraryType::MatePairs*/ && cfg::get().ds.reads[i].data().mean_insert_size > 0.0) {
            PairedInfoLibrary* lib = MakeNewLib(gp.g, gp.clustered_indices, i);
            if (use_auto_threshold) {
                lib->SetSingleThreshold(cfg::get().ds.reads[i].data().pi_threshold);
                INFO("Threshold for library #" << i << " is " << cfg::get().ds.reads[i].data().pi_threshold);
            }
            rr_libs.push_back(lib);
        }
    }
    for (size_t i = 0; i < gp.paired_indices.size(); ++i) {
		if (cfg::get().ds.reads[i].type()
				== io::LibraryType::MatePairs&& cfg::get().ds.reads[i].data().mean_insert_size > 0.0) {
			PairedInfoLibrary* lib = MakeNewLib(gp.g, gp.paired_indices/*paired_index*/, i);
			mp_libs.push_back(lib);
		}
	}
    std::sort(rr_libs.begin(), rr_libs.end(), InsertSizeCompare);
    if (cfg::get().use_scaffolder
            && cfg::get().pe_params.param_set.scaffolder_options.on) {
        for (size_t i = 0; i < gp.scaffolding_indices.size(); ++i) {
            if ((cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd)
             && cfg::get().ds.reads[i].data().mean_insert_size > 0.0) {
                PairedInfoLibrary* lib = MakeNewLib(gp.g, gp.scaffolding_indices, i);
                scaff_libs.push_back(lib);
            }
        }
    }
    std::sort(scaff_libs.begin(), scaff_libs.end(), InsertSizeCompare);
    ResolveRepeatsManyLibs(gp, rr_libs, scaff_libs, mp_libs, long_reads, output_dir,
                           contigs_name, traverseLoops, broken_contigs);
    DeleteLibs(rr_libs);
    DeleteLibs(mp_libs);
    DeleteLibs(scaff_libs);
}


} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
