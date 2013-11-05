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
#include "pe_io.hpp"
#include "path_visualizer.hpp"
#include "loop_traverser.hpp"
#include "long_read_storage.hpp"
#include "split_graph_pair_info.hpp"
#include "mate_pair_scaffolding.hpp"
#include "next_path_searcher.hpp"


namespace path_extend {

using namespace debruijn_graph;

size_t FindMaxOverlapedLen(const vector<PairedInfoLibraries>& libes) {
    size_t max = 0;
    for (size_t i = 0; i < libes.size(); ++i) {
        for (size_t j = 0; j < libes[i].size(); ++j) {
            size_t overlap = libes[i].at(j)->GetISMax();
            if (overlap > max) {
                max = overlap;
            }
        }
    }
    return max;
}

size_t GetMinInsertSize(const vector<PairedInfoLibraries>& libs) {
    size_t min = 0;
    size_t index = 0;
    while (index < libs.size() && libs[index].size() == 0) {
        index++;
    }
    if (index == libs.size()) {
        return 0;
    }
    min = libs[index][0]->GetISMax();
    for (size_t i = 0; i < libs.size(); ++i) {
        for (size_t j = 0; j < libs[i].size(); ++j) {
            int overlap = (int) libs[i][j]->GetISMax() - (int) libs[i][j]->GetReadSize() / 2;
            if (overlap < (int) min && overlap > 0) {
                min = (size_t) overlap;
            }
        }
    }
    return min;
}

string GetEtcDir(const std::string& output_dir) {
    return output_dir + cfg::get().pe_params.etc_dir + "/";
}

void DebugOutputPaths(const ContigWriter& writer, const conj_graph_pack& gp,
                      const std::string& output_dir, PathContainer& paths,
                      const string& name) {
    PathInfoWriter path_writer;
    PathVisualizer visualizer(gp.g.k());
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

void DebugOutputEdges(const ContigWriter& writer, const conj_graph_pack& gp,
                      const std::string& output_dir, const string& name) {
    if (!cfg::get().pe_params.debug_output) {
        return;
    }
    PathVisualizer visualizer(gp.g.k());
    string etcDir = GetEtcDir(output_dir);
    DEBUG("debug_output " << cfg::get().pe_params.debug_output);
    if (cfg::get().pe_params.viz.print_paths) {
        writer.writeEdges(etcDir + name + ".fasta");
        visualizer.writeGraphSimple(gp, etcDir + name + ".dot", name);
    }
}

double GetWeightThreshold(const PairedInfoLibraries& lib,
                          const pe_config::ParamSetT& pset) {
    return lib[0]->IsMp() ?
            pset.mate_pair_options.weight_threshold :
            pset.extension_options.weight_threshold;
}

double GetSingleThreshold(const PairedInfoLibraries& lib,
                          const pe_config::ParamSetT& pset) {
    return lib[0]->IsMp() ?
            pset.mate_pair_options.single_threshold :
            pset.extension_options.single_threshold;
}

double GetPriorityCoeff(const PairedInfoLibraries& lib,
                        const pe_config::ParamSetT& pset) {
    return lib[0]->IsMp() ?
            pset.mate_pair_options.priority_coeff :
            pset.extension_options.priority_coeff;
}

string MakeNewName(const std::string& contigs_name,
		const std::string& subname){
	return contigs_name.substr(0, contigs_name.rfind(".fasta")) + "_" + subname + ".fasta";
}


void OutputBrokenScaffolds(PathContainer& paths, int k,
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
    writer.writePaths(breaker.container(), filename);
}

vector<SimpleExtender *> MakeExtenders(const conj_graph_pack& gp,
                                       const pe_config::ParamSetT& pset,
                                       vector<PairedInfoLibraries>& libs,
                                       bool investigate_loops) {
    vector<WeightCounter*> wcs;
    for (size_t i = 0; i < libs.size(); ++i) {
        wcs.push_back(
                new PathCoverWeightCounter(gp.g, libs[i],
                                           GetWeightThreshold(libs[i], pset),
                                           GetSingleThreshold(libs[i], pset)));
    }
    vector<SimpleExtender *> usualPEs;
    for (size_t i = 0; i < libs.size(); ++i) {
        wcs[i]->setNormalizeWeight(pset.normalize_weight);
        double priory_coef = GetPriorityCoeff(libs[i], pset);
        SimpleExtensionChooser * extension = new SimpleExtensionChooser(
                gp.g, wcs[i], priory_coef);
        usualPEs.push_back(
                new SimpleExtender(gp.g, pset.loop_removal.max_loops, extension,
                                   investigate_loops));
    }
    return usualPEs;
}

void AddPathsToContainer(const conj_graph_pack& gp,
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

        //FIXME: why weights are int/double in different path structures?
        new_path->SetWeight((double) path.getWeight());
        conj_path->SetWeight((double) path.getWeight());
        result.AddPair(new_path, conj_path);
    }
    DEBUG("Long reads paths " << result.size() << " == ");
}

vector<SimpleExtender*> MakeLongReadsExtender(
        const conj_graph_pack& gp,
        const vector<PathStorageInfo<Graph> >& long_reads, size_t max_loops) {
    vector<SimpleExtender*> result;
    for (size_t i = 0; i < long_reads.size(); ++i) {
        if (long_reads[i].GetPaths().size() == 0) {
            continue;
        }
        PathContainer paths;
        AddPathsToContainer(gp, long_reads[i].GetPaths(), 1, paths);
        ExtensionChooser * longReadEC = new LongReadsExtensionChooser(
                gp.g, paths, long_reads[i].GetFilteringThreshold(),
                long_reads[i].GetWeightPriorityThreshold(),
                long_reads[i].GetUniqueEdgePriorityThreshold());
        SimpleExtender * longReadExtender = new SimpleExtender(gp.g, max_loops,
                                                               longReadEC,
                                                               true);
        result.push_back(longReadExtender);
    }
    return result;
}
void TestIdealInfo(conj_graph_pack& gp) {
    map<int, size_t> distr;
    distr[220] = 1;
    IdealPairInfoCounter counter(gp.g, 220, 221, 100, distr);
    /*for (int i = 0; i < 220; ++i){
        size_t edge_len = 100;
        counter.IdealPairedInfo(edge_len, edge_len, edge_len + i);
    }*/
    size_t edge_len = 100;
    double w1 = counter.IdealPairedInfo(edge_len, edge_len, edge_len);
    size_t edge_len1_1 = 50;
    size_t edge_len1_2 = 50;
    double w2_1 = counter.IdealPairedInfo(edge_len, edge_len1_1, edge_len);
    double w2_2 = counter.IdealPairedInfo(edge_len, edge_len1_2, edge_len + edge_len1_2);
    size_t edge_len2_1 = 20;
    size_t edge_len2_2 = 30;
    size_t edge_len2_3 = 50;
    double w3_1 = counter.IdealPairedInfo(edge_len, edge_len2_1, edge_len);
    double w3_2 = counter.IdealPairedInfo(edge_len, edge_len2_2, edge_len + edge_len2_1);
    double w3_3 = counter.IdealPairedInfo(edge_len, edge_len2_3, edge_len + edge_len2_1 + edge_len2_2);
    DEBUG("TEST " << w1 << " " << w2_1 + w2_2 << " " << w3_1 + w3_2 + w3_3);
}
void ResolveRepeatsManyLibs(conj_graph_pack& gp,
		vector<PairedInfoLibraries>& libs,
		vector<PairedInfoLibraries>& scafolding_libs,
		vector<PairedInfoLibraries>& mp_libs,
		const vector<PathStorageInfo<Graph> >& long_reads,
		const std::string& output_dir,
		const std::string& contigs_name,
		bool traversLoops,
		boost::optional<std::string> broken_contigs) {

	INFO("Path extend repeat resolving tool started");
	TestIdealInfo(gp);
	make_dir(output_dir);
	if (cfg::get().developer_mode) {
	    make_dir(GetEtcDir(output_dir));
	}
	const pe_config::ParamSetT& pset = cfg::get().pe_params.param_set;

	INFO("Using " << libs.size() << " paired lib(s)");
	INFO("Using " << scafolding_libs.size() << " scaffolding libs");
	INFO("Scaffolder is " << (pset.scaffolder_options.on ? "on" : "off"));

	ContigWriter writer(gp.g);
	DebugOutputEdges(writer, gp, output_dir, "before_resolve");
    vector<WeightCounter*> scaf_wcs;
    GapJoiner * gapJoiner = new HammingGapJoiner(gp.g,
				pset.scaffolder_options.min_gap_score,
				(int) (pset.scaffolder_options.max_must_overlap * (double) gp.g.k()),
				(int) (pset.scaffolder_options.max_can_overlap * (double) gp.g.k()),
				pset.scaffolder_options.short_overlap);

	vector<ScaffoldingPathExtender*> scafPEs;
	for (size_t i = 0; i < scafolding_libs.size(); ++i){
		scaf_wcs.push_back(new ReadCountWeightCounter(gp.g, scafolding_libs[i]));
        ScaffoldingExtensionChooser * scaff_chooser =
                new ScaffoldingExtensionChooser(
                        gp.g, scaf_wcs[i], GetPriorityCoeff(libs[i], pset));
		scafPEs.push_back(new ScaffoldingPathExtender(gp.g, pset.loop_removal.max_loops, scaff_chooser, gapJoiner));
	}

	PathExtendResolver resolver(gp.g);
    auto seeds = resolver.makeSimpleSeeds();
    vector<SimpleExtender *> usualPEs = MakeExtenders(gp, pset, libs, false);
    vector<SimpleExtender*> long_reads_extenders = MakeLongReadsExtender(
            gp, long_reads, pset.loop_removal.max_loops);
    vector<SimpleExtender *> shortLoopPEs = MakeExtenders(gp, pset, libs, true);
    vector<PathExtender *> all_libs(usualPEs.begin(), usualPEs.end());
    all_libs.insert(all_libs.end(), long_reads_extenders.begin(),
                    long_reads_extenders.end());
    all_libs.insert(all_libs.end(), shortLoopPEs.begin(), shortLoopPEs.end());
    all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
    CoveringPathExtender * mainPE = new CompositeExtender(
            gp.g, pset.loop_removal.max_loops, all_libs);

	seeds.SortByLength();
	seeds.ResetPathsId();
	INFO("Growing paths");
	auto paths = resolver.extendSeeds(seeds, *mainPE);
	if (cfg::get().pe_params.output.write_overlaped_paths) {
		writer.writePaths(paths, GetEtcDir(output_dir) + "overlaped_paths.fasta");
	}

	DebugOutputPaths(writer, gp, output_dir, paths, "overlaped_paths");
    size_t max_over = FindMaxOverlapedLen(libs);
    /*NextPathSearcher searcher (gp.g, mainPE->GetCoverageMap(), 4000);
    for (size_t i = 0; i < paths.size(); ++i) {
    	searcher.FindNextPaths(paths.Get(i));
    }*/
   /* if (pset.remove_overlaps) {
        resolver.removeOverlaps(paths, mainPE->GetCoverageMap(), max_over,
                                writer, output_dir);
    }*/
    paths.FilterEmptyPaths();
	paths.CheckSymmetry();
    resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());
	paths.SortByLength();
    DebugOutputPaths(writer, gp, output_dir, paths, "final_paths");
    if (broken_contigs.is_initialized()) {
        OutputBrokenScaffolds(paths, (int) gp.g.k(), writer, output_dir + broken_contigs.get());
    }
    writer.writePaths(paths, output_dir + "final_paths.fasta");


    //MP
    DEBUG("mate pair path-extend started");
    vector<SimpleExtender *> mpPEs;
	for (size_t i = 0; i < mp_libs.size(); ++i) {
		MatePairExtensionChooser* chooser = new MatePairExtensionChooser(gp.g,
				*mp_libs[i][0], mainPE->GetCoverageMap());
		SimpleExtender* mp_extender = new SimpleExtender(gp.g,
				pset.loop_removal.max_loops, chooser);
		mpPEs.push_back(mp_extender);
	}
	ExtensionChooser * longReadEC =
			new LongReadsExtensionChooser(
					gp.g,
					paths,
					cfg::get().pe_params.long_reads.coverage_base_rr.filtering,
					cfg::get().pe_params.long_reads.single_reads.weight_priority,
					cfg::get().pe_params.long_reads.coverage_base_rr.unique_edge_priority);
    SimpleExtender * longReadExtender = new SimpleExtender(gp.g,
			pset.loop_removal.max_loops, longReadEC, true);
	all_libs.clear();
	all_libs.push_back(longReadExtender);
	all_libs.insert(all_libs.end(), long_reads_extenders.begin(),
			long_reads_extenders.end());
	all_libs.insert(all_libs.end(), shortLoopPEs.begin(), shortLoopPEs.end());
	all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
	all_libs.insert(all_libs.end(), mpPEs.begin(), mpPEs.end());
	CompositeExtender* mp_main_pe = new CompositeExtender(gp.g, pset.loop_removal.max_loops, all_libs);
	INFO("Growing mp paths");

	auto result_paths = resolver.extendSeeds(paths, *mp_main_pe);
	if (pset.remove_overlaps) {
		resolver.removeOverlaps(result_paths, mp_main_pe->GetCoverageMap(), max_over,
				writer, output_dir);
		resolver.RemoveMatePairEnds(result_paths, max_over);
	}
	result_paths.FilterEmptyPaths();
	result_paths.CheckSymmetry();
	resolver.addUncoveredEdges(result_paths, mp_main_pe->GetCoverageMap());
	result_paths.SortByLength();
	DebugOutputPaths(writer, gp, output_dir, result_paths, "mp_final_paths");


    INFO("End mp libs");
    writer.writePaths(result_paths, output_dir + "final_paths_with_mp.fasta");
	//MP end

    //pe again
    all_libs.clear();
    all_libs.push_back(longReadExtender);
    all_libs.insert(all_libs.end(), long_reads_extenders.begin(),
                    long_reads_extenders.end());
    all_libs.insert(all_libs.end(), shortLoopPEs.begin(), shortLoopPEs.end());
    all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
    CompositeExtender* last_extender = new CompositeExtender(gp.g, pset.loop_removal.max_loops, all_libs);
    auto last_paths = resolver.extendSeeds(result_paths, *last_extender);
    resolver.removeOverlaps(last_paths, last_extender->GetCoverageMap(), max_over,
                    writer, output_dir);
    last_paths.FilterEmptyPaths();
    resolver.addUncoveredEdges(last_paths, last_extender->GetCoverageMap());
    last_paths.SortByLength();

    if (traversLoops) {
        INFO("Traversing tandem repeats");
        LoopTraverser loopTraverser(gp.g, last_extender->GetCoverageMap(), last_extender);
        loopTraverser.TraverseAllLoops();
        last_paths.SortByLength();
    }
    DebugOutputPaths(writer, gp, output_dir, last_paths, "last_paths");
    writer.writePaths(last_paths, output_dir + contigs_name);


    INFO("Path extend repeat resolving tool finished");
    for (size_t i = 0; i < all_libs.size(); ++i){
        delete all_libs[i];
    }
    delete mp_main_pe;
}

PairedInfoLibrary* MakeNewLib(conj_graph_pack::graph_t& g,
		vector<PairedIndexT*>& paired_index, size_t initial_index,
		size_t index) {

	size_t read_length = cfg::get().ds.reads[initial_index].data().read_length;
	double is_min = cfg::get().ds.reads[initial_index].data().insert_size_left_quantile;
	double is_max = cfg::get().ds.reads[initial_index].data().insert_size_right_quantile;
	double var =  cfg::get().ds.reads[initial_index].data().insert_size_deviation;
	bool is_mp = cfg::get().ds.reads[initial_index].type() == io::LibraryType::MatePairs;
	PairedInfoLibrary* lib =new PairedInfoLibrary(
			cfg::get().K, g, read_length, is_min > 0.0 ? size_t(is_min) : 0,
					is_max > 0.0 ? size_t(is_max) : 0, var, *paired_index[index],
							is_mp,
							cfg::get().ds.reads[initial_index].data().insert_size_distribution);

	return lib;
}

void DeleteLibs(vector<PairedInfoLibraries>& libs) {
    for (size_t j = 0; j < libs.size(); ++j) {
        for (size_t i = 0; i < libs[j].size(); ++i) {
            delete libs[j].at(i);
        }
    }
}

bool InsertSizeCompare(const PairedInfoLibraries& lib1,
                       const PairedInfoLibraries& lib2) {
    if (lib1.size() < 1 or lib2.size() < 1) {
        return true;
    }
    return lib1[0]->GetISMax() < lib2[0]->GetISMax();
}

void AnalyzeMatePairInfo(const conj_graph_pack& gp, const PairedInfoLibrary& lib, PathContainer& paths){
    PathsPairInfoContainer container(gp, paths, lib);
    container.FillPairInfo();
}

void ResolveRepeatsPe(conj_graph_pack& gp, vector<PairedIndexT*>& paired_index,
                      vector<PairedIndexT*>& scaff_index,
                      vector<size_t>& indexs, vector<PairedIndexT*>& paired_indexes_not_clust,
                      const vector<PathStorageInfo<Graph> >& long_reads,
                      const std::string& output_dir,
                      const std::string& contigs_name, bool traverseLoops,
                      boost::optional<std::string> broken_contigs,
                      bool use_auto_threshold = true) {
    vector<PairedInfoLibraries> rr_libs;
    vector<PairedInfoLibraries> mp_libs;
    vector<PairedInfoLibraries> scaff_libs;
    for (size_t i = 0; i < paired_index.size(); ++i) {
        if (cfg::get().ds.reads[indexs[i]].type() == io::LibraryType::PairedEnd
                /*|| cfg::get().ds.reads[indexs[i]].type()
                        == io::LibraryType::MatePairs*/) {
            PairedInfoLibrary* lib = MakeNewLib(gp.g, paired_index, indexs[i], i);
            if (use_auto_threshold) {
                lib->SetSingleThreshold(cfg::get().ds.reads[indexs[i]].data().pi_threshold);
                INFO("Threshold for library # " << indexs[i] << " is " << cfg::get().ds.reads[indexs[i]].data().pi_threshold);
            }
            PairedInfoLibraries libs;
            libs.push_back(lib);
            rr_libs.push_back(libs);
        }
    }
    for (size_t i = 0; i < paired_index.size(); ++i) {
		if (cfg::get().ds.reads[indexs[i]].type()
				== io::LibraryType::MatePairs) {
			PairedInfoLibrary* lib = MakeNewLib(gp.g, paired_indexes_not_clust/*paired_index*/,
					indexs[i], i);
			PairedInfoLibraries libs;
			libs.push_back(lib);
			mp_libs.push_back(libs);
		}
	}
    std::sort(rr_libs.begin(), rr_libs.end(), InsertSizeCompare);
    for (size_t i = 0; i < scaff_index.size(); ++i) {
        if (cfg::get().ds.reads[indexs[i]].type() == io::LibraryType::PairedEnd
               /* || cfg::get().ds.reads[indexs[i]].type()
                        == io::LibraryType::MatePairs*/) {
            PairedInfoLibrary* lib = MakeNewLib(gp.g, scaff_index, indexs[i], i);
            PairedInfoLibraries libs;
            libs.push_back(lib);
            scaff_libs.push_back(libs);
        }
    }
    std::sort(scaff_libs.begin(), scaff_libs.end(), InsertSizeCompare);
    ResolveRepeatsManyLibs(gp, rr_libs, scaff_libs, mp_libs, long_reads, output_dir,
                           contigs_name, traverseLoops, broken_contigs);
    /*for (size_t i = 0; i < paired_indexes_not_clust.size(); ++i) {
        if (cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs) {
            PairedInfoLibrary* lib = MakeNewLib(gp.g, paired_indexes_not_clust,
                                                i, i);
            AnalyzeMatePairInfo(gp, *lib, resolved_paths);
            delete lib;
        }
    }*/
    DeleteLibs(rr_libs);
    DeleteLibs(scaff_libs);
}


} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
