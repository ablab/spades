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
#include "single_threshold_finder.hpp"
#include "long_read_storage.hpp"
#include "split_graph_pair_info.hpp"

namespace path_extend {

using namespace debruijn_graph;

size_t FindMaxOverlapedLen(const vector<PairedInfoLibraries>& libes) {
    size_t max = 0;
    for (size_t i = 0; i < libes.size(); ++i) {
        for (size_t j = 0; j < libes[i].size(); ++j) {
            size_t overlap = libes[i].at(j)->insert_size_
                    + libes[i].at(j)->is_variation_;
            if (overlap > max) {
                max = overlap;
            }
        }
    }
    return max;
}

size_t GetMinInsertSize(const vector<PairedInfoLibraries>& libs) {
    int min = 0;
    size_t index = 0;
    while (index < libs.size() && libs[index].size() == 0) {
        index++;
    }
    if (index == libs.size()) {
        return 0;
    }
    min = libs[index][0]->insert_size_ - libs[index][0]->read_size_ / 2;
    for (size_t i = 0; i < libs.size(); ++i) {
        for (size_t j = 0; j < libs[i].size(); ++j) {
            int overlap = libs[i][j]->insert_size_ - libs[i][j]->read_size_ / 2;
            if (overlap < min && overlap > 0) {
                min = overlap;
            }
        }
    }
    return (size_t) min;
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

double GetWeightThreshold(const PairedInfoLibraries& lib){
	return lib[0]->is_mate_pair_ ?
	        cfg::get().pe_params.param_set.mate_pair_options.select_options.weight_threshold :
	        cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold;
}

double GetSingleThreshold(const PairedInfoLibraries& lib){
    return lib[0]->is_mate_pair_ ?
            *cfg::get().pe_params.param_set.mate_pair_options.select_options.single_threshold :
            *cfg::get().pe_params.param_set.extension_options.select_options.single_threshold;
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
		const pe_config::ParamSetT& pset, vector<PairedInfoLibraries>& libs,
		bool investigateShortLoops) {
	vector<WeightCounter*> wcs;
	for (size_t i = 0; i < libs.size(); ++i) {
		wcs.push_back(new PathCoverWeightCounter(gp.g, libs[i],
				GetWeightThreshold(libs[i]), GetSingleThreshold(libs[i])));
	}
	vector<SimpleExtender *> usualPEs;
	for (size_t i = 0; i < libs.size(); ++i) {
		wcs[i]->setNormalizeWeight(pset.normalize_weight);
		wcs[i]->setNormalizeWightByCoverage(pset.normalize_by_coverage);
		double
				priory_coef =
						libs[i].at(0)->is_mate_pair_ ? pset.mate_pair_options.select_options.priority_coeff
								: pset.extension_options.select_options.priority_coeff;
		SimpleExtensionChooser * extensionChooser = new SimpleExtensionChooser(
				gp.g, wcs[i], priory_coef);
		usualPEs.push_back(new SimpleExtender(gp.g,
				pset.loop_removal.max_loops, extensionChooser,
				investigateShortLoops));
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
        BidirectionalPath* conj_path = new BidirectionalPath(
                new_path->Conjugate());
        new_path->SetWeight(path.getWeight());
        conj_path->SetWeight(path.getWeight());
        result.AddPair(new_path, conj_path);
    }
    INFO("== Long reads paths " << result.size() << " == ");
    for (size_t index = 0; index < result.size(); ++index) {
        DEBUG("Long contig " << index);
        result.Get(index)->Print();
    }INFO("==== ");
}

vector<SimpleExtender*> MakeLongReadsExtender(const conj_graph_pack& gp,
                           const vector<PathStorageInfo<Graph> >& long_reads,
                           size_t max_loops) {
    vector<SimpleExtender*> result;
    for (size_t i = 0; i < long_reads.size(); ++i) {
        PathContainer paths;
        AddPathsToContainer(gp, long_reads[i].GetPaths(), 1, paths);
        ExtensionChooser * longReadEC = new LongReadsExtensionChooser(
                gp.g, paths, long_reads[i].GetFilteringThreshold(),
                long_reads[i].GetPriorityThreshold());
        SimpleExtender * longReadPathExtender = new SimpleExtender(
                gp.g, max_loops, longReadEC, true);
        result.push_back(longReadPathExtender);
    }
    return result;
}

void ResolveRepeatsManyLibs(conj_graph_pack& gp,
		vector<PairedInfoLibraries>& libs,
		vector<PairedInfoLibraries>& scafolding_libs,
		const vector<PathStorageInfo<Graph> >& long_reads,
		const std::string& output_dir,
		const std::string& contigs_name,
		bool traversLoops,
		boost::optional<std::string> broken_contigs) {

	INFO("Path extend repeat resolving tool started");
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
				(int) (pset.scaffolder_options.max_must_overlap * gp.g.k()),
				(int) (pset.scaffolder_options.max_can_overlap * gp.g.k()),
				pset.scaffolder_options.short_overlap);

	vector<ScaffoldingPathExtender*> scafPEs;
	for (size_t i = 0; i < scafolding_libs.size(); ++i){
		scaf_wcs.push_back(new ReadCountWeightCounter(gp.g, scafolding_libs[i]));
		ScaffoldingExtensionChooser * scafExtensionChooser = new ScaffoldingExtensionChooser(gp.g, scaf_wcs[i],
							scafolding_libs[i].at(0)->is_mate_pair_? pset.mate_pair_options.select_options.priority_coeff : pset.extension_options.select_options.priority_coeff);
		scafPEs.push_back(new ScaffoldingPathExtender(gp.g, pset.loop_removal.max_loops, scafExtensionChooser, gapJoiner));
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
	INFO("Extending seeds");
	auto paths = resolver.extendSeeds(seeds, *mainPE);
	if (cfg::get().pe_params.output.write_overlaped_paths) {
		writer.writePaths(paths, GetEtcDir(output_dir) + "overlaped_paths.fasta");
	}

	DebugOutputPaths(writer, gp, output_dir, paths, "overlaped_paths");
    size_t max_over = FindMaxOverlapedLen(libs);
    if (pset.filter_options.remove_overlaps) {
        resolver.removeOverlaps(paths, mainPE->GetCoverageMap(), max_over,
                                writer, output_dir);
    }
    paths.FilterEmptyPaths();
	paths.CheckSymmetry();
    resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());
	paths.SortByLength();
    DebugOutputPaths(writer, gp, output_dir, paths, "final_paths");

	if (broken_contigs.is_initialized()) {
	    OutputBrokenScaffolds(paths, gp.g.k(), writer, output_dir + broken_contigs.get());
	}

    if (traversLoops) {
        INFO("Traversing tandem repeats");
        LoopTraverser loopTraverser(gp.g, paths, mainPE->GetCoverageMap(), mainPE);
        loopTraverser.TraverseAllLoops();
        paths.SortByLength();
    }

    writer.writePaths(paths, output_dir + contigs_name);

    INFO("Path extend repeat resolving tool finished");
    for (size_t i = 0; i < all_libs.size(); ++i){
        delete all_libs[i];
    }
}

PairedInfoLibrary* MakeNewLib(conj_graph_pack::graph_t& g,
		vector<PairedIndexT*>& paired_index, vector<size_t>& indexs,
		size_t index) {

	size_t read_length = cfg::get().ds.reads[indexs[index]].data().read_length;
	double is = cfg::get().ds.reads[indexs[index]].data().mean_insert_size;
	double var = cfg::get().ds.reads[indexs[index]].data().insert_size_deviation;
	bool is_mp = cfg::get().ds.reads[indexs[index]].type() == io::LibraryType::MatePairs;
	PairedInfoLibrary* lib = new PairedInfoLibrary(cfg::get().K, g, read_length,
			is, var, *paired_index[index], is_mp);
	return lib;
}

void DeleteLibs(vector<PairedInfoLibraries>& libs) {
    for (size_t j = 0; j < libs.size(); ++j) {
        for (size_t i = 0; i < libs[j].size(); ++i) {
            delete libs[j].at(i);
        }
    }
}

void SetOldThreshold(PairedInfoLibrary* lib, size_t index, size_t split_edge_length) {
	INFO("Searching for paired info threshold for lib #"
						<< index << " (IS = " << lib->insert_size_ << ",  DEV = " << lib->is_variation_ << ")");

	SingleThresholdFinder finder(lib->insert_size_ - 2 * lib->is_variation_, lib->insert_size_ + 2 * lib->is_variation_, lib->read_size_);
	double threshold = finder.find_threshold(index);

	INFO("Paired info threshold is " << threshold);
	lib->SetSingleThreshold(threshold);
}

void SetNewThreshold(conj_graph_pack& gp, PairedInfoLibrary* lib, size_t index,
                     size_t split_edge_length) {
    SplitGraphPairInfo split_graph(gp, *lib, index, 99);
    INFO("Calculating paired info threshold");
    split_graph.ProcessReadPairs();
    double threshold = split_graph.FindThreshold(
            split_edge_length, lib->insert_size_ - 2 * lib->is_variation_,
            lib->insert_size_ + 2 * lib->is_variation_);
    lib->SetSingleThreshold(threshold);
    //double tr = *cfg::get().pe_params.param_set.extension_options.select_options.single_threshold;
    //INFO("threshold taken from config - "  << tr);
    //lib->SetSingleThreshold(tr);
}

bool InsertSizeCompare(const PairedInfoLibraries& lib1,
                       const PairedInfoLibraries& lib2) {
    if (lib1.size() < 1 or lib2.size() < 1) {
        return true;
    }
    return lib1[0]->insert_size_ < lib2[0]->insert_size_;
}

void ResolveRepeatsPe(conj_graph_pack& gp, vector<PairedIndexT*>& paired_index,
                      vector<PairedIndexT*>& scaff_index,
                      vector<size_t>& indexs,
                      const vector<PathStorageInfo<Graph> >& long_reads,
                      const std::string& output_dir,
                      const std::string& contigs_name, bool traverseLoops,
                      boost::optional<std::string> broken_contigs,
                      bool use_auto_threshold = true) {

    const pe_config::ParamSetT& pset = cfg::get().pe_params.param_set;
    vector<PairedInfoLibraries> rr_libs;
    vector<PairedInfoLibraries> scaff_libs;
    for (size_t i = 0; i < paired_index.size(); ++i) {
        if (cfg::get().ds.reads[indexs[i]].type() == io::LibraryType::PairedEnd
                || cfg::get().ds.reads[indexs[i]].type()
                        == io::LibraryType::MatePairs) {
            PairedInfoLibrary* lib = MakeNewLib(gp.g, paired_index, indexs, i);
            if (use_auto_threshold) {
                SetNewThreshold(gp, lib, indexs[i], pset.split_edge_length);
            }
            PairedInfoLibraries libs;
            libs.push_back(lib);
            rr_libs.push_back(libs);
        }
    }
    std::sort(rr_libs.begin(), rr_libs.end(), InsertSizeCompare);
    for (size_t i = 0; i < scaff_index.size(); ++i) {
        if (cfg::get().ds.reads[indexs[i]].type() == io::LibraryType::PairedEnd
                || cfg::get().ds.reads[indexs[i]].type()
                        == io::LibraryType::MatePairs) {
            PairedInfoLibrary* lib = MakeNewLib(gp.g, scaff_index, indexs, i);
            PairedInfoLibraries libs;
            libs.push_back(lib);
            scaff_libs.push_back(libs);
        }
    }
    std::sort(scaff_libs.begin(), scaff_libs.end(), InsertSizeCompare);
    ResolveRepeatsManyLibs(gp, rr_libs, scaff_libs, long_reads, output_dir,
                           contigs_name, traverseLoops, broken_contigs);
    for (size_t i = 0; i < rr_libs.size(); ++i) {
        INFO("Insert size for lib# " << i << " = " << rr_libs[i][0]->insert_size_);
    }
    DeleteLibs(rr_libs);
    DeleteLibs(scaff_libs);
}


} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
