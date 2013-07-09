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
#include "../include/omni/edges_position_handler.hpp"
#include "split_graph_pair_info.hpp"
#include "long_read_mapper.hpp"

namespace path_extend {

using namespace debruijn_graph;

size_t find_max_overlaped_len(vector<PairedInfoLibraries>& libes) {
	size_t max = 0;
	for (size_t i = 0; i < libes.size(); ++i) {
		for (size_t j = 0; j < libes[i].size(); ++j) {
			size_t overlap = libes[i][j]->insert_size_ + libes[i][j]->is_variation_;
			if (overlap > max) {
				max = overlap;
			}
		}
	}
	return max;
}

size_t getMinInsertSize(const vector<PairedInfoLibraries>& libs) {

	int min = 0;
	size_t index = 0;
	while (index < libs.size() && libs[index].size() == 0){
		index++;
	}
	if (index == libs.size()){
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

string get_etc_dir(const std::string& output_dir){
	return output_dir + cfg::get().pe_params.etc_dir + "/";
}

void debug_output_paths(ContigWriter& writer, conj_graph_pack& gp,
		const std::string& output_dir, PathContainer& paths,
		const string& name) {
	PathInfoWriter path_writer;
	PathVisualizer visualizer(gp.g.k());
    string etcDir = get_etc_dir(output_dir);
	if (!cfg::get().pe_params.debug_output){
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

void debug_output_edges(ContigWriter& writer, conj_graph_pack& gp,
		const std::string& output_dir, const string& name) {
	if (!cfg::get().pe_params.debug_output){
        return;
    }
	PathVisualizer visualizer(gp.g.k());
	string etcDir = get_etc_dir(output_dir);
	DEBUG("debug_output " << cfg::get().pe_params.debug_output);
	if (cfg::get().pe_params.viz.print_paths) {
		writer.writeEdges(etcDir + name + ".fasta");
		visualizer.writeGraphSimple(gp, etcDir + name + ".dot", name);
	}
}

double get_weight_threshold(PairedInfoLibraries& lib){
	return lib[0]->is_mate_pair_ ?
	        cfg::get().pe_params.param_set.mate_pair_options.select_options.weight_threshold :
	        cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold;
}

double get_single_threshold(PairedInfoLibraries& lib){
    return lib[0]->is_mate_pair_ ?
            *cfg::get().pe_params.param_set.mate_pair_options.select_options.single_threshold :
            *cfg::get().pe_params.param_set.extension_options.select_options.single_threshold;
}

string make_new_name(const std::string& contigs_name,
		const std::string& subname){
	return contigs_name.substr(0, contigs_name.rfind(".fasta")) + "_" + subname + ".fasta";
}


void output_broken_scaffolds(PathContainer& paths, int k, ContigWriter& writer, const std::string& filename) {
    if (!cfg::get().pe_params.param_set.scaffolder_options.on or
            !cfg::get().use_scaffolder or
            cfg::get().pe_params.obs == obs_none) {
        return;
    }

    int min_gap = cfg::get().pe_params.obs == obs_break_all ? k / 2 : k;

    ScaffoldBreaker breaker(min_gap);
    breaker.Split(paths);
    breaker.container().SortByLength();
    writer.writePaths(breaker.container(), filename);
}

vector<SimplePathExtender *> make_extenders(conj_graph_pack& gp,
		const pe_config::ParamSetT& pset, vector<PairedInfoLibraries>& libs,
		bool investigateShortLoops) {
	vector<WeightCounter*> wcs;
	for (size_t i = 0; i < libs.size(); ++i) {
		wcs.push_back(new PathCoverWeightCounter(gp.g, libs[i],
				get_weight_threshold(libs[i]), get_single_threshold(libs[i])));
	}
	vector<SimplePathExtender *> usualPEs;
	for (size_t i = 0; i < libs.size(); ++i) {
		wcs[i]->setNormalizeWeight(pset.normalize_weight);
		wcs[i]->setNormalizeWightByCoverage(pset.normalize_by_coverage);
		double
				priory_coef =
						libs[i][0]->is_mate_pair_ ? pset.mate_pair_options.select_options.priority_coeff
								: pset.extension_options.select_options.priority_coeff;
		SimpleExtensionChooser * extensionChooser = new SimpleExtensionChooser(
				gp.g, wcs[i], priory_coef);
		usualPEs.push_back(new SimplePathExtender(gp.g,
				pset.loop_removal.max_loops, extensionChooser,
				investigateShortLoops));
	}
	return usualPEs;
}

void resolve_repeats_pe_many_libs(conj_graph_pack& gp,
		vector<PairedInfoLibraries>& libs,
		vector<PairedInfoLibraries>& scafolding_libs,
		PathContainer& true_paths,
		size_t RL_true_paths,
		const std::string& output_dir,
		const std::string& contigs_name,
		bool traversLoops,
		boost::optional<std::string> broken_contigs) {

	INFO("Path extend repeat resolving tool started");
	make_dir(output_dir);
	if (cfg::get().developer_mode) {
	    make_dir(get_etc_dir(output_dir));
	}
	const pe_config::ParamSetT& pset = cfg::get().pe_params.param_set;

	INFO("Using " << libs.size() << " paired lib(s)");
	INFO("Using " << scafolding_libs.size() << " scaffolding libs");
	INFO("Scaffolder is " << (pset.scaffolder_options.on ? "on" : "off"));

	ContigWriter writer(gp.g);
	debug_output_edges(writer, gp, output_dir, "before_resolve");

	PathContainer supportingContigs;
	if (FileExists(cfg::get().pe_params.additional_contigs)) {
		INFO("Reading additional contigs " << cfg::get().pe_params.additional_contigs);
		supportingContigs = CreatePathsFromContigs(gp, cfg::get().pe_params.additional_contigs);
		writer.writePaths(supportingContigs, get_etc_dir(output_dir) + "supporting_paths.fasta");
	}

    vector<WeightCounter*> wcs;
    vector<WeightCounter*> scaf_wcs;
    for (size_t i = 0; i < libs.size(); ++i){
        wcs.push_back(new PathCoverWeightCounter(gp.g, libs[i], get_weight_threshold(libs[i]), get_single_threshold(libs[i])));
    }

	vector<SimplePathExtender *> usualPEs = make_extenders(gp, pset, libs, false);

	GapJoiner * gapJoiner = new HammingGapJoiner(gp.g,
				pset.scaffolder_options.min_gap_score,
				(int) (pset.scaffolder_options.max_must_overlap * gp.g.k()),
				(int) (pset.scaffolder_options.max_can_overlap * gp.g.k()),
				pset.scaffolder_options.short_overlap);

	vector<ScaffoldingPathExtender*> scafPEs;
	for (size_t i = 0; i < scafolding_libs.size(); ++i){
		scaf_wcs.push_back(new ReadCountWeightCounter(gp.g, scafolding_libs[i]));
		ScaffoldingExtensionChooser * scafExtensionChooser = new ScaffoldingExtensionChooser(gp.g, scaf_wcs[i],
							scafolding_libs[i][0]->is_mate_pair_? pset.mate_pair_options.select_options.priority_coeff : pset.extension_options.select_options.priority_coeff);
		scafPEs.push_back(new ScaffoldingPathExtender(gp.g, pset.loop_removal.max_loops, scafExtensionChooser, gapJoiner));
	}

	PathExtendResolver resolver(gp.g);
	auto seeds = resolver.makeSimpleSeeds();

	ExtensionChooser * longReadEC = new LongReadsExtensionChooser(gp.g, true_paths, RL_true_paths);
	INFO("Long Reads supporting contigs " << true_paths.size());
	SimplePathExtender * longReadPathExtender = new SimplePathExtender(gp.g, pset.loop_removal.max_loops, longReadEC, false);
	ExtensionChooser * pdEC = new LongReadsExtensionChooser(gp.g, supportingContigs, RL_true_paths);
	SimplePathExtender * pdPE = new SimplePathExtender(gp.g, pset.loop_removal.max_loops, pdEC, false);
	vector<SimplePathExtender *> shortLoopPEs = make_extenders(gp, pset, libs, true);

	vector<PathExtender *> all_libs(usualPEs.begin(), usualPEs.end());
	all_libs.push_back(pdPE);
	all_libs.push_back(longReadPathExtender);
	all_libs.insert(all_libs.end(), shortLoopPEs.begin(), shortLoopPEs.end());
	all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
	CoveringPathExtender * mainPE = new CompositePathExtender(gp.g, pset.loop_removal.max_loops, all_libs);

	seeds.SortByLength();
	INFO("Extending seeds");
	auto paths = resolver.extendSeeds(seeds, *mainPE);
	if (cfg::get().pe_params.output.write_overlaped_paths) {
		writer.writePaths(paths, get_etc_dir(output_dir) + "overlaped_paths.fasta");
	}

	debug_output_paths(writer, gp, output_dir, paths, "overlaped_paths");
    size_t max_over = find_max_overlaped_len(libs);
    resolver.removeOverlaps(paths, mainPE->GetCoverageMap(), max_over, writer, output_dir);
    paths.FilterEmptyPaths();
	paths.CheckSymmetry();
    resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());
	paths.SortByLength();
    debug_output_paths(writer, gp, output_dir, paths, "final_paths");

	if (broken_contigs.is_initialized()) {
	    output_broken_scaffolds(paths, gp.g.k(), writer, output_dir + broken_contigs.get());
	}

    if (traversLoops) {
        INFO("Traversing tandem repeats");
        LoopTraverser loopTraverser(gp.g, paths, mainPE->GetCoverageMap(), mainPE);
        loopTraverser.TraverseAllLoops();
        paths.SortByLength();
    }

    writer.writePaths(paths, output_dir + contigs_name);
    INFO("Reads");
    true_paths.print();
    debug_output_paths(writer, gp, output_dir, true_paths, "long_reads");

    INFO("Path extend repeat resolving tool finished");
}

void add_not_empty_lib(vector<PairedInfoLibraries>& libs, const PairedInfoLibraries& lib){
	if (lib.size() > 0){
		libs.push_back(lib);
	}
}

PairedInfoLibrary* add_lib(conj_graph_pack::graph_t& g,
		vector<PairedIndexT*>& paired_index, vector<size_t>& indexs,
		size_t index, PairedInfoLibraries& libs) {

	size_t read_length = cfg::get().ds.reads[indexs[index]].data().read_length;
	double is = cfg::get().ds.reads[indexs[index]].data().mean_insert_size;
	double var = cfg::get().ds.reads[indexs[index]].data().insert_size_deviation;
	bool is_mp = cfg::get().ds.reads[indexs[index]].type() == io::LibraryType::MatePairs;
	PairedInfoLibrary* lib = new PairedInfoLibrary(cfg::get().K, g, read_length,
			is, var, *paired_index[index], is_mp);
	libs.push_back(lib);
	return lib;
}

void delete_libs(PairedInfoLibraries& libs){
	for (size_t i = 0; i < libs.size(); ++i){
		delete libs[i];
	}
}

void set_threshold(PairedInfoLibrary* lib, size_t index, size_t split_edge_length) {
	INFO("Searching for paired info threshold for lib #"
						<< index << " (IS = " << lib->insert_size_ << ",  DEV = " << lib->is_variation_ << ")");

	SingleThresholdFinder finder(lib->insert_size_ - 2 * lib->is_variation_, lib->insert_size_ + 2 * lib->is_variation_, lib->read_size_);
	double threshold = finder.find_threshold(index);

	INFO("Paired info threshold is " << threshold);
	lib->SetSingleThreshold(threshold);
}

void find_new_threshold(conj_graph_pack& gp, PairedInfoLibrary* lib, size_t index, size_t split_edge_length){
	//SplitGraphPairInfo splitGraph(gp, *lib, index, 99);
	//INFO("Calculating paired info threshold");
	//splitGraph.ProcessReadPairs();
	//double threshold = splitGraph.FindThreshold(split_edge_length, lib->insert_size_ - 2 * lib->is_variation_, lib->insert_size_ + 2 * lib->is_variation_);
	//lib->SetSingleThreshold(threshold);
	double tr = *cfg::get().pe_params.param_set.extension_options.select_options.single_threshold;
	INFO("threshold taken from config - "  << tr);
	lib->SetSingleThreshold(tr);

}

void add_paths_to_container(conj_graph_pack& gp, const std::vector<PathInfo<Graph> >& paths, PathContainer& supportingContigs,
		size_t size_threshold, size_t weight_threshold){

	for (size_t i = 0; i < paths.size(); ++i){
		PathInfo<Graph> path = paths[i];
		if (path.getPath().size() <= size_threshold or path.getWeight() <= weight_threshold){
			continue;
		}
		vector<EdgeId> edges = path.getPath();
		BidirectionalPath* new_path = new BidirectionalPath(gp.g, edges);
		BidirectionalPath* conj_path = new BidirectionalPath(new_path->Conjugate(gp.g.int_id(gp.g.conjugate(edges.back()))));
		new_path->setWeight(path.getWeight());
		conj_path->setWeight(path.getWeight());
		supportingContigs.AddPair(new_path, conj_path);
	}
}

size_t get_max_read_length(){
	size_t RL = 0;
	for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i){
		RL = std::max(RL, cfg::get().ds.reads[i].data().read_length);
	}
	return RL;
}

void resolve_repeats_pe(conj_graph_pack& gp,
		vector<PairedIndexT*>& paired_index, vector<PairedIndexT*>& scaff_index,
		vector<size_t>& indexs, const std::vector<PathInfo<Graph> >& true_paths,
		const std::string& output_dir, const std::string& contigs_name, bool traverseLoops, boost::optional<std::string> broken_contigs, bool use_auto_threshold = true) {

    const pe_config::ParamSetT& pset = cfg::get().pe_params.param_set;

	PairedInfoLibraries paired_end_libs;
	PairedInfoLibraries mate_pair_libs;
	PairedInfoLibraries pe_scaf_libs;
	PairedInfoLibraries mp_scaf_libs;

	for (size_t i = 0; i < paired_index.size(); ++i) {
		if (cfg::get().ds.reads[indexs[i]].type()
				== io::LibraryType::PairedEnd) {
			PairedInfoLibrary* lib = add_lib(gp.g, paired_index, indexs, i, paired_end_libs);

			if (use_auto_threshold){
				find_new_threshold(gp, lib, indexs[i], pset.split_edge_length);
			}
		}
		else if (cfg::get().ds.reads[indexs[i]].type()
				== io::LibraryType::MatePairs) {
			PairedInfoLibrary* lib = add_lib(gp.g, paired_index, indexs, i, mate_pair_libs);
			if (use_auto_threshold){
				find_new_threshold(gp, lib, indexs[i], pset.split_edge_length);
			}
		}
	}

    for (size_t i = 0; i < scaff_index.size(); ++i) {
        if (cfg::get().ds.reads[indexs[i]].type()
                == io::LibraryType::PairedEnd) {
            add_lib(gp.g, scaff_index, indexs, i, pe_scaf_libs);
        }
        else if (cfg::get().ds.reads[indexs[i]].type()
                == io::LibraryType::MatePairs) {
            add_lib(gp.g, scaff_index, indexs, i, mp_scaf_libs);
        }
    }

	vector<PairedInfoLibraries> rr_libs;
	add_not_empty_lib(rr_libs, paired_end_libs);
	add_not_empty_lib(rr_libs, mate_pair_libs);
	vector<PairedInfoLibraries> scaff_libs;
	add_not_empty_lib(scaff_libs, pe_scaf_libs);
	add_not_empty_lib(scaff_libs, mp_scaf_libs);
	PathContainer supportingContigs;
	add_paths_to_container(gp, true_paths, supportingContigs, 1, 1);

	INFO("== Long reads paths == ");
	for (auto it = true_paths.begin(); it != true_paths.end(); ++it) {
	    INFO("Path " << gp.g.int_id(it->path[0]));
	    for (auto e = it->path.begin(); e != it->path.end(); ++e) {
	        INFO(gp.g.int_id(*e) << ", length " << gp.g.length(*e));
	    }
	}
	INFO("==== ");



	resolve_repeats_pe_many_libs(gp, rr_libs, scaff_libs, supportingContigs, get_max_read_length(), output_dir, contigs_name, traverseLoops, broken_contigs);

	delete_libs(paired_end_libs);
	delete_libs(mate_pair_libs);
	delete_libs(pe_scaf_libs);
	delete_libs(mp_scaf_libs);
}


} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
