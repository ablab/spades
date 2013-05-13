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
#include "path_validator.hpp"
#include "loop_traverser.hpp"
#include "single_threshold_finder.hpp"
#include "long_read_storage.hpp"
#include "../include/omni/edges_position_handler.hpp"

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

string get_etc_dir(const std::string& output_dir){
	return output_dir + "path_extend/";
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
	double pe_weight_threshold = cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold;
	double mp_weight_threshold = cfg::get().pe_params.param_set.mate_pair_options.select_options.weight_threshold;
	return lib[0]->is_mate_pair_ ? mp_weight_threshold :pe_weight_threshold;

}

void resolve_repeats_pe_many_libs(size_t k, conj_graph_pack& gp,
		vector<PairedInfoLibraries>& libes,
		vector<PairedInfoLibraries>& scafolding_libes,
		PathContainer& true_paths,
		const std::string& output_dir, const std::string& contigs_name) {

	INFO("Path extend repeat resolving tool started");
	if (!cfg::get().developer_mode) {
	    cfg::get_writable().pe_params.viz.DisableAll();
	    cfg::get_writable().pe_params.output.DisableAll();
	}
	make_dir(output_dir);
	std::string etcDir = get_etc_dir(output_dir);
	make_dir(etcDir);

	DEBUG("Using " << cfg::get().pe_params.param_set.metric << " metric");
	INFO("Using " << libes.size() << " paird libs");
	INFO("Scaffolder is " << (cfg::get().pe_params.param_set.scaffolder_options.on ? "on" : "off"));
	ContigWriter writer(gp, gp.g.k());
	debug_output_edges(writer, gp, output_dir, "before_resolve");
	vector<WeightCounter*> wcs;
	vector<WeightCounter*> scaf_wcs;
	PathContainer supportingContigs;
	if (FileExists(cfg::get().pe_params.additional_contigs)) {
		INFO("Reading additional contigs " << cfg::get().pe_params.additional_contigs);
		supportingContigs = CreatePathsFromContigs(gp, cfg::get().pe_params.additional_contigs);
		writer.writePaths(supportingContigs, etcDir + "supporting_paths.fasta");
	}
	if (cfg::get().pe_params.param_set.metric == "read_count") {
		for (size_t i = 0; i < libes.size(); ++i){
			wcs.push_back(new ReadCountWeightCounter(gp.g, libes[i], get_weight_threshold(libes[i])));
		}
	} else if (cfg::get().pe_params.param_set.metric == "path_cover") {
		for (size_t i = 0; i < libes.size(); ++i){
			double pe_single_threshold = *cfg::get().pe_params.param_set.extension_options.select_options.single_threshold;
			double mp_single_threshold = *cfg::get().pe_params.param_set.mate_pair_options.select_options.single_threshold;
			double single_threshold = libes[i][0]->is_mate_pair_ ? mp_single_threshold :pe_single_threshold;
			wcs.push_back(new PathCoverWeightCounter(gp.g, libes[i], get_weight_threshold(libes[i]), single_threshold));
		}
	} else {
		WARN("Unknown metric");
		return;
	}
	vector<SimplePathExtender *> usualPEs;
	for (size_t i = 0; i < wcs.size(); ++i){
		wcs[i]->setNormalizeWeight(cfg::get().pe_params.param_set.normalize_weight);
		wcs[i]->setNormalizeWightByCoverage(cfg::get().pe_params.param_set.normalize_by_coverage);
		double mp_pr_coef = cfg::get().pe_params.param_set.mate_pair_options.select_options.priority_coeff;
		double pe_pr_coef = cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff;
		double priory_coef = libes[i][0]->is_mate_pair_ ? mp_pr_coef : pe_pr_coef;
		SimpleExtensionChooser * extensionChooser = new SimpleExtensionChooser(gp.g, wcs[i], priory_coef);
		usualPEs.push_back( new SimplePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, extensionChooser));
	}
	DEBUG("usualPE " << usualPEs.size());
	GapJoiner * gapJoiner = new HammingGapJoiner(gp.g,
				cfg::get().pe_params.param_set.scaffolder_options.min_gap_score,
				(int) (cfg::get().pe_params.param_set.scaffolder_options.max_must_overlap
						* gp.g.k()),
				(int) (cfg::get().pe_params.param_set.scaffolder_options.max_can_overlap
						* gp.g.k()),
				cfg::get().pe_params.param_set.scaffolder_options.short_overlap);
	vector<ScaffoldingPathExtender*> scafPEs;
	for (size_t i = 0; i < scafolding_libes.size(); ++i){
		scaf_wcs.push_back(new ReadCountWeightCounter(gp.g, scafolding_libes[i]));
		ScaffoldingExtensionChooser * scafExtensionChooser =
					new ScaffoldingExtensionChooser(gp.g, scaf_wcs[i],
							scafolding_libes[i][0]->is_mate_pair_?cfg::get().pe_params.param_set.mate_pair_options.select_options.priority_coeff : cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);
		scafPEs.push_back(new ScaffoldingPathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, scafExtensionChooser, gapJoiner));
	}
	PathExtendResolver resolver(gp.g, k);
	auto seeds = resolver.makeSimpleSeeds();
	ExtensionChooser * longReadEC = new LongReadsExtensionChooser(gp.g, true_paths);
	INFO("Long Reads supporting contigs " << true_paths.size());
	SimplePathExtender * longReadPathExtender = new SimplePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, longReadEC);
	ExtensionChooser * pdEC = new LongReadsExtensionChooser(gp.g, supportingContigs);
	SimplePathExtender * pdPE = new SimplePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, pdEC);
	vector<PathExtender *> all_libs(usualPEs.begin(), usualPEs.end());
	all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
	all_libs.push_back(pdPE);
	all_libs.push_back(longReadPathExtender);
	CoveringPathExtender * mainPE = new CompositePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, all_libs);

	seeds.SortByLength();
	INFO("Extending seeds");
	auto paths = resolver.extendSeeds(seeds, *mainPE);
	if (cfg::get().pe_params.output.write_overlaped_paths) {
		writer.writePaths(paths, etcDir + "overlaped_paths.fasta");
	}
	debug_output_paths(writer, gp, output_dir, paths, "overlaped_paths");
    size_t max_over = find_max_overlaped_len(libes);
    resolver.removeOverlaps(paths, mainPE->GetCoverageMap(), max_over);
	paths.CheckSymmetry();
    resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());
	paths.SortByLength();
//    for (size_t i = 0; i < paths.size(); ++i){
//        paths.Get(i)->Print();
//    }
    paths.SortByLength();
    writer.writePaths(paths, output_dir + contigs_name);
    debug_output_paths(writer, gp, output_dir, paths, "final_paths");
	INFO("Traversing tandem repeats");
    LoopTraverser loopTraverser(gp.g, paths, mainPE->GetCoverageMap(), mainPE);
	loopTraverser.TraverseAllLoops();
	paths.SortByLength();
	INFO("Found " << paths.size() << " contigs");
	writer.writePaths(paths, output_dir + contigs_name.substr(0, contigs_name.rfind(".fasta")) + "_loop_tr.fasta");
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

void set_threshold(PairedInfoLibrary* lib, size_t index) {
	INFO("Searching for paired info threshold for lib #"
						<< index << " (IS = " << lib->insert_size_ << ",  DEV = " << lib->is_variation_);

	SingleThresholdFinder finder(lib->insert_size_ - 2 * lib->is_variation_, lib->insert_size_ + 2 * lib->is_variation_, 99);
	double threshold = finder.find_threshold(index);

	INFO("Paired info threshold is " << threshold);
	lib->SetSingleThreshold(threshold);
}

void add_paths_to_container(conj_graph_pack& gp, const std::vector<LongReadInfo<Graph> >& paths, PathContainer& supportingContigs){
	for (size_t i = 0; i < paths.size(); ++i){
		LongReadInfo<Graph> path = paths[i];
		vector<EdgeId> edges = path.getPath();
		BidirectionalPath* new_path = new BidirectionalPath(gp.g, edges);
		BidirectionalPath* conj_path = new BidirectionalPath(new_path->Conjugate());
		new_path->setWeight(path.getWeight());
		conj_path->setWeight(path.getWeight());
		supportingContigs.AddPair(new_path, conj_path);
	}
}

void resolve_repeats_pe(size_t k, conj_graph_pack& gp,
		vector<PairedIndexT*>& paired_index, vector<PairedIndexT*>& scaff_index,
		vector<size_t>& indexs, const std::vector<LongReadInfo<Graph> >& true_paths,
		const std::string& output_dir, const std::string& contigs_name) {

	PairedInfoLibraries libs;
	PairedInfoLibraries libs_mp;
	PairedInfoLibraries scaf_libs;
	PairedInfoLibraries scaf_libs_mp;

	for (size_t i = 0; i < paired_index.size(); ++i) {
		if (cfg::get().ds.reads[indexs[i]].type()
				== io::LibraryType::PairedEnd) {
			PairedInfoLibrary* lib = add_lib(gp.g, paired_index, indexs, i, libs);
			//set_threshold(lib, indexs[i]);
			add_lib(gp.g, scaff_index, indexs, i, scaf_libs);
		} else if (cfg::get().ds.reads[indexs[i]].type()
				== io::LibraryType::MatePairs) {
			PairedInfoLibrary* lib = add_lib(gp.g, paired_index, indexs, i, libs_mp);
			//set_threshold(lib, indexs[i]);
			add_lib(gp.g, scaff_index, indexs, i, scaf_libs_mp);
		}
	}
	vector<PairedInfoLibraries> libes;
	add_not_empty_lib(libes, libs);
	add_not_empty_lib(libes, libs_mp);
	vector<PairedInfoLibraries> scaff_libes;
	add_not_empty_lib(scaff_libes, scaf_libs);
	add_not_empty_lib(scaff_libes, scaf_libs_mp);
	PathContainer supportingContigs;
	add_paths_to_container(gp, true_paths, supportingContigs);
	resolve_repeats_pe_many_libs(k, gp, libes, scaff_libes, supportingContigs,
			output_dir, contigs_name);
	delete_libs(libs);
	delete_libs(libs_mp);
	delete_libs(scaf_libs);
	delete_libs(scaf_libs_mp);
}


} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
