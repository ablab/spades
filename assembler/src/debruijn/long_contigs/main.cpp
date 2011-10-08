/*
 * main.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: andrey
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../config_struct.hpp"
#include "simple_tools.hpp"
#include "lc_config_struct.hpp"

#include "lc_common.hpp"
#include "lc_io.hpp"
#include "seeds.hpp"
#include "path_utils.hpp"
#include "paths.hpp"
#include "quality.hpp"
#include "visualize.hpp"

DECL_PROJECT_LOGGER("d")

int main() {
	cfg::create_instance(debruijn_graph::cfg_filename);
	using namespace long_contigs;
	lc_cfg::create_instance(lc_cfg_filename);

	checkFileExistenceFATAL(lc_cfg_filename);
	checkFileExistenceFATAL(debruijn_graph::cfg_filename);

	Graph g(K);
	EdgeIndex<K + 1, Graph> index(g);
	IdTrackHandler<Graph> intIds(g);
	PairedInfoIndex<Graph> pairedIndex(g, 0);
	PairedInfoIndices pairedInfos;
	KmerMapper<K+1, Graph> mapper(g);
	Sequence sequence("");

	std::vector<BidirectionalPath> seeds;
	std::vector<BidirectionalPath> rawSeeds;
	std::vector<BidirectionalPath> filteredSeeds;
	std::vector<BidirectionalPath> paths;

	string output_dir = cfg::get().output_dir;

	PathStopHandler stopHandler(g);

	if (!lc_cfg::get().from_file) {
		INFO("Load from file");
		return -1;
	}
	else {
		LoadFromFile(lc_cfg::get().ds.graph_file, &g, &intIds, sequence, &mapper);
	}
	mkdir(cfg::get().output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

	Path<Graph::EdgeId> path1 = FindGenomePath<K> (sequence, g, index);
	Path<Graph::EdgeId> path2 = FindGenomePath<K> (!sequence, g, index);

	if (cfg::get().etalon_info_mode) {
		AddEtalonInfo<K>(g, index, sequence, pairedInfos);
	} else {
		pairedInfos.clear();
		AddRealInfo<K>(g, index, intIds, pairedInfos, mapper, lc_cfg::get().use_new_metrics);

		if (!cfg::get().etalon_info_mode && lc_cfg::get().write_real_paired_info) {
			SavePairedInfo(g, pairedInfos, intIds, output_dir + lc_cfg::get().paired_info_file_prefix);
		}
	}

	FindSeeds(g, rawSeeds, &pairedInfos);
	INFO("Seeds found");

	ResolveUnequalComplement(g, rawSeeds, lc_cfg::get().ss.sym.cut_tips, lc_cfg::get().ss.sym.min_conjugate_len);

	if (lc_cfg::get().fo.remove_sefl_conjugate) {
		RemoveWrongConjugatePaths(g, rawSeeds, &filteredSeeds);
	} else {
		filteredSeeds = rawSeeds;
	}

	//RemoveUnagreedPaths();

	std::vector<int> seedPairs;
	std::vector<double> seedQuality;
	FilterComplement(g, filteredSeeds, &seedPairs, &seedQuality);

	RemoveSubpaths(g, filteredSeeds, seeds);
	INFO("Sub seeds removed");

	FilterComplement(g, seeds, &seedPairs, &seedQuality);

	if (lc_cfg::get().rs.research_mode && lc_cfg::get().rs.fiter_by_edge) {
		FilterEdge(g, seeds, lc_cfg::get().rs.edge_length);
	}

	double MIN_COVERAGE = lc_cfg::get().ss.min_coverage;
	FilterLowCovered(g, seeds, MIN_COVERAGE);
	INFO("Seeds filtered");

	FilterComplement(g, seeds, &seedPairs, &seedQuality);

	size_t found = PathsInGenome<K>(g, index, sequence, seeds, path1, path2);
	INFO("Good seeds found " << found << " in total " << seeds.size());
	INFO("Seed coverage " << PathsCoverage(g, seeds));
	INFO("Path length coverage " << PathsLengthCoverage(g, seeds));

	if (lc_cfg::get().write_seeds) {
		WriteGraphWithPathsSimple(output_dir + "seeds.dot", "seeds", g, seeds, path1, path2);

		OutputPathsAsContigsNoComplement(g, seeds, seedPairs, output_dir + "seeds.contigs", std::set<int>());
	}

	if (lc_cfg::get().total_symmetric_mode) {
		FindPathsSymmetric(g, seeds, pairedInfos,  stopHandler, seedPairs);
		paths.resize(seeds.size());
		std::copy(seeds.begin(), seeds.end(), paths.begin());
	} else {
		FindPaths(g, seeds, pairedInfos, paths, stopHandler);
	}

	ResolveUnequalComplement(g, paths, lc_cfg::get().fo.sym.cut_tips, lc_cfg::get().fo.sym.min_conjugate_len);

	std::vector<int> pairs;
	std::vector<double> quality;
	FilterComplement(g, filteredSeeds, &pairs, &quality);

	for (int i = 0; i < (int) pairs.size(); ++i) {
		DETAILED_INFO("Paired result: " << i << " - " << pairs[i] << " : " << quality[i]);
	}

	std::vector<BidirectionalPath> filteredPaths;;
	if (lc_cfg::get().fo.remove_sefl_conjugate) {
		RemoveWrongConjugatePaths(g, paths, &filteredPaths);
	} else {
		filteredPaths = paths;
	}

	std::vector<BidirectionalPath> result;
	std::vector<double> pathQuality;
	if (lc_cfg::get().fo.remove_subpaths || lc_cfg::get().fo.remove_overlaps) {
		RemoveSubpaths(g, filteredPaths, result, &pathQuality);
		INFO("Subpaths removed");
	}
	else if (lc_cfg::get().fo.remove_duplicates) {
		RemoveDuplicate(g, filteredPaths, result, &pathQuality);
		INFO("Duplicates removed");
	}
	else {
		result = filteredPaths;
		SortPathsByLength(g, result);
		pathQuality.resize(result.size(), 1.0);
	}

	if (lc_cfg::get().write_overlaped_paths) {
		WriteGraphWithPathsSimple(output_dir + "overlaped_paths.dot", "overlaped_paths", g, result, path1, path2);
	}

	found = PathsInGenome<K>(g, index, sequence, result, path1, path2, &pathQuality);
	INFO("Good paths found " << found << " in total " << result.size());
	INFO("Path coverage " << PathsCoverage(g, result));
	INFO("Path length coverage " << PathsLengthCoverage(g, result));

	if (lc_cfg::get().print_stats) {
		stopHandler.print();
	}

	if (lc_cfg::get().write_paths) {
		WriteGraphWithPathsSimple(output_dir + "final_paths.dot", "final_paths", g, result, path1, path2);
	}

	if (lc_cfg::get().write_contigs) {
		OutputPathsAsContigs(g, result, output_dir + "all_paths.contigs");
		OutputContigsNoComplement(g, output_dir + "edges.contigs");

		FilterComplement(g, result, &pairs, &quality);

		std::set<int> toRemove;
		if (lc_cfg::get().fo.remove_overlaps) {
			RemoveOverlaps(g, result);
			DETAILED_INFO("Removed overlaps");

			if (lc_cfg::get().fo.remove_similar) {
				RemoveSimilar(g, result, pathQuality, toRemove);
			}
			DETAILED_INFO("Removed similar");

		}
		OutputPathsAsContigsNoComplement(g, result, pairs, output_dir + "paths.contigs", toRemove);
		INFO("All contigs written");
	}

	if (lc_cfg::get().write_graph) {
		SaveGraph(g, intIds, output_dir + "graph");
	}

	PrintPath(g, path1);
	PrintPath(g, path2);

	INFO("Tool finished");
	DeleteAdditionalInfo(pairedInfos);
	return 0;
}

