//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * lc_launch.hpp
 *
 *  Created on: Dec 1, 2011
 *      Author: andrey
 */

#ifndef LC_LAUNCH_HPP_
#define LC_LAUNCH_HPP_

#include "lc_config_struct.hpp"

#include "lc_common.hpp"
#include "lc_io.hpp"
#include "seeds.hpp"
#include "path_utils.hpp"
#include "paths.hpp"
#include "quality.hpp"
#include "visualize.hpp"
#include <boost/optional.hpp>

namespace long_contigs {

using namespace debruijn_graph;

void resolve_repeats_ml(const conj_graph_pack& gp, PairedInfoIndices& pairedInfos, std::string output_dir, const lc_config::lc_params& p,
		boost::optional<const PairedInfoIndex<Graph>&> jump_index_opt = boost::none) {
	INFO("Multilayer resolving tool started");

//	for (auto it  = pairedInfos.begin(); it != pairedInfos.end(); ++it) {
//		cout << "Index size " << it->pairedInfoIndex->size() << endl;
//	}

	make_dir(output_dir);

	params = p;

//	EdgeIndex<K + 1, Graph> index(g);

	std::vector<BidirectionalPath> seeds;
	std::vector<BidirectionalPath> rawSeeds;
	std::vector<BidirectionalPath> filteredSeeds;
	std::vector<BidirectionalPath> lowCoveredSeeds;

	const Graph& g = gp.g;

	PathStopHandler stopHandler(g);

	std::vector<int> seedPairs;
	std::vector<double> seedQuality;

	FindSeeds(g, rawSeeds, &pairedInfos);
	CheckIds(g, rawSeeds);

	ResolveUnequalComplement(g, rawSeeds, params.ps.ss.sym.cut_tips, params.ps.ss.sym.min_conjugate_len);
	CheckIds(g, rawSeeds);

	std::vector<BidirectionalPath> goodSeeds;
	RemoveUnagreedPaths(g, rawSeeds, pairedInfos, params.ps.fo.agreed_coeff, &goodSeeds);
	CheckIds(g, goodSeeds);

	if (params.ps.fo.remove_sefl_conjugate) {
		RemoveWrongConjugatePaths(g, goodSeeds, params.ps.ss.sym.min_conjugate_len, &lowCoveredSeeds);
	} else {
		lowCoveredSeeds = goodSeeds;
	}
	//RemoveUnagreedPaths();
	INFO("Wrong conjugate removed");
	CheckIds(g, lowCoveredSeeds);


	FilterLowCovered(g, lowCoveredSeeds, params.ps.ss.min_coverage, &filteredSeeds);
	INFO("Seeds filtered");
	CheckIds(g, filteredSeeds);

	RemoveSubpaths(g, filteredSeeds, seeds);
	INFO("Sub seeds removed");
	CheckIds(g, seeds);

//	if (params.rs.research_mode && params.rs.fiter_by_edge) {
//		FilterEdge(g, seeds, params.rs.edge_length);
//	}


	FilterComplement(g, seeds, &seedPairs, &seedQuality);

	//PrintPathsShort(g, seeds);
//	size_t found = PathsInGenome<K>(g, index, sequence, seeds, path1, path2);
//	INFO("Good seeds found " << found << " in total " << seeds.size());
//	INFO("Seed coverage " << PathsCoverage(g, seeds));
//	INFO("Path length coverage " << PathsLengthCoverage(g, seeds));
//
	if (params.write_seeds) {
        INFO("=== Final seeds ===");
        for (size_t i = 0; i < seeds.size(); ++i) {
            PrintPath(g, seeds[i]);
        }

		WriteGraphWithPathsSimple(gp, output_dir + "seeds.dot", "seeds", seeds);
		OutputPathsAsContigsNoComplement(g, seeds, output_dir + "seeds.fasta", std::set<int>());
	}

//	if (params.total_symmetric_mode) {
//		FindPathsSymmetric(g, seeds, pairedInfos,  stopHandler, seedPairs);
//		paths.resize(seeds.size());
//		std::copy(seeds.begin(), seeds.end(), paths.begin());
//	} else {
	PairedInfoIndex<Graph> empty_info(g);
	if (!jump_index_opt) {
//		cout << "FAIL!!!" << endl;
//		jump_index_opt = in_place<const PairedInfoIndex<Graph>&>(cref(g));
		jump_index_opt.reset(empty_info);
	}
	FindPaths(g, seeds, pairedInfos, stopHandler, *jump_index_opt);
//	}
	std::vector<BidirectionalPath> & paths = seeds;
	CheckIds(g, paths);

	ResolveUnequalComplement(g, paths, params.ps.fo.sym.cut_tips, params.ps.fo.sym.min_conjugate_len);
	CheckIds(g, paths);
	PrintPathsShort(g, paths);
	std::vector<int> pairs;
	std::vector<double> quality;
	INFO("After unequal")
	FilterComplement(g, paths, &pairs, &quality);

	std::vector<BidirectionalPath> goodPaths;
	RemoveUnagreedPaths(g, paths, pairedInfos, params.ps.fo.agreed_coeff, &goodPaths);
	CheckIds(g, goodPaths);


	INFO("After unagreed")
	FilterComplement(g, goodPaths, &pairs, &quality);

	std::vector<BidirectionalPath> filteredPaths;;
	if (params.ps.fo.remove_sefl_conjugate) {
		RemoveWrongConjugatePaths(g, goodPaths, params.ps.fo.sym.min_conjugate_len, &filteredPaths);
	} else {
		filteredPaths = goodPaths;
	}

	CheckIds(g, filteredPaths);
	INFO("After wrong conjugate")
	FilterComplement(g, filteredPaths, &pairs, &quality);

	std::vector<BidirectionalPath> result;
	std::vector<double> pathQuality;
	if (params.ps.fo.remove_subpaths || params.ps.fo.remove_overlaps) {
		RemoveSubpaths(g, filteredPaths, result, &pathQuality);
		INFO("Subpaths removed");
	}
	else if (params.ps.fo.remove_duplicates) {
		RemoveDuplicate(g, filteredPaths, result, &pathQuality);
		INFO("Duplicates removed");
	}
	else {
		result = filteredPaths;
		SortPathsByLength(g, result);
		pathQuality.resize(result.size(), 1.0);
	}

	if (params.write_overlaped_paths) {
		WriteGraphWithPathsSimple(gp, output_dir + "overlaped_paths.dot", "overlaped_paths", result);
	}

	CheckIds(g, result);
	INFO("After subpaths")
	FilterComplement(g, result, &pairs, &quality);

	if (params.print_stats) {
		stopHandler.print();
	}

//	if (params.write_paths) {
//		WriteGraphWithPathsSimple(gp, output_dir + "final_paths.dot", "final_paths", result);
//	}

	if (params.write_path_loc) {
		make_dir(output_dir + "path_loc/");
		WritePathLocalities(gp, output_dir + "path_loc/", result);
	}

	if (params.write_contigs) {
		OutputPathsAsContigs(g, result, output_dir + "all_paths.fasta");
		OutputContigsNoComplement(g, output_dir + "complement_filtered.fasta");

		std::set<int> toRemove;
		std::vector<BidirectionalPath> noOverlaps;

		if (params.ps.fo.remove_overlaps) {
			RemoveOverlaps(g, result);
			INFO("Removed overlaps");
			CheckIds(g, result);
			//TODO: remove subseeds again
		}
		if (params.ps.fo.remove_similar) {
			RemoveSimilar(g, result, pathQuality, &noOverlaps);
		} else {
			noOverlaps = result;
		}
		INFO("Removed similar");

		if (params.ps.fo.write_uncovered_edges) {
		    AddUncoveredEdges(g, noOverlaps);
		}

//		found = PathsInGenome<K>(g, index, sequence, noOverlaps, path1, path2, &pathQuality);
//		INFO("Good paths found " << found << " in total " << noOverlaps.size());
//		INFO("Path coverage " << PathsCoverage(g, noOverlaps));
//		INFO("Path length coverage " << PathsLengthCoverage(g, noOverlaps));

		PrintPathsShort(g, noOverlaps);

		OutputPathsAsContigsNoComplement(g, noOverlaps, output_dir + "resolved.fasta", toRemove);
		INFO("All contigs written");

		INFO("=== Final paths ===");
		for (size_t i = 0; i < noOverlaps.size(); ++i) {
		    PrintPath(g, noOverlaps[i]);
		}
	}

//	PrintPath(gp.g, FindGenomePath<K>(gp.genome, gp.g, gp.index));
//	PrintPath(gp.g, FindGenomePath<K>(!gp.genome, gp.g, gp.index));

	INFO("Tool finished");
	DeleteAdditionalInfo(pairedInfos);
}

void resolve_repeats_ml(const conj_graph_pack& gp, PairedInfoIndex<Graph>& paired_index, const std::string& output_dir, const lc_config::lc_params& p,
		boost::optional<const PairedInfoIndex<Graph>&> jump_index_opt = boost::none) {
    PairedInfoIndices pairedInfos;
    pairedInfos.push_back(PairedInfoIndexLibrary(gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, 2, size_t(*cfg::get().ds.is_var), 5, &paired_index));

    resolve_repeats_ml(gp, pairedInfos, output_dir, p, jump_index_opt);
}

} /* long_contigs */



#endif /* LC_LAUNCH_HPP_ */
