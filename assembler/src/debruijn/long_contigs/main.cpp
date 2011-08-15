/*
 * main.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: andrey
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "lc_common.hpp"
#include "lc_io.hpp"
#include "seeds.hpp"
#include "paths.hpp"
#include "quality.hpp"
#include "visualize.hpp"

#include "simple_tools.hpp"

namespace {

std::string MakeLaunchTimeDirName() {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%m.%d_%H_%M", timeinfo);
	return string(buffer);
}
}

DECL_PROJECT_LOGGER("d")

int main() {
	using namespace long_contigs;

	checkFileExistenceFATAL(LC_CONFIG_FILENAME);
	checkFileExistenceFATAL(CONFIG_FILENAME);

	Graph g(K);
	EdgeIndex<K + 1, Graph> index(g);
	IdTrackHandler<Graph> intIds(g);
	PairedInfoIndex<Graph> pairedIndex(g, 0);
	PairedInfoIndices pairedInfos;
	Sequence sequence("");

	std::vector<BidirectionalPath> seeds;
	std::vector<BidirectionalPath> rawSeeds;
	std::vector<BidirectionalPath> paths;

	std::string dataset = CONFIG.read<string>("dataset");
	string output_root = CONFIG.read<string>("output_dir");
	string output_dir_suffix = MakeLaunchTimeDirName()+ "." + dataset + "/";
	string output_dir = output_root + output_dir_suffix;

	if (!LC_CONFIG.read<bool>("from_file")) {
		INFO("Load from file");
		return -1;
	}
	else {
		LoadFromFile<K>(LC_CONFIG.read<std::string>("graph_file_" + dataset), g, pairedIndex, index, intIds, sequence);
	}
	mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

	Path<Graph::EdgeId> path1 = FindGenomePath<K> (sequence, g, index);
	Path<Graph::EdgeId> path2 = FindGenomePath<K> (!sequence, g, index);

	PairedInfoIndexLibrary basicPairedLib(LC_CONFIG.read<size_t>("read_size"), LC_CONFIG.read<size_t>("insert_size"), &pairedIndex);
	pairedInfos.push_back(basicPairedLib);

	if (CONFIG.read<bool>("etalon_info_mode")) {
		AddEtalonInfo<K>(g, index, sequence, pairedInfos);
	} else {
		AddRealInfo<K>(g, index, intIds, pairedInfos);
	}

	FindSeeds(g, rawSeeds);
	INFO("Seeds found");
	RemoveSubpaths(g, rawSeeds, seeds);
	INFO("Sub seeds removed");

	double MIN_COVERAGE = LC_CONFIG.read<double>("min_coverage");
	FilterLowCovered(g, seeds, MIN_COVERAGE);
	INFO("Seeds filtered");

	size_t found = PathsInGenome<K>(g, index, sequence, seeds, path1, path2, true);
	INFO("Good seeds found " << found << " in total " << seeds.size());
	INFO("Seed coverage " << PathsCoverage(g, seeds));
	INFO("Path length coverage " << PathsLengthCoverage(g, seeds));

	if (LC_CONFIG.read<bool>("write_seeds")) {
		WriteGraphWithPathsSimple(output_dir + "seeds.dot", "seeds", g, seeds, path1, path2);
	}

	FindPaths(g, seeds, pairedInfos, paths);

	std::vector<BidirectionalPath> result;
	if (LC_CONFIG.read<bool>("remove_duplicates_only") || LC_CONFIG.read<bool>("remove_overlaps")) {
		RemoveSubpaths(g, paths, result);
		INFO("Duplicates removed");
	}
	else {
		RemoveSubpaths(g, paths, result);
		INFO("Subpaths removed");
	}

	found = PathsInGenome<K>(g, index, sequence, result, path1, path2);
	INFO("Good paths found " << found << " in total " << result.size());
	INFO("Path coverage " << PathsCoverage(g, result));
	INFO("Path length coverage " << PathsLengthCoverage(g, result));

	if (LC_CONFIG.read<bool>("write_overlaped_paths")) {
		WriteGraphWithPathsSimple(output_dir + "overlaped_paths.dot", "overlaped_paths", g, result, path1, path2);
	}

	if (LC_CONFIG.read<bool>("remove_overlaps")) {
		RemoveOverlaps(result);
	}

	if (LC_CONFIG.read<bool>("write_paths")) {
		WriteGraphWithPathsSimple(output_dir + "final_paths.dot", "final_paths", g, result, path1, path2);
	}

	if (LC_CONFIG.read<bool>("write_contigs")) {
		OutputPathsAsContigs(g, result, output_dir + "paths.contigs");
	}

	if (!CONFIG.read<bool>("etalon_info_mode") && LC_CONFIG.read<bool>("write_real_paired_info")) {
		SavePairedInfo(g, pairedInfos, intIds, output_dir + LC_CONFIG.read<std::string>("paired_info_file_prefix"));
	}

	DeleteAdditionalInfo(pairedInfos);
	return 0;
}

