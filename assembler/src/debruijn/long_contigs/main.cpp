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

namespace long_contigs {

template<size_t k, class ReadStream>
void BuildDeBruijnGraph(ReadStream& stream,
		const Sequence& genome, bool paired_mode, bool rectangle_mode,
		bool etalon_info_mode, bool from_saved, size_t insert_size,
		size_t max_read_length, const string& output_folder,
		const string& work_tmp_dir,
		Graph& g, PairedInfoIndex<Graph>& paired_index, EdgeIndex<k + 1, Graph>& index) {

	INFO("Edge graph construction tool started");
	INFO("Paired mode: " << (paired_mode ? "Yes" : "No") );
	INFO("Etalon paired info mode: " << (etalon_info_mode ? "Yes" : "No"))
	INFO("From file: " << (from_saved ? "Yes" : "No"))
	mkdir(work_tmp_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

	IdTrackHandler<Graph> IntIds(g);
	EdgesPositionHandler<Graph> EdgePos(g);
	EdgesPosGraphLabeler<Graph> EdgePosLab(g, EdgePos);

	if (!from_saved) {

		if (paired_mode) {
			if (etalon_info_mode) {
				ConstructGraphWithEtalonPairedInfo<k, ReadStream> (g, index,
						paired_index, stream, insert_size, max_read_length,
						genome);
			} else {
				ConstructGraphWithPairedInfo<k, ReadStream> (g, index,
						paired_index, stream);
			}
		} else {
			typedef io::ConvertingReaderWrapper UnitedStream;
			UnitedStream united_stream(&stream);
			ConstructGraphWithCoverage<k, UnitedStream> (g, index,
					united_stream);
		}

		ProduceInfo<k> (g, index, genome, output_folder + "edge_graph.dot",
				"edge_graph");


		FillEdgesPos<k>(g, index, genome, EdgePos);
		omnigraph::WriteSimple(
				output_folder + "before_simplification_pos.dot",
				"no_repeat_graph", g, EdgePosLab);


		SimplifyGraph<k> (g, index, 3, genome, output_folder);
//		MapPairedReads<k, ReadStream, Graph>(g, stream, index);

		ProduceInfo<k> (g, index, genome,
				output_folder + "simplified_graph.dot", "simplified_graph");

		WriteGraphComponents<k> (g, index, genome,
				output_folder + "graph_components" + "/", "graph.dot",
				"graph_component", insert_size);
//		if (paired_mode) {
//			ProducePairedInfoStats(g, insert_size, max_read_length, paired_index, output_folder);
//		}
	}
	//	if (paired_mode) {
	//		paired_index.OutputData(output_folder + "edges_dist.txt");
	//
	//		SimpleOfflineClusterer<Graph> clusterer(paired_index);
	//		PairedInfoIndex<Graph> clustered_paired_index(g);
	//		clusterer.cluster(clustered_paired_index);
	//	}

	omnigraph::WriteSimple(
			output_folder + "repeats_resolved_before_poslab.dot",
			"no_repeat_graph", g, EdgePosLab);
	omnigraph::WriteSimple(work_tmp_dir + "repeats_resolved_before_poslab.dot",
			"no_repeat_graph", g, EdgePosLab);
	PairedInfoIndex<Graph> clustered_index(g);
	if(paired_mode) {
		DistanceEstimator<Graph> estimator(g, paired_index, insert_size, max_read_length, 10, 10, 75);
		estimator.Estimate(clustered_index);
	}

	INFO("Building de Bruijn graph finished");

	INFO("Writing graphs");
	RealIdGraphLabeler<Graph> IdTrackLabelerBefore(g, IntIds);
	omnigraph::WriteSimple(
			output_folder + "repeats_resolved_before.dot",
			"no_repeat_graph", g, IdTrackLabelerBefore);
	printGraph(g, IntIds, work_tmp_dir + "graph", clustered_index, EdgePos);
	printGraph(g, IntIds, output_folder + "graph", clustered_index, EdgePos);
}



template<size_t k>
void BuildDeBruijnGraph(Graph& g,  PairedInfoIndex<Graph>& paired_index, EdgeIndex<k + 1, Graph>& index,
		Sequence& sequence) {
	// check config.hpp parameters
	if (K % 2 == 0) {
		FATAL("K in config.hpp must be odd!\n");
	}
	checkFileExistenceFATAL(CONFIG_FILENAME);

	// read configuration file (dataset path etc.)
	string input_dir = CONFIG.read<string>("input_dir");
	string dataset = CONFIG.read<string>("dataset");
	string output_root = CONFIG.read<string>("output_dir");
	string output_dir_suffix = MakeLaunchTimeDirName()+ "." + dataset + "/";
	string output_dir = output_root + output_dir_suffix;
	string work_tmp_dir = output_root + "tmp/";
//	std::cout << "here " << mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH| S_IWOTH) << std::endl;
	mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	string genome_filename = input_dir
			+ CONFIG.read<string>("reference_genome");
	string reads_filename1 = input_dir + CONFIG.read<string>(dataset + "_1");
	string reads_filename2 = input_dir + CONFIG.read<string>(dataset + "_2");
	checkFileExistenceFATAL(genome_filename);
	checkFileExistenceFATAL(reads_filename1);
	checkFileExistenceFATAL(reads_filename2);

	size_t insert_size = CONFIG.read<size_t>(dataset + "_IS");
	size_t max_read_length = 100; //CONFIG.read<size_t> (dataset + "_READ_LEN");
	int dataset_len = CONFIG.read<int>(dataset + "_LEN");
	bool paired_mode = CONFIG.read<bool>("paired_mode");
    bool rectangle_mode  = CONFIG.read<bool>("rectangle_mode");
	bool etalon_info_mode = CONFIG.read<bool>("etalon_info_mode");
	bool from_saved = CONFIG.read<bool>("from_saved_graph");
	// typedefs :)
	//typedef MateReader<Read, ireadstream>::type ReadStream;
	//typedef PairedReader<ireadstream> PairedReadStream;
	//typedef RCReaderWrapper<PairedReadStream, PairedRead> RCStream;
	// read data ('reads')
	//const string reads[2] = {reads_filename1, reads_filename2};
	//ReadStream reader(reads);
	//PairedReadStream pairStream(reader, insert_size);
	//RCStream rcStream(pairStream);

	// typedefs :)
  typedef io::Reader<io::SingleRead> ReadStream;
  typedef io::Reader<io::PairedRead> PairedReadStream;
  typedef io::RCReaderWrapper<io::PairedRead> RCStream;

	// read data ('reads')
  PairedReadStream pairStream(std::pair<std::string, 
                              std::string>(reads_filename1,
                                           reads_filename2), 
                              insert_size);
	RCStream rcStream(&pairStream);


	// read data ('genome')
	std::string genome;
	{
		ReadStream genome_stream(genome_filename);
		io::SingleRead full_genome;
		genome_stream >> full_genome;
		genome = full_genome.GetSequenceString().substr(0, dataset_len); // cropped
	}
	sequence = Sequence(genome);
	// assemble it!

	BuildDeBruijnGraph<K, RCStream>(rcStream, sequence, paired_mode, rectangle_mode, etalon_info_mode,
			from_saved, insert_size, max_read_length, output_dir, work_tmp_dir, g, paired_index, index);
}

template<size_t k>
void LoadFromFile(std::string fileName, Graph& g,  PairedInfoIndex<Graph>& paired_index, EdgeIndex<k + 1, Graph>& index,
		Sequence& sequence) {

	string input_dir = CONFIG.read<string>("input_dir");
	string dataset = CONFIG.read<string>("dataset");
	string genome_filename = input_dir
			+ CONFIG.read<string>("reference_genome");
	checkFileExistenceFATAL(genome_filename);
	int dataset_len = CONFIG.read<int>(dataset + "_LEN");

	  typedef io::Reader<io::SingleRead> ReadStream;
	  typedef io::Reader<io::PairedRead> PairedReadStream;
	  typedef io::RCReaderWrapper<io::PairedRead> RCStream;

	// read data ('genome')
	std::string genome;
	{
		ReadStream genome_stream(genome_filename);
		io::SingleRead full_genome;
		genome_stream >> full_genome;
		genome = full_genome.GetSequenceString().substr(0, dataset_len); // cropped
	}
	sequence = Sequence(genome);

	INFO("Reading graph");
	//IdTrackHandler<Graph> IntIds(g);
	IdTrackHandler<Graph> conj_IntIds(g);

	scanConjugateGraph(g, conj_IntIds,	fileName, paired_index);
}

template<size_t k>
void AddEtalonInfo(Graph& g, EdgeIndex<k+1, Graph>& index, const Sequence& genome, PairedInfoIndices& pairedInfos) {
	size_t libCount = LC_CONFIG.read<size_t>("lib_count");

	for (size_t i = 1; i <= libCount; ++i) {
		std::string num = ToString<size_t>(i);
		size_t insertSize = LC_CONFIG.read<size_t>("insert_size_" + num);
		size_t readSize = LC_CONFIG.read<size_t>("read_size_" + num);
		INFO("Generating info with read size " << readSize << ", insert size " << insertSize);

		pairedInfos.push_back(PairedInfoIndexLibrary(readSize, insertSize, new PairedInfoIndex<Graph>(g, 0)));
		FillEtalonPairedIndex<k> (g, *pairedInfos.back().pairedInfoIndex, index, insertSize, readSize, genome);
	}
}

void DeleteEtalonInfo(PairedInfoIndices& pairedInfos) {
	while (pairedInfos.size() > 1) {
		delete pairedInfos.back().pairedInfoIndex;
		pairedInfos.pop_back();
	}
}

} // namespace long_contigs

int main() {
	using namespace long_contigs;

	checkFileExistenceFATAL(LC_CONFIG_FILENAME);
	checkFileExistenceFATAL(CONFIG_FILENAME);

	double MIN_COVERAGE = LC_CONFIG.read<double>("min_coverage");

	Graph g(K);
	EdgeIndex<K + 1, Graph> index(g);
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
		BuildDeBruijnGraph<K>(g, pairedIndex, index, sequence);
	}
	else {
		mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

		LoadFromFile<K>(LC_CONFIG.read<std::string>("graph_file_" + dataset), g, pairedIndex, index, sequence);
	}

	Path<Graph::EdgeId> path1 = FindGenomePath<K> (sequence, g, index);
	Path<Graph::EdgeId> path2 = FindGenomePath<K> (!sequence, g, index);

	PairedInfoIndexLibrary basicPairedLib(LC_CONFIG.read<size_t>("read_size"), LC_CONFIG.read<size_t>("insert_size"), &pairedIndex);
	pairedInfos.push_back(basicPairedLib);

	if (LC_CONFIG.read<bool>("etalon_info_mode")) {
		AddEtalonInfo<K>(g, index, sequence, pairedInfos);
	}

	FindSeeds(g, rawSeeds);
	INFO("Seeds found");
//	PrintPathsShort(g, rawSeeds);
	WriteGraphWithPathsSimple(output_dir + "raw_seeds.dot", "raw_seeds", g, rawSeeds, path1, path2);

	RemoveSubpaths(g, rawSeeds, seeds);
	INFO("Removed subseeds");
//	PrintPathsShort(g, seeds);
	WriteGraphWithPathsSimple(output_dir + "no_dupl_seeds.dot", "no_dupl_seeds", g, seeds, path1, path2);

	FilterLowCovered(g, seeds, MIN_COVERAGE);
	INFO("Seeds filtered");
//	PrintPathsShort(g, seeds);

	size_t found = PathsInGenome<K>(g, index, sequence, seeds, path1, path2, true);
	INFO("Good seeds found " << found << " in total " << seeds.size());
	INFO("All seed coverage " << PathsCoverage(g, seeds));

	WriteGraphWithPathsSimple(output_dir + "seeds.dot", "seeds", g, seeds, path1, path2);

	MakeKeyPause();
	FindPaths(g, seeds, pairedInfos, paths);
	//INFO("Final paths");
	//PrintPathsShort(g, paths);

	std::vector<BidirectionalPath> result;
	RemoveSubpaths(g, paths, result);
	INFO("Removed subpaths");
	//PrintPathsShort(g, result);

	found = PathsInGenome<K>(g, index, sequence, result, path1, path2);
	INFO("Good paths found " << found << " in total " << result.size());
	INFO("Path coverage " << PathsCoverage(g, result));
	INFO("Path length coverage " << PathsLengthCoverage(g, result));

	WriteGraphWithPathsSimple(output_dir + "paths.dot", "paths", g, result, path1, path2);

//	INFO("Genome paths");
//	PrintPathWithVertices(g, path1);
//	PrintPathWithVertices	(g, path2);

	DeleteEtalonInfo(pairedInfos);
	return 0;
}

