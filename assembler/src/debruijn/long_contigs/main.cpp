/*
 * main.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: andrey
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../launch.hpp"
#include "../config.hpp"
#include "../../common/logging.hpp"
#include "../../common/simple_tools.hpp"

#include "long_contigs.hpp"

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
			typedef SimpleReaderWrapper<ReadStream> UnitedStream;
			UnitedStream united_stream(stream);
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
		ireadstream genome_stream(genome_filename);
		Read full_genome;
		genome_stream >> full_genome;
		genome = full_genome.getSequenceString().substr(0, dataset_len); // cropped
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

	// read data ('genome')
	std::string genome;
	{
		ireadstream genome_stream(genome_filename);
		Read full_genome;
		genome_stream >> full_genome;
		genome = full_genome.getSequenceString().substr(0, dataset_len); // cropped
	}
	sequence = Sequence(genome);

	INFO("Reading graph");
	//IdTrackHandler<Graph> IntIds(g);
	IdTrackHandler<Graph> conj_IntIds(g);

	scanConjugateGraph(g, conj_IntIds,	fileName, paired_index);
}


int main() {
	using namespace long_contigs;

	LoadLCConstants();

	Graph g(K);
	EdgeIndex<K + 1, Graph> index(g);
	PairedInfoIndex<Graph> paired_index(g, 0);
	Sequence sequence("");

	std::vector<BidirectionalPath> seeds;
	std::vector<BidirectionalPath> paths;

	if (!LC_CONFIG.read<bool>("from_file")) {
		BuildDeBruijnGraph<K>(g, paired_index, index, sequence);
	}
	else {
		LoadFromFile<K>(LC_CONFIG.read<std::string>("graph_file"), g, paired_index, index, sequence);
	}

	FindSeeds(g, seeds);

	for(auto iter = seeds.begin(); iter != seeds.end(); ++iter) {
		PrintPath(g, *iter);
	}

	INFO("Seeds length");
	PrintPathLengthStats(g ,seeds);
	INFO("Seeds edges count");
	PrintPathEdgeLengthStats(seeds);
	INFO("Seeds coverage");
	PrintPathCoverage(g, seeds);

	FindPaths(g, seeds, paired_index, paths);

	INFO("Final paths");
	for(auto path = paths.begin(); path != paths.end(); ++path) {
		PrintPath(g, *path);
	}

	size_t found = PathsInGenome<K>(g, index, sequence, paths);
	INFO("Good paths found " << found << " in total " << paths.size());
	INFO("Path coverage " << PathsCoverage(g, paths));
	return 0;
}

