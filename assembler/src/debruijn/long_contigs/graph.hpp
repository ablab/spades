/*
 * graph.hpp
 *
 *  Created on: Jul 11, 2011
 *      Author: andrey
 */

#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include "../launch.hpp"

namespace long_contigs {


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
		const Sequence& genome, bool paired_mode, bool etalon_info_mode,
		bool from_saved, size_t insert_size, size_t max_read_length,
		const string& output_folder, const string& work_tmp_dir,
		Graph& g) {

	INFO("Building de Bruijn graph started");
	INFO("Paired mode: " << (paired_mode ? "Yes" : "No") );
	INFO("Etalon paired info mode: " << (etalon_info_mode ? "Yes" : "No"))
	INFO("From file: " << (from_saved ? "Yes" : "No"))
	mkdir(work_tmp_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	EdgeIndex<k + 1, Graph> index(g);
	IdTrackHandler<Graph> IntIds(g);
	// if it's not paired_mode, then it'll be just unused variable -- takes O(1) to initialize from graph
	PairedInfoIndex<Graph> paired_index(g);

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

		SimplifyGraph<k> (g, index, 3, genome, output_folder);

		ProduceInfo<k> (g, index, genome,
				output_folder + "simplified_graph.dot", "simplified_graph");

		WriteGraphComponents<k> (g, index, genome,
				output_folder + "graph_components" + "/", "graph.dot",
				"graph_component", insert_size);
		if (paired_mode) {
			ProducePairedInfo(g, insert_size, max_read_length, paired_index,
					output_folder);
		}
	}
	INFO("Building de Bruijn graph started");
}

void BuildDeBruijnGraph(Graph& g) {
	// check config.hpp parameters
	if (K % 2 == 0) {
		FATAL("K in config.hpp must be odd!\n");
	}
	std::cout << CONFIG_FILENAME << std::endl;
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
	bool etalon_info_mode = CONFIG.read<bool>("etalon_info_mode");
	bool from_saved = CONFIG.read<bool>("from_saved_graph");
	// typedefs :)
	typedef MateReader<Read, ireadstream>::type ReadStream;
	typedef PairedReader<ireadstream> PairedReadStream;
	typedef RCReaderWrapper<PairedReadStream, PairedRead> RCStream;

	// read data ('reads')
	const string reads[2] = {reads_filename1, reads_filename2};
	ReadStream reader(reads);
	PairedReadStream pairStream(reader, insert_size);
	RCStream rcStream(pairStream);

	// read data ('genome')
	std::string genome;
	{
		ireadstream genome_stream(genome_filename);
		Read full_genome;
		genome_stream >> full_genome;
		genome = full_genome.getSequenceString().substr(0, dataset_len); // cropped
	}
	// assemble it!

	BuildDeBruijnGraph<K, RCStream>(rcStream, Sequence(genome), paired_mode, etalon_info_mode, from_saved, insert_size, max_read_length, output_dir, work_tmp_dir, g);
}



}

#endif /* GRAPH_HPP_ */
