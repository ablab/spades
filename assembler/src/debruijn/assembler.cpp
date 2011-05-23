/*
 * Assembler Main
 */
#include "launch.hpp"
#include "config.hpp"
#include "libs/ConfigFile/ConfigFile.h"

void RunEdgeGraphTool() {

	// config parsing... dataset path etc.
	ConfigFile config(CONFIG_FILE);
	string input_dir = config.read<string>("input_dir");
	string debruijn_dir = config.read<string>("debruijn_dir");
	string dataset = config.read<string>("dataset");
	string reference_genome = input_dir + "/" + config.read<string>("reference_genome");
	string reads1 = input_dir + "/" + config.read<string>(dataset + "_1");
	string reads2 = input_dir + "/" + config.read<string>(dataset + "_2");
	int insert_size = config.read<int>(dataset + "_IS");
	int dataset_len = config.read<int>(dataset + "_LEN");

	// typedefs :)
	typedef StrobeReader<2, Read, ireadstream> ReadStream;
	typedef PairedReader<ireadstream> PairedStream;
	typedef RCReaderWrapper<PairedStream, PairedRead> RCStream;

	// read data
	const string reads[2] = {reads1, reads2};
	ReadStream reader(reads);
	PairedStream pairStream(reader, insert_size);
	RCStream rcStream(pairStream);

	// assemble
	ireadstream genome_stream(reference_genome);
	Read genome;
	genome_stream >> genome;
	genome_stream.close();
	debruijn_graph::EdgeGraphTool<K, RCStream>(rcStream, genome.getSequenceString().substr(0, dataset_len), debruijn_dir);
	reader.close();
}

int main() {
	RunEdgeGraphTool();
	return 0;
}
