/*
 * Assembler Main
 */
#include "launch.hpp"
#include "config.hpp"

void RunEdgeGraphTool() {
	typedef StrobeReader<2, Read, ireadstream> ReadStream;
	typedef PairedReader<ireadstream> PairedStream;
	typedef RCReaderWrapper<PairedStream> RCStream;

	const string reads[2] = {tr1::get<0>(INPUT), tr1::get<1>(INPUT)};
	ReadStream reader(reads);
	PairedStream pairStream(reader, tr1::get<2>(INPUT));
	RCStream rcStream(pairStream);

	ireadstream genome_stream(ECOLI_FILE);
	Read genome;
	genome_stream >> genome;
	edge_graph::EdgeGraphTool<K, RCStream>(rcStream, genome.getSequenceString().substr(0, tr1::get<3>(INPUT)), DE_BRUIJN_DATA_FOLDER);
	reader.close();
	genome_stream.close();
}

int main() {
	RunEdgeGraphTool();
	return 0;
}
