
/////////////////
//for read generator
//#define SUBSTR_LENGTH 10000
//#define COVERAGE 30
//#define R 35
/////////////////

#include "config.hpp"
#include "visualization_utils.hpp"
#include "ireadstream.hpp"

void RunEdgeGraphTool() {
	typedef StrobeReader<2, Read, ireadstream> ReadStream;
	typedef PairedReader<ireadstream> PairedStream;
	typedef RCReaderWrapper<PairedStream> RCStream;

	//	const tr1::tuple<string, string, int> input = ;
	const string reads[2] = {tr1::get<0>(INPUT), tr1::get<1>(INPUT)};
	ReadStream reader(reads);
	PairedStream pairStream(reader, tr1::get<2>(INPUT));
	RCStream rcStream(pairStream);

	ireadstream genome_stream(ECOLI_FILE);
	Read genome;
	genome_stream >> genome;
	edge_graph::EdgeGraphTool(rcStream,  genome.getSequenceString().substr(0, tr1::get<2>(INPUT)));
	reader.close();
	genome_stream.close();
}

int main() {
	RunEdgeGraphTool();
	return 0;
}
