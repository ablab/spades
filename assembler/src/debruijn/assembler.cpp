/*
 * Main of Assembler!
 */

#include "config.hpp"
#include "visualization_utils.hpp"
#include "ireadstream.hpp"

#include "paired_info.hpp"
#include "edge_graph_constructor.hpp"
#include "tip_clipper.hpp"
#include "bulge_remover.hpp"
#include "coverage_handler.hpp"


/*
 * EdgeGraphTool
 */

namespace edge_graph {

typedef de_bruijn::DeBruijnPlus<K + 1, EdgeId> DeBruijn;
typedef de_bruijn::DeBruijnPlus<K + 1, EdgeId> Index;
typedef de_bruijn::Path<EdgeId> Path;
typedef de_bruijn::PairedInfoIndex<EdgeGraph> PairedIndex;

void CountStats(const EdgeGraph& g) {
	INFO("Counting stats");
	de_bruijn::DFS<EdgeGraph> dfs(g);
	de_bruijn::SimpleStatCounter<EdgeGraph> stat_c;
	dfs.Traverse(&stat_c);
	INFO("Vertex count=" << stat_c.v_count() << "; Edge count="
			<< stat_c.e_count());
	INFO("Stats counted");
}

void WriteToDotFile(const EdgeGraph &g, const string& file_name,
		string graph_name,
		de_bruijn::Path<EdgeId> path = de_bruijn::Path<EdgeId>()) {
	INFO("Writing graph '" << graph_name << "' to file " << file_name);
	WriteToFile(DE_BRUIJN_DATA_FOLDER + file_name, graph_name, g, path);
	INFO("Graph '" << graph_name << "' written to file " << file_name);
}

Path FindGenomePath(const string &genome, const EdgeGraph& g,
		const Index& index) {
	de_bruijn::SimpleSequenceMapper<K, EdgeGraph> srt(g, index);
	return srt.MapSequence(Sequence(genome));
}

void ProduceInfo(const EdgeGraph& g, const Index& index, const string& genome,
		const string& file_name, const string& graph_name) {
	CountStats(g);
	Path path = FindGenomePath(genome, g, index);
	WriteToDotFile(g, file_name, graph_name, path);
}

/*template<class ReadStream>
 void ConstructUncondensedGraph(DeBruijn& debruijn, ReadStream& stream) {
 INFO("Constructing DeBruijn graph");
 debruijn.ConstructGraphFromStream(stream);
 INFO("DeBruijn graph constructed");
 }*/

template<class ReadStream>
void CondenseGraph(DeBruijn& debruijn, EdgeGraph& g, Index& index,
		ReadStream& stream, const string& genome) {
	INFO("Condensing graph");
	EdgeGraphConstructor<K> g_c(debruijn);
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");

	INFO("Counting coverage");
	CoverageCounter<K, EdgeGraph> cc(g, index);
	typedef SimpleReaderWrapper<ReadStream> SimpleRCStream;
	SimpleRCStream unitedStream(stream);
	cc.CountCoverage(unitedStream);
	INFO("Coverage counted");

	ProduceInfo(g, index, genome, "edge_graph.dot", "edge_graph");
}

void ClipTips(EdgeGraph &g, Index &index, const string& genome,
		string dotFileName) {
	INFO("Clipping tips");
	TipComparator<EdgeGraph> comparator(g);
	TipClipper<EdgeGraph, TipComparator<EdgeGraph> > tc(comparator);
	tc.ClipTips(g);
	INFO("Tips clipped");

	ProduceInfo(g, index, genome, dotFileName, "no_tip_graph");
}

void RemoveBulges(EdgeGraph &g, Index &index, const string& genome,
		string dotFileName) {
	INFO("Removing bulges");
	de_bruijn::BulgeRemover<EdgeGraph> bulge_remover;
	bulge_remover.RemoveBulges(g);
	INFO("Bulges removed");

	ProduceInfo(g, index, genome, dotFileName, "no_bulge_graph");
}

template<class ReadStream>
void FillPairedIndex(PairedIndex& paired_info_index, ReadStream& stream,
		Index& index) {
	INFO("Counting paired info");
	paired_info_index.FillIndex<K, ReadStream> (index, stream);
	INFO("Paired info counted");
}

template<class ReadStream>
void EdgeGraphTool(ReadStream& stream, const string& genome) { // actually, it's assembler :)
	INFO("Edge graph construction tool started");

	INFO("Constructing DeBruijn graph");
	SimpleReaderWrapper<ReadStream> unitedStream(stream);
	DeBruijn debruijn(unitedStream);
	INFO("DeBruijn graph constructed");

	EdgeGraph g(K);
	de_bruijn::DeBruijnPlus<K + 1, EdgeId> &index = debruijn;
	EdgeHashRenewer<K + 1, EdgeGraph> index_handler(g, index);
	g.AddActionHandler(&index_handler);

	de_bruijn::CoverageHandler<EdgeGraph> coverageHandler(g);
	g.AddActionHandler(&coverageHandler);

	stream.reset();
	CondenseGraph<ReadStream> (debruijn, g, index, stream, genome);

	stream.reset();
	PairedIndex paired_info_index(g, I);

	FillPairedIndex(paired_info_index, stream, index);
	ClipTips(g, index, genome, "tips_clipped.dot");

//	paired_info_index.OutputData();

	RemoveBulges(g, index, genome, "bulges_removed.dot");

	g.RemoveActionHandler(&index_handler);
	g.RemoveActionHandler(&coverageHandler);

	INFO("Tool finished")
}

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
	edge_graph::EdgeGraphTool(rcStream,  genome.getSequenceString().substr(0, tr1::get<3>(INPUT)));
	reader.close();
	genome_stream.close();
}

} // namespace edge_graph


/*
 * Assembler Main
 */

int main() {
	edge_graph::RunEdgeGraphTool();
	return 0;
}
