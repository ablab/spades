#include "common.hpp"

#include "constructHashTable.hpp"
#include "graphConstruction.hpp"
#include "graphSimplification.hpp"
#include "graphio.hpp"
#include "readTracing.hpp"
#include "sequence.hpp"
#include "readsReformatter.hpp"
#include "read_generator.hpp"

using namespace paired_assembler;

LOGGER("p.main");

PairedGraph graph;

void init() {
	initConstants(ini_file);
	initGlobal();
	freopen(error_log.c_str(), "w",stderr);
	INFO("Constants inited...");
	cerr << l << " " << k;
}

void run() {
	char str[100];
//	forgetQualityPairedData("I:/bioinf/eas20_8/s_6_1.cor.fastq.gz", "I:/bioinf/eas20_8/s_6_2.cor.fastq.gz", "data/paireddebruijn/reads_100_200_corr.txt" );
//	forgetQualityPairedData("data/paireddebruijn/s_6_1.fastq.gz", "data/paireddebruijn/s_6_2.fastq.gz", "/media/605005E05005BDB2/data/realreads.txt" );
//	LOG_ASSERT(1 == 0, "Something wrong");

	if (needPairs) {
		cerr << endl << " constructing pairs" << endl;
		readsToPairs(parsed_reads, parsed_k_l_mers);
	}
	if (needLmers) {
		cerr << endl << " constructing Lmers" << endl;
		pairsToLmers(parsed_k_l_mers, parsed_l_mers);
	}
	if (needSequences) {
		cerr << endl << " constructing Sequences" << endl;
		pairsToSequences(parsed_k_l_mers, parsed_l_mers, parsed_k_sequence);
	}
	//	map<>sequencesToMap(parsed_k_sequence);

	if (needGraph) {
		cerr << endl << " constructing Graph" << endl;
		constructGraph(graph);
		sprintf(str, "data/paireddebruijn/graph.txt");
		save(str,graph);
		outputLongEdges(graph.longEdges, graph, "data/paireddebruijn/beforeExpand.dot");
	}

	if (useExpandDefinite){
		INFO("Expand definite...");
		if (!needGraph){
			sprintf(str, "data/paireddebruijn/graph.txt");
			load(str,graph);
			graph.removeLowCoveredEdges(graph.longEdges, 3);
			graph.RebuildVertexMap();
			graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
		}
//		expandDefinite(graph.longEdges, graph, graph.VertexCount, false);
		expandDefinite(graph.longEdges, graph, graph.VertexCount, true);
		outputLongEdges(graph.longEdges, graph, "data/paireddebruijn/afterExpand.dot");
//		outputLongEdgesThroughGenome(graph, "data/paireddebruijn/afterExpand_g.dot");
		sprintf(str, "data/paireddebruijn/expandedGraph.txt");
		save(str,graph);
	}

	if (useTraceReads){
		INFO("Trace reads...");
		if (!useExpandDefinite){
			sprintf(str, "data/paireddebruijn/expandedGraph.txt");
			INFO("Loading graph...");
			load(str,graph);
			INFO("Graph loaded!");
			graph.RebuildVertexMap();
			graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
		}
		traceReads(graph.verts, graph.longEdges, graph, graph.VertexCount, graph.EdgeId);
		outputLongEdges(graph.longEdges,"data/paireddebruijn/ReadsTraced.dot");
	//	outputLongEdgesThroughGenome(graph, "data/paireddebruijn/ReadsTraced_g.dot");
		sprintf(str, "data/paireddebruijn/tracedGraph.txt");
		save(str,graph);
	}

	if (useProcessLower){
		INFO("Process lowers");

		if (!useTraceReads){
			sprintf(str, "data/paireddebruijn/tracedGraph.txt");
			INFO("Load");
			load(str,graph);
			INFO("Rebuild");
			graph.RebuildVertexMap();
		}

		graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
		while (processLowerSequence(graph.longEdges, graph, graph.VertexCount))
		{
			graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
			expandDefinite(graph.longEdges , graph, graph.VertexCount);
			INFO("one more");
		}
		outputLongEdges(graph.longEdges,"data/paireddebruijn/afterLowers.dot");
		outputLongEdgesThroughGenome(graph, "data/paireddebruijn/afterLowers_g.dot");

		graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
		outputLongEdges(graph.longEdges, graph, "data/paireddebruijn/afterLowers_info.dot");
		sprintf(str, "data/paireddebruijn/afterLowerGraph.txt");
		save(str,graph);
	}
	if (useExtractDefinite){
		INFO("extractDefinite RIGHT Start");
		if (!useProcessLower){
			sprintf(str, "data/paireddebruijn/afterLowerGraph.txt");
			load(str,graph);
			graph.RebuildVertexMap();
			DEBUG("loaded");
		}
		graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
		extractDefinite(graph, RIGHT);
		outputLongEdges(graph.longEdges, graph,  "data/paireddebruijn/afterExtractDefinite1.dot");
	//	outputLongEdgesThroughGenome(graph, "data/paireddebruijn/afterExtractDefinite1_g.dot");

		INFO ("extractDefinite LEFT Start");
			graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
			extractDefinite(graph, LEFT);
			outputLongEdges(graph.longEdges, graph,  "data/paireddebruijn/afterExtractDefinite2.dot");
			//outputLongEdgesThroughGenome(graph, "data/paireddebruijn/afterExtractDefinite2_g.dot");


	}
	cerr << "\n Finished";
	INFO("Finished");
}

void generateReads(string fileName, string genomeFileName, int insertLength, int coverage) {
	ofstream os;
	os.open(fileName.c_str());
	Sequence genome = readGenomeFromFile(genomeFileName);
	stringstream ss;
	ss << genome;
	ReadGenerator<100, 2, int, RandomPositionChooser> gen(ss.str(), coverage, insertLength);
	gen.setErrorProbability(0);
	gen.setMaxInsertLengthError(10);
	strobe_read<100, 2> readPair;
	while(!gen.eof()) {
		 gen >> readPair;
		 os << readPair[0] << " " << readPair[1] << endl;
	}
	os.close();
}

int main() {
//	generateReads("data/paireddebruijn/someFile.txt", "data/input/MG1655-K12_cut.fasta", 20, 20);
	init();
	run();
	return 0;
}
