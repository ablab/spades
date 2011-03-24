#include "common.hpp"

#include "constructHashTable.hpp"
#include "graphConstruction.hpp"
#include "graphSimplification.hpp"
#include "graphio.hpp"
#include "readTracing.hpp"
#include "sequence.hpp"

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
		sprintf(str, "data/graph.txt");
		save(str,graph);
		outputLongEdges(graph.longEdges, graph, "data/beforeExpand.dot");
	}

	if (useExpandDefinite){
		INFO("Expand definite...");
		if (!needGraph){
			sprintf(str, "data/graph.txt");
			load(str,graph);
			graph.RebuildVertexMap();
			graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);

		}
		expandDefinite(graph.longEdges, graph, graph.VertexCount, true);
		outputLongEdges(graph.longEdges, graph, "data/afterExpand.dot");
		outputLongEdgesThroughGenome(graph, "data/afterExpand_g.dot");
		sprintf(str, "data/expandedGraph.txt");
		save(str,graph);
	}

	if (useTraceReads){
		INFO("Trace reads...");
		if (!useExpandDefinite){
			sprintf(str, "data/expandedGraph.txt");
			load(str,graph);
			graph.RebuildVertexMap();
			graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
		}
		traceReads(graph.verts, graph.longEdges, graph, graph.VertexCount, graph.EdgeId);
		outputLongEdges(graph.longEdges,"data/ReadsTraced.dot");
		outputLongEdgesThroughGenome(graph, "data/ReadsTraced_g.dot");
		sprintf(str, "data/tracedGraph.txt");
		save(str,graph);
	}

	if (useProcessLower){
		INFO("Process lowers");

		if (!useTraceReads){
			sprintf(str, "data/tracedGraph.txt");
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
		outputLongEdges(graph.longEdges,"data/afterLowers.dot");
		outputLongEdgesThroughGenome(graph, "data/afterLowers_g.dot");

		graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
		outputLongEdges(graph.longEdges, graph, "data/afterLowers_info.dot");
		sprintf(str, "data/afterLowerGraph.txt");
		save(str,graph);
	}
	cerr << "\n Finished";
	INFO("Finished");
}



int main() {
	init();
/*	freopen(error_log.c_str(), "w",stderr);
	char str[100];
	sprintf(str, "data/tracedGraph.txt");
	load(str,graph);
	graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
	PairThreader pg(graph,1);
	for (longEdgesMap::iterator it = graph.longEdges.begin(); it != graph.longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			vector<pair<int, Edge *> > vp = pg.threadLower(it->second);
			forn(i,vp.size())
				cerr<<"edge "<<it->first<<" may jump into "<<vp[i].second->EdgeId<<" dist "<<vp[i].first<<endl;

		}
	}


*/
		run();
	return 0;
}
