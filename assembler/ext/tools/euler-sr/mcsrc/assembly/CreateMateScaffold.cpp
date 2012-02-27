/***************************************************************************
 * Title:          CreateMateScaffold.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "MateLibrary.h"
#include "PathLib.h"
#include <vector>
#include "IntegralTupleStatic.h"


typedef std::vector<std::list<ssize_t > > ScaffoldList;


	
ssize_t SearchContainedScaffold(ScaffoldList &scaffoldList, ssize_t src, ssize_t dest, ssize_t maxDepth) {
	if (maxDepth == 0)
		return 0;

	if (std::find(scaffoldList[src].begin(), 
								scaffoldList[src].end(), dest) != scaffoldList[src].end()) {
		return 1;
	}
	else {
		std::list<ssize_t>::iterator it;
		for (it = scaffoldList[src].begin(); it != scaffoldList[src].end(); ++it) {
			if (SearchContainedScaffold(scaffoldList, *it, dest, maxDepth - 1)) {
				// return this path early 
				return 1;
			}
		}
		return 0;
	}
}

void RemoveContainedScaffolds(ScaffoldList &scaffoldList) {
	ssize_t e;
	std::list<ssize_t>::iterator it, erased;
	for (e = 0; e < scaffoldList.size(); e++) {
		it = scaffoldList[e].begin();
		while (it != scaffoldList[e].end()) {
			if (*it == e) {
				erased = it;
				++it;
				scaffoldList[e].erase(erased);
			}
			else {
				++it;
			}
		}
	}
	for (e = 0; e < scaffoldList.size(); e++) { 
		std::list<ssize_t>::iterator it1, it2, erasedIt;

		if (scaffoldList[e].size() > 1) {
			it1 = scaffoldList[e].begin();
			while (it1 != scaffoldList[e].end()) {
				it2 = scaffoldList[e].begin();
				while(it2 != scaffoldList[e].end()) {
					if (*it1 == *it2) {
						++it2;
						continue;
					}

					//UNUSED// ssize_t removedIt1 = 0;
					if (SearchContainedScaffold(scaffoldList, *it2, *it1, 5)) {
						erasedIt = it2;
						++it2;
						std::cout << "erasing " << e << " " << *erasedIt << " ---> " << *it1 << std::endl;
						scaffoldList[e].erase(erasedIt);
					}
					else {
						++it2;
					}
				}
				++it1;
			}
		}
	}
}


void PrintUsage() {
	std::cout << " usage: createMateScaffold graphIn mateTable ruleFile graphOut " << std::endl
						<< "  [-minPathCount c]  Remove paths with count less than 'c' " << std::endl
						<< "  [-minMatePairCount c] Remove paths with mate count less than 'c'" << std::endl;
}

int main(int argc, char* argv[]) {
	
	std::string graphFileName, mateTableName, graphOutName, ruleFileName;
	if (argc < 5) {
		PrintUsage();
		exit(0);
	}
	graphFileName = argv[1];
	mateTableName = argv[2];
	ruleFileName  = argv[3];
	graphOutName  = argv[4];
	ssize_t minPathCount = 2;
	ssize_t minMatePairCount = 2;
	int argi = 5;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minPathCount") == 0) {
			minPathCount = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minMatePairCount") == 0) {
			minMatePairCount = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(0);
		}
		++argi;
	}

	IntervalGraph graph;
	ReadMateList  mateList;
	std::cout << "Reading interval graph." << std::endl;
	int vertexSize;
	ReadIntervalGraph(graphFileName, graph, vertexSize);
	std::cout << "Reading mate table." << std::endl;
	ReadMateTable(mateTableName, mateList);

	RuleList rules;
	ParseRuleFile(ruleFileName, rules);

	//UNUSED// PathIntervalList &paths       = graph.paths;
	//UNUSED// PathLengthList   &pathLengths = graph.pathLengths;
	TEdgeList        &edges       = graph.edges;
	TVertexList      &vertices    = graph.vertices;

	/*
	 * Now collect the mate-edges from each edge to see 
	 * if there is a unique edge to scaffold it with.
	 *
	 */
	TEdgeList connectingEdges;
	ssize_t newEdgeIndex = edges.size();
	ssize_t newEdgeCount = 0;
	ssize_t v;
	// Now disconnect all vertices.  They will be scaffolded later on.
	TVertexList newVertices;
	for (v = 0; v < vertices.size(); v++ ){
		ssize_t vi, vo;
		for (vi = 0; vi < vertices[v].in.size(); vi++ ){
			/*			if (vertices[v].in[vi] != -1) {
				edges[vertices[v].in[vi]].dest = -1;
			*/
				vertices[v].in[vi] = -1;
				//			}
		}
		for (vo = 0; vo < vertices[v].out.size(); vo++ ){
			/*			if (vertices[v].out[vo] != -1) {
				edges[vertices[v].out[vo]].src = -1;
			*/
			vertices[v].out[vo] = -1;
				//			}
		}
	}

	vertices.resize(edges.size() * 2);
	v = 0;
	std::vector<std::list<ssize_t> > scaffold;
	scaffold.resize(edges.size());
	ssize_t e;
	for (e = 0; e < edges.size(); e++) {
		vertices[v].AddOutEdge(e);
		vertices[v+1].AddInEdge(e);
		edges[e].src = v;
		edges[e].dest = v + 1;
		vertices[v].vertexSize = vertexSize;
		vertices[v+1].vertexSize = vertexSize;
		v+=2;
	}
	ssize_t curRule;
	ssize_t expMateSep;
	double stddevMateSep;
	cout << "num rules: " << rules.size() << endl;
	for (curRule = 0; curRule < rules.size(); curRule++) {

		ComputeMatePairLengthDistribution(graph, mateList, curRule, expMateSep, stddevMateSep); 
		cout << "cur rule: " << curRule << " exp sep: " << expMateSep << " stddev " << stddevMateSep << endl;


		for (e = 0; e < edges.size(); e++) {
			MateEdgeMap mateEdges;
			CollectMateEdges(graph, mateList, e, mateEdges, curRule);
			
			MateEdgeMap::iterator matesIt, matesEnd, deletedIt;
			matesEnd =  mateEdges.end();
		
			// Remove low count mate edges.
			std::cout << "for edge " << e << " (" << edges[e].index << ", " 
								<< edges[e].length  << ") collected: " << mateEdges.size() << " edges."
								<< std::endl;
			for (matesIt = mateEdges.begin(); matesIt != mateEdges.end(); ++matesIt) {
				ssize_t mateSep = ((edges[e].length - (*matesIt).second.avgStartPos) + 
											 (*matesIt).second.avgEndPos);
				std::cout << (*matesIt).first << " (" << (*matesIt).second.count << ", " 
									<< mateSep << " " << (*matesIt).second.avgStartPos << " " << (*matesIt).second.avgEndPos << " "
									<< edges[(*matesIt).first].length << ") ";
				cout << "1: " << edges[e].length << ", " << expMateSep << " 2: " << edges[(*matesIt).first].length
						 << ", " << expMateSep << " 3: " << mateSep << ", " << expMateSep + stddevMateSep << endl;
				if (edges[e].length > expMateSep and 
						edges[(*matesIt).first].length > expMateSep and
						mateSep < expMateSep + stddevMateSep ) {
					cout << "adding " << (*matesIt).first << endl;
					scaffold[e].push_back((*matesIt).first);
				}
			}
			std::cout << std::endl;
		}

		ssize_t nscaffolds = 0;
		for (e = 0; e < edges.size(); e++)  {
			if (scaffold[e].size() > 0) {
				std::list<ssize_t>::iterator sit;
				std::cout << "scaffold " << e << " " << edges[e].length << " : ";
				for (sit = scaffold[e].begin(); sit != scaffold[e].end(); ++sit) {
					std::cout << *sit << " " ;
					++nscaffolds;
				}
				std::cout << std::endl;
			}
		}
		std::cout << "before transitive edge removal found: " << nscaffolds << " scaffolds."
							<< std::endl;

		RemoveContainedScaffolds(scaffold);
	
		//UNUSED// ssize_t nremaining = 0;
		for (e = 0; e < edges.size(); e++)  {
			if (scaffold[e].size() > 0) {
				std::list<ssize_t>::iterator sit;
				std::cout << "scaffold " << e << " " << edges[e].length << " : ";
				for (sit = scaffold[e].begin(); sit != scaffold[e].end(); ++sit) {
					std::cout << *sit << " " << edges [*sit].length << endl;
				}
				std::cout << std::endl;
			}
		}

		for (e = 0; e < edges.size(); e++) {
			MateEdgeMap mateEdges;

			// These two edges should be joined by a new edge.
			TEdge newEdge;
			//		if (mateEdges.size() == 1 or mateEdges.size() == 2) {

			//UNUSED// ssize_t numSpanning, avgSrcEnd, avgDestBegin;

			std::list<ssize_t>::iterator pairedEdgeIt;
			ssize_t srcEdge, destEdge;
			srcEdge = e;
			if (scaffold[e].size() == 1){
				for (pairedEdgeIt = scaffold[e].begin();
						 pairedEdgeIt != scaffold[e].end();
						 ++pairedEdgeIt) {
					destEdge = *pairedEdgeIt;
					cout << "Joining: " << e << " (" << edges[e].length << ") with " << destEdge
							 << " (" << edges[destEdge].length << ")" << endl;
					ssize_t gapLength = 100;
					connectingEdges.push_back(newEdge);
					connectingEdges[newEdgeCount].src = edges[srcEdge].dest;
					connectingEdges[newEdgeCount].dest = edges[destEdge].src;
					connectingEdges[newEdgeCount].Init();
					// Initialize the intervals to NULL.
					connectingEdges[newEdgeCount].length = vertices[edges[srcEdge].dest].vertexSize + gapLength;
					ssize_t newEdgeLength = connectingEdges[newEdgeCount].length;				
					connectingEdges[newEdgeCount].seq.seq = new unsigned char[newEdgeLength];
					ssize_t s;
				
				
					for (s = 0; s < newEdgeLength; s++ )
						connectingEdges[newEdgeCount].seq.seq[s] = 'N';
					connectingEdges[newEdgeCount].seq.length = newEdgeLength;
					vertices[edges[srcEdge].dest].AddOutEdge(newEdgeIndex);
					vertices[edges[destEdge].src].AddInEdge(newEdgeIndex);
					++newEdgeIndex;
					++newEdgeCount;
				}
			}
		}
	}
	
	std::cout << "old number of edges: " << edges.size() << " new: " << connectingEdges.size() << std::endl;
	ssize_t n;
	for (n =0; n < connectingEdges.size(); n++ ) {
		edges.push_back(connectingEdges[n]);
		edges[edges.size() - 1].flagged = GraphEdge::Marked;
	}



	//	std::string coloredGvzName = graphOutName + ".colored.dot";
	//  GVZPrintBGraph(graph.vertices, graph.edges, coloredGvzName);
	
	graph.CondenseSimplePaths();
	
	std::string bGraphOutName = graphOutName + ".bgraph";
	std::string intvOutName = graphOutName + ".intv";
	std::string gvzOutName  = graphOutName + ".dot";
	std::string pathOutName = graphOutName + ".path";
	std::string edgeOutName = graphOutName + ".edge";
	std::string euGraphOutName = graphOutName + ".graph";
	CheckEdges(graph.vertices, graph.edges);
  PrintEdges(graph.vertices, graph.edges, edgeOutName);

	graph.PrintIntervalGraph(bGraphOutName, intvOutName);


  PrintGraph(graph.vertices,graph.edges, euGraphOutName);
  GVZPrintBGraph(graph.vertices, graph.edges, gvzOutName);
	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths);

}
