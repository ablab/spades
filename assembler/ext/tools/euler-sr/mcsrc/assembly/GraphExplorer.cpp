/***************************************************************************
 * Title:          GraphExplorer.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h" 
#include "MateLibrary.h"
#include "PathLib.h"
#include <iostream>
#include <set>
#include "IntegralTupleStatic.h"
#include "compatibility.h"

 
void PrintEdge(IntervalGraph &g, ssize_t edge) {
	std::cout << "E " << edge << "  src " << g.edges[edge].src << " dest " << g.edges[edge].dest << " length "
						<< g.edges[edge].length << " #intv  " << g.edges[edge].intervals->size() << " bal "
						<< g.edges[edge].balancedEdge << std::endl;
}

void PrintEdgeTraces(IntervalGraph &g, PathTraceList &pathTraces, 
										 TraceMapMatrix &traceMap, ssize_t edge) {
	PrintEdge(g, edge);
	PrintEdgeTraces(pathTraces, traceMap, edge);
}

void PrintVertex(IntervalGraph &g, ssize_t vertex) {
	std::cout << "V " << vertex << std::endl;
	std::cout << " in: " << std::endl;
	if (0 > vertex or vertex >= g.vertices.size()) {
		std::cout << vertex << " is out of bounds." << std::endl;
		return; // 0;
	}
			
	ssize_t i;
	ssize_t edge;
	for (i = 0; i < g.vertices[vertex].in.size(); i++) {
		edge = g.vertices[vertex].in[i];
		if (edge >= 0) {
			std::cout << "  " << edge;
			std::cout << " length " << g.edges[edge].length 
								<< " src "    << g.edges[edge].src 
								<< " #intv "  << g.edges[edge].intervals->size()
								<< " bal "    << g.edges[edge].balancedEdge;
			std::cout << std::endl;
		}
	}
	std::cout << " out: " << std::endl;
	ssize_t o;
	for (o = 0; o < g.vertices[vertex].out.size(); o++)  {
		edge = g.vertices[vertex].out[o];
		if (edge >= 0) {
			std::cout << "  " << edge  << " ";
			std::cout << " length "<< g.edges[edge].length << " " 
								<< " dest "   << g.edges[edge].dest << " "
								<< " #intv " << g.edges[edge].intervals->size() << " " 
								<< " bal "   << g.edges[edge].balancedEdge;
			std::cout << std::endl;
		}
	}
}

void PrintCommands() {
	std::cout <<"   Navigate the graph with the following commands: " << std::endl
						<< "   [e E]  Set the current edge to E" << std::endl
						<< "   [v V]  Set the current vertex to V" << std::endl
						<< "   [f]    Set the current vertex to the dest of the current edge" << std::endl
						<< "   [b]    Set the current vertex to the source of the current edge" << std::endl
						<< "   [i I]    Current edge is the I'th in edge from current vertex" << std::endl
						<< "   [o O]  Current edge is the O'th out edge from current vertex" << std::endl
						<< "   [t T]  Print traces through edge 'T'" << std::endl;
	std::cout << "   [m E t]  Print mate pairs of type 't' that are connected to 'E'. " << std::endl;

}

int main(int argc, char* argv[]) {

	std::string graphName;
	if (argc < 2) {
		std::cout <<" usage: graphExplorer graphName [matesFile]" <<std::endl;
		PrintCommands();
		exit(0);
	}
	graphName = argv[1];
	int argi = 2;
	_INT_ useMates = 0;
	ReadMateList mateList;
	if (argi < argc) {
		std::string mateTableName = argv[argi];
		ReadMateTable(mateTableName, mateList);
		useMates = 1;
	}
		
	IntervalGraph graph;
	int vertexSize;
	ReadIntervalGraph(graphName, graph, vertexSize, !useMates);

	PathBranch pathTree, removedPathTree;
	
	CollectPathTree(graph.paths, graph.pathLengths, pathTree);
	PathTraceList pathTraces;
	PathTreeToPathList(pathTree, pathTraces);
	TraceMapMatrix traceMaps;
	traceMaps.resize(graph.edges.size());
	StoreTraceMaps(pathTraces, traceMaps);

	std::string command, opt;
	ssize_t edge, vertex;
	while(std::cin){
		std::cout << ">";
		std::cout.flush();
		std::cin >> command;
		if (command == "") {
			std::cout << std::endl;
			continue;
		}
		if (command == "t") {
			std::cin >> edge;
			if (0 > edge or edge >= graph.edges.size()) {
				std::cout << "There are only " << graph.edges.size() << " edges." << std::endl;
			continue;
			}
			PrintEdgeTraces(pathTraces, traceMaps, edge);
		}
		else if (command == "e") {
			std::cin >> edge;
			if (0 > edge or edge >= graph.edges.size()) {
				std::cout << "There are only " << graph.edges.size() << " edges." << std::endl;
			continue;
			}
			PrintEdge(graph, edge);
			std::cout << edge << " src: " << std::endl;
			PrintVertex(graph, graph.edges[edge].src);
			std::cout << edge << " dest: " << std::endl;
			PrintVertex(graph, graph.edges[edge].dest);

		}
		else if (command == "v") {
			std::cin >> vertex;
			PrintVertex(graph, vertex);
		}
		else if (command == "f") {
			vertex = graph.edges[edge].dest;
			ssize_t forEdge;
			std::cin >> forEdge;
			ssize_t o;
			if (forEdge < 0) {
				std::cout << " enter a real edge." << std::endl;
				continue;
			}
			for (o = 0; graph.vertices[vertex].out.size() > o; o++ ) {
				if (graph.vertices[vertex].out[o] == forEdge) {
					edge = graph.vertices[vertex].out[o];
					break;
				}
			}
			if (o >= graph.vertices[vertex].out.size()){
				std::cout << "edge: " << forEdge << " not found." << std::endl;
				continue;
			}
			edge = forEdge;
			vertex = graph.edges[edge].dest;
			std::cout << edge << " src: " << std::endl;
			PrintVertex(graph, graph.edges[edge].src);
			std::cout << edge << " dest: " << std::endl;
			PrintVertex(graph, graph.edges[edge].dest);
		}
		else if (command == "b") {
			vertex = graph.edges[edge].src;
			ssize_t backEdge;
			std::cin >> backEdge;
			ssize_t i;
			if (backEdge < 0) {
				std::cout << " enter a real edge." << std::endl;
				continue;
			}
			for (i = 0; i < graph.vertices[vertex].in.size(); i++ ) {
				if (graph.vertices[vertex].in[i] == backEdge) {
					edge = graph.vertices[vertex].in[i];
					break;
				}
			}
			if (i >= graph.vertices[vertex].in.size()){
				std::cout << "edge: " << backEdge << " not found." << std::endl;
				continue;
			}
			vertex = graph.edges[backEdge].src;
			edge = backEdge;
			std::cout << edge << " src: " << std::endl;
			PrintVertex(graph, graph.edges[edge].src);
			std::cout << edge << " dest: " << std::endl;
			PrintVertex(graph, graph.edges[edge].dest);
		}
		else if (command == "i") {
			ssize_t in;
			std::cin >> in;
			if (in < 0 or in >= graph.vertices[vertex].in.size()) {
					std::cout << in << " is out of bounds." << std::endl;
			}
			else {
				if (graph.vertices[vertex].in[in] != -1) {
						edge = graph.vertices[vertex].in[in];
						PrintEdge(graph, edge);
				}
				else {
					std::cout << "No in edge at slot: " << in << std::endl;
				}
			}
			
		}
		else if (command == "o") {
			ssize_t out;
			std::cin >> out;
			if ((out < 0) or (out >= graph.vertices[vertex].out.size())) {
				std::cout << out << " is too large." << std::endl;
			}
			else {
				if (graph.vertices[vertex].out[out] != -1) {
					edge = graph.vertices[vertex].out[out];
					PrintEdge(graph, edge);
				}
				else {
					std::cout << "No out edge exists at position " << out << std::endl;
				}
			}
		}
		else if (command == "m") {
			ssize_t edge;
			ssize_t type = -1;
			std::cin >> edge >> type;
			if (edge > 0 and edge <= graph.edges.size()) {
				MateEdgeMap mateEdges;
				CollectMateEdges(graph, mateList, edge, mateEdges, type);
				MateEdgeMap::iterator edgeIt, endIt;
				endIt = mateEdges.end();
				std::cout << std::endl;
				for (edgeIt = mateEdges.begin(); edgeIt != endIt; ++edgeIt) {
					std::cout << (*edgeIt).first << " (" << (*edgeIt).second.count << ") ";
				}
				std::cout << std::endl;
			}
		}
		else PrintCommands();
	}
	std::cout << std::endl;
}

