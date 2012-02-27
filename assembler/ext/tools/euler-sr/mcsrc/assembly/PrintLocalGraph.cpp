/***************************************************************************
 * Title:          PrintLocalGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/03/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "compatibility.h"
#include "IntervalGraph.h"
#include "DeBruijnGraph.h"
#include "utils.h"
#include "IntegralTupleStatic.h"


void RecPrintGraph(IntervalGraph &G, ssize_t curVertex, ssize_t depth, 
									 ssize_t maxRevDepth, ssize_t maxDepth, std::ofstream &graphOut);


int main(int argc, char* argv[]) {

	std::string b, g, e, i, p, out;
	ssize_t vertex, forwardBreadth, reverseBreadth;
	
	if (argc != 6) {
		std::cout << "usage: printLocalGraph graphbase gvzOut vertex reverseBreath forwardBreath " << std::endl;
		std::cout << "       Do a dfs search of graph, and print the traversed graph " << std::endl;
		exit(1);
	}
			
	b = argv[1];
	g = b + ".bgraph";
	e = b + ".edge";
	i = b + ".intv";
	p = b + ".path";
	out = argv[2];
	vertex = atosz(argv[3]);
	reverseBreadth = atoi(argv[4]);
	forwardBreadth = atoi(argv[5]);
	std::cout << "v: " << vertex << " b: " << reverseBreadth << " " << forwardBreadth << std::endl;
	IntervalGraph G;
	G.ReadIntervalGraph(g, i, p);
	std::ofstream gvzOut;
	openck(out, gvzOut, std::ios::out);

	G.Unmark();
  gvzOut << "digraph G {"<<std::endl
					 << "\tsize=\"8,8\";"<<std::endl;
	RecPrintGraph(G, vertex, 0, reverseBreadth, forwardBreadth, gvzOut);
	gvzOut << "}" << std::endl;
	
}


void RecPrintGraph(IntervalGraph &G, ssize_t v, ssize_t depth, 
									 ssize_t maxRevDepth, ssize_t maxForDepth, std::ofstream &graphOut) {

	if (G.vertices[v].marked == GraphVertex::Marked) 
		return;

	ssize_t e;
	ssize_t edge;
	G.vertices[v].marked = GraphVertex::Marked;
	if (maxForDepth > 0) {
		for (e = 0; e < G.vertices[v].out.size(); e++ ){ 
			if (G.vertices[v].out[e] != -1) {
				edge = G.vertices[v].out[e];
				if (G.edges[edge].marked != GraphEdge::Marked) {
					graphOut << v << " -> " << G.edges[edge].dest << " [label=\""
									 << edge << " {" << G.edges[edge].index << "} (" << G.edges[edge].length 
									 << " " << G.edges[edge].intervals->size() << ")\"];" << std::endl;
					G.edges[edge].marked = GraphEdge::Marked;
					RecPrintGraph(G, G.edges[edge].dest, depth+1, maxRevDepth, maxForDepth - 1, graphOut);
				}
			}
		}
	}
	if (maxRevDepth > 0) {
		for (e = 0; e < G.vertices[v].in.size(); e++ ){ 
			if (G.vertices[v].in[e] != -1) {
				edge = G.vertices[v].in[e];
				if (G.edges[edge].marked != GraphEdge::Marked) {
					graphOut << G.edges[edge].src << " -> " << v << " [label=\""
									 << edge << " {" << G.edges[edge].index << "} (" << G.edges[edge].length 
									 << " " << G.edges[edge].intervals->size() << ")\"];" << std::endl;
					G.edges[edge].marked = GraphEdge::Marked;
					RecPrintGraph(G, G.edges[edge].src, depth+1, maxRevDepth - 1, maxForDepth, graphOut);
				}
			}
		}
	}
	graphOut << " " << v << " [label=\"" << v <<" (" << G.vertices[v].index
					 << ")\"];" << std::endl;
}
	
				
			

