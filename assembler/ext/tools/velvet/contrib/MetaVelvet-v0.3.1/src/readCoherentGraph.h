/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#ifndef _READCOHERENTGRAPH_H_
#define _READCOHERENTGRAPH_H_

// Original
boolean isUniqueSolexaSubgraph(Node * node, double expCovSubgraph);
int identifyAndSeparateInterRepeats(Graph * argGraph, double * expCovMulti,
				    double repeatNodeCovSD);
void identifyUniqueNodesSubgraph(Graph * graph, int * subgraphMask, 
				 boolean(*isUniqueSubgraph) (Node *, double), 
				 double expCovSubgraph);
void readCoherentSubgraph(Graph * inGraph, double expCovSubgraph, 
			  ReadSet * reads, int * subgraphMask);
// Original

void readCoherentGraph(Graph * graph, boolean(*isUnique) (Node * node),
		       double coverage, ReadSet * reads);

boolean isUniqueBasic(Node * node);

boolean isUniqueSolexa(Node * node);

void setBaseCoverage(double coverage);

void setMultiplicityCutoff(int value);
#endif
