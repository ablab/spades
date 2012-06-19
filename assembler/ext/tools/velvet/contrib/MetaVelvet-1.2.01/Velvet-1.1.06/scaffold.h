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
#ifndef _SSCAFFOLD_H_
#define _SCAFFOLD_H_

typedef struct connection_st Connection;

//General scaffold function
void buildScaffold(Graph * graph,
		   ReadSet * reads,
		   boolean * dubious,
		   boolean * shadows);
Connection *createNewConnection(IDnum nodeID, IDnum node2ID,
				       IDnum direct_count,
				       IDnum paired_count,
				       Coordinate distance,
				       double variance);
void readjustConnection(Connection * connect, Coordinate distance,
			       double variance, IDnum direct_count,
			       IDnum paired_count);
void destroyConnection(Connection * connect, IDnum nodeID);

void cleanScaffoldMemory();

void setUnreliableConnectionCutoff(int val);
void setPairedExpFraction(double x);

// Connection handlers
Connection * getConnection(Node * node);

Node * getConnectionDestination(Connection * connect);
Coordinate getConnectionDistance(Connection * connect);
Connection * getNextConnection(Connection * connect);
Connection * getTwinConnection(Connection * connect);
double getConnectionVariance(Connection * connect);
IDnum getConnectionDirectCount(Connection * connect);
IDnum getConnectionPairedCount(Connection * connect);

void incrementConnectionDistance(Connection * connect, Coordinate increment);
void printConnections(ReadSet * reads, boolean * shadows);
void printScaffold(Graph * graph,
		   ReadSet * reads,
		   boolean * dubious,
		   boolean * shadows);

//
void computeConnections( Graph* argGraph, ReadSet* reads, boolean* dubious, boolean* shadows );

#endif
