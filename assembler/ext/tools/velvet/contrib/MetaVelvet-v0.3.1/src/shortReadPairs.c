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
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "globals.h"
#include "graph.h"
#include "concatenatedGraph.h"
#include "recycleBin.h"
#include "locallyCorrectedGraph.h"
#include "passageMarker.h"
#include "readSet.h"
#include "utility.h"
#include "scaffold.h"

#define BLOCK_SIZE  100000
#define LN2 1.4
#define BACKTRACK_CUTOFF 100

typedef struct miniConnection_st MiniConnection;

struct miniConnection_st {
	Coordinate distance;
	double variance;
	Connection *frontReference;
	Connection *backReference;
	NodeList *nodeList;
};


// Global pointers
static Graph *graph;
static NodeList *markedNodes;
static RecycleBin *nodeListMemory = NULL;
static MiniConnection *localScaffold = NULL;

static NodeList *allocateNodeList()
{
	if (nodeListMemory == NULL)
		nodeListMemory =
		    newRecycleBin(sizeof(NodeList), BLOCK_SIZE);

	return allocatePointer(nodeListMemory);
}

static void deallocateNodeList(NodeList * nodeList)
{
	deallocatePointer(nodeListMemory, nodeList);
}

static NodeList *recordNode(Node * node)
{
	NodeList *nodeList = allocateNodeList();
	nodeList->node = node;
	nodeList->next = markedNodes;
	nodeList->previous = NULL;

	if (markedNodes != NULL)
		markedNodes->previous = nodeList;

	markedNodes = nodeList;

	return nodeList;
}

static void destroyNodeList(NodeList * nodeList)
{
	//printf("Destroy NL  %p > %p > %p\n", nodeList->previous, nodeList, nodeList->next);

	if (nodeList->previous != NULL)
		nodeList->previous->next = nodeList->next;
	else
		markedNodes = nodeList->next;

	if (nodeList->next != NULL)
		nodeList->next->previous = nodeList->previous;

	nodeList->previous = nodeList->next = NULL;

	deallocateNodeList(nodeList);
}

static Node *popNodeRecord()
{
	MiniConnection *localConnect;

	NodeList *nodeList = markedNodes;
	Node *node;

	if (markedNodes == NULL)
		return NULL;

	node = nodeList->node;
	markedNodes = nodeList->next;
	if (markedNodes != NULL)
		markedNodes->previous = NULL;

	localConnect =
	    &localScaffold[getNodeID(nodeList->node) + nodeCount(graph)];
	localConnect->nodeList = NULL;

	deallocateNodeList(nodeList);
	return node;
}

void detachImprobablePairs(ReadSet * sequences)
{
	IDnum index, nodeIndex;
	IDnum maxNodeIndex = 2 * nodeCount(graph) + 1;
	ShortReadMarker *nodeArray, *shortMarker;
	Node *node;
	IDnum nodeReadCount;
	IDnum seqID, pairID;
	IDnum *mateReads = sequences->mateReads;
	Category *cats = sequences->categories;

	for (nodeIndex = 0; nodeIndex < maxNodeIndex; nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex - nodeCount(graph));
		if (node == NULL)
			continue;

		nodeArray = getNodeReads(node, graph);
		nodeReadCount = getNodeReadCount(node, graph);

		for (index = 0; index < nodeReadCount; index++) {
			shortMarker =
			    getShortReadMarkerAtIndex(nodeArray, index);

			seqID = getShortReadMarkerID(shortMarker);
			if (mateReads[seqID] == -1)
				continue;

			if (getNodeLength(node) -
			    getShortReadMarkerPosition(shortMarker) >
			    2 * getInsertLength(graph, cats[seqID])) {
				pairID = mateReads[seqID];

				if (pairID != -1) {
					mateReads[seqID] = -1;
					mateReads[pairID] = -1;
				}
			}
		}
	}
}

static void resetMiniConnection(Node * node, MiniConnection * localConnect,
				Coordinate distance, double variance,
				Connection * frontReference,
				Connection * backReference, boolean status)
{
	setSingleNodeStatus(node, status);
	localConnect->distance = distance;
	localConnect->variance = variance;
	localConnect->frontReference = frontReference;
	localConnect->backReference = backReference;
	localConnect->nodeList = recordNode(node);
}

static void setEmptyMiniConnection(Node * node)
{
	MiniConnection *localConnect =
	    &localScaffold[getNodeID(node) + nodeCount(graph)];
	localConnect->distance = 0;
	localConnect->variance = 1;
	localConnect->frontReference = NULL;
	localConnect->backReference = NULL;
	localConnect->nodeList = recordNode(node);
	setSingleNodeStatus(node, true);
}

static void readjustMiniConnection(Node * node,
				   MiniConnection * localConnect,
				   Coordinate distance,
				   Coordinate min_distance,
				   double variance,
				   Connection * frontReference,
				   Connection * backReference)
{

	localConnect->distance =
	    (variance * localConnect->distance +
	     distance * localConnect->variance) / (variance +
						   localConnect->variance);
	localConnect->variance =
	    (variance *
	     localConnect->variance) / (variance + localConnect->variance);

	if (frontReference != NULL)
		localConnect->frontReference = frontReference;
	if (backReference != NULL)
		localConnect->backReference = backReference;

	if (localConnect->distance > min_distance)
		setSingleNodeStatus(node, 1);
	else
		setSingleNodeStatus(node, -1);
}

static void integrateDerivativeDistances(Connection * connect,
					 Coordinate min_distance,
					 boolean direction)
{
	Node *reference = getConnectionDestination(connect);
	Node *destination;
	IDnum destinationID;
	Coordinate distance, baseDistance;
	double variance, baseVariance;
	Connection *connect2;
	MiniConnection *localConnect;

	// debug 
	IDnum counter = 0;

	if (!getUniqueness(reference))
		return;

	//printf("Opposite node %li length %li at %li ± %f\n", getNodeID(reference), getNodeLength(reference), getConnectionDistance(connect), getConnectionVariance(connect));

	baseDistance = getConnectionDistance(connect);
	baseVariance = getConnectionVariance(connect);

	for (connect2 = getConnection(reference);
	     connect2 != NULL; connect2 = getNextConnection(connect2)) {
		// Avoid null derivative
		if (connect2 == getTwinConnection(connect))
			continue;

		destination = getConnectionDestination(connect2);

		// Beware of directionality
		if (!direction)
			destination = getTwinNode(destination);

		// Derivate values
		destinationID = getNodeID(destination);
		// Beware of directionality (bis)
		if (direction)
			distance = baseDistance - getConnectionDistance(connect2);
		else
			distance = getConnectionDistance(connect2) - baseDistance;
		variance = getConnectionVariance(connect2) + baseVariance;
		localConnect =
		    &localScaffold[destinationID + nodeCount(graph)];

		// Avoid over-projection
		if (distance < min_distance) {
			//printf("Node %li not at distance %li± %f (min %li)\n", destinationID, distance, variance, min_distance);
			continue;
		}

		counter++;

		if (getNodeStatus(destination)) {
			readjustMiniConnection(destination, localConnect,
					       distance, min_distance,
					       variance, NULL, NULL);
		} else
			resetMiniConnection(destination, localConnect,
					    distance, variance, NULL, NULL,
					    true);

		//printf("Node %li now at distance %li\n", destinationID, localConnect->distance);
	}

	//printf("%li secondary distances added\n", counter);
}

static void markInterestingNodes(Node * node)
{
	Connection *connect;
	Node *destination;
	MiniConnection *localConnect;
	Coordinate min_distance =
	    getNodeLength(node) / 2 - BACKTRACK_CUTOFF;

	// Mark own node
	setEmptyMiniConnection(node);

	// Loop thru primary scaffold
	for (connect = getConnection(node); connect != NULL;
	     connect = getNextConnection(connect)) {
		destination = getTwinNode(getConnectionDestination(connect));

		localConnect =
		    &localScaffold[getNodeID(destination) +
				   nodeCount(graph)];

		if (getNodeStatus(destination)) {
			readjustMiniConnection(destination, localConnect,
					       getConnectionDistance(connect),
					       min_distance,
					       getConnectionVariance(connect), connect,
					       NULL);
			localConnect->backReference = NULL;
		} else {
			resetMiniConnection(destination, localConnect,
					    getConnectionDistance(connect),
					    getConnectionVariance(connect), connect,
					    NULL, true);
		}

		integrateDerivativeDistances(connect, min_distance, true);
	}

	// Loop thru twin's primary scaffold
	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect)) {
		destination = getConnectionDestination(connect);
		localConnect =
		    &localScaffold[getNodeID(destination) +
				   nodeCount(graph)];

		if (getNodeStatus(destination))
			readjustMiniConnection(destination, localConnect,
					       -getConnectionDistance(connect),
					       min_distance,
					       getConnectionVariance(connect), NULL,
					       connect);
		else
			resetMiniConnection(destination, localConnect,
					    -getConnectionDistance(connect),
					    getConnectionVariance(connect), NULL,
					    connect, -1);

		integrateDerivativeDistances(connect, min_distance, false);
	}
}

void unmarkNode(Node * node, MiniConnection * localConnect)
{
	if (localConnect->frontReference != NULL
	    || localConnect->backReference != NULL) {
		if (getNodeStatus(node) > 0)
			setSingleNodeStatus(node, 10);
		else
			setSingleNodeStatus(node, -10);
	} else {
		setSingleNodeStatus(node, false);
		destroyNodeList(localConnect->nodeList);
		localConnect->frontReference = NULL;
		localConnect->backReference = NULL;
		localConnect->nodeList = NULL;
	}
}

void handicapNode(Node * node)
{
	if (getNodeStatus(node) > 0)
		setSingleNodeStatus(node, 10);
	else
		setSingleNodeStatus(node, -10);
}

static void absorbExtension(Node * node, Node * extension)
{
	Arc *arc;

	appendNodeGaps(node, extension, graph);
	appendDescriptors(node, extension);

	// Destroy old nodes    
	while (getArc(node) != NULL)
		destroyArc(getArc(node), graph);

	// Create new
	for (arc = getArc(extension); arc != NULL; arc = getNextArc(arc))
		createAnalogousArc(node, getDestination(arc), arc, graph);
}

NodeList *getMarkedNodeList()
{
	return markedNodes;
}

static void absorbExtensionInScaffold(Node * node, Node * source)
{
	IDnum nodeID = getNodeID(node);
	IDnum sourceID = getNodeID(source);
	IDnum sourceIndex = sourceID + nodeCount(graph);
	Node *twinSource = getTwinNode(source);
	IDnum twinSourceIndex = getNodeID(twinSource) + nodeCount(graph);
	Connection *connect, *original;
	Node *destination;
	IDnum destinationID;
	Coordinate distance_shift =
	    (getNodeLength(node) - getNodeLength(source)) / 2;
	Coordinate min_distance =
	    getNodeLength(node) / 2 - BACKTRACK_CUTOFF;
	MiniConnection *localConnect;
	Coordinate distance;
	double variance;
	IDnum direct_count;
	IDnum paired_count;

	while ((connect = getConnection(source))) {
		destination = getTwinNode(getConnectionDestination(connect));

		if (destination == getTwinNode(node)) {
			localConnect = &localScaffold[twinSourceIndex];
			localConnect->frontReference = NULL;
			unmarkNode(twinSource, localConnect);
			destroyConnection(connect, sourceID);
			continue;
		}
		if (destination == node) {
			localConnect = &localScaffold[sourceIndex];
			localConnect->backReference = NULL;
			unmarkNode(source, localConnect);
			destroyConnection(connect, sourceID);
			continue;
		}

		destinationID = getNodeID(destination);
		localConnect =
		    &localScaffold[destinationID + nodeCount(graph)];
		incrementConnectionDistance(connect, distance_shift);
		distance = getConnectionDistance(connect);
		variance = getConnectionVariance(connect);
		direct_count = getConnectionDirectCount(connect);
		paired_count = getConnectionPairedCount(connect);

		if (getNodeStatus(destination)) {
			readjustMiniConnection(destination, localConnect,
					       distance, min_distance,
					       variance, NULL, NULL);
			if ((original = localConnect->frontReference))
				readjustConnection(original, distance,
						   variance, direct_count,
						   paired_count);
			else
				localConnect->frontReference =
				    createNewConnection(nodeID,
							-destinationID,
							direct_count,
							paired_count,
							distance,
							variance);
		} else
			resetMiniConnection(destination, localConnect,
					    distance, variance,
					    createNewConnection(nodeID,
								-destinationID,
								direct_count,
								paired_count,
								distance,
								variance),
					    NULL, true);

		integrateDerivativeDistances(connect, min_distance, true);

		destroyConnection(connect, sourceID);
	}

	// Loop thru twin's primary scaffold
	while ((connect = getConnection(getTwinNode(source)))) {
		destination = getConnectionDestination(connect);

		if (destination == node) {
			localConnect = &localScaffold[sourceIndex];
			localConnect->frontReference = NULL;
			unmarkNode(source, localConnect);
			destroyConnection(connect, -sourceID);
			continue;
		}
		if (destination == getTwinNode(node)) {
			localConnect = &localScaffold[twinSourceIndex];
			localConnect->backReference = NULL;
			unmarkNode(twinSource, localConnect);
			destroyConnection(connect, -sourceID);
			continue;
		}

		destinationID = getNodeID(destination);

		localConnect =
		    &localScaffold[destinationID + nodeCount(graph)];
		incrementConnectionDistance(connect, -distance_shift);
		distance = getConnectionDistance(connect);
		variance = getConnectionVariance(connect);
		direct_count = getConnectionDirectCount(connect);
		paired_count = getConnectionPairedCount(connect);

		if (distance > min_distance && getNodeStatus(destination) < 0) {
			readjustMiniConnection(destination, localConnect,
					       -distance, min_distance,
					       variance, NULL, NULL);
			if ((original = localConnect->backReference))
				readjustConnection(original, distance,
						   variance, direct_count,
						   paired_count);
		} else if (getNodeStatus(destination) < 0) {
			if ((original = localConnect->backReference)) {
				destroyConnection(original, -nodeID);
				localConnect->backReference = NULL;
			}
			unmarkNode(destination, localConnect);
		} else if (getNodeStatus(destination) > 0) {
			if ((original = localConnect->frontReference)) {
				destroyConnection(original, nodeID);
				localConnect->frontReference = NULL;
			}
			unmarkNode(destination, localConnect);
		} else if (distance > min_distance) {
			resetMiniConnection(destination, localConnect,
					    -distance, variance, NULL,
					    createNewConnection(-nodeID,
								destinationID,
								direct_count,
								paired_count,
								distance,
								variance),
					    -1);
			integrateDerivativeDistances(connect, min_distance, true);
		}

		destroyConnection(connect, -sourceID);
	}
}

static void recenterNode(Node * node, Coordinate oldLength)
{
	IDnum nodeID = getNodeID(node);
	Connection *connect, *next;
	Coordinate distance_shift = (getNodeLength(node) - oldLength) / 2;
	Coordinate min_distance =
	    getNodeLength(node) / 2 - BACKTRACK_CUTOFF;
	MiniConnection *localConnect;

	//puts("Recentering node");

	for (connect = getConnection(node); connect != NULL;
	     connect = next) {
		next = getNextConnection(connect);
		incrementConnectionDistance(connect, -distance_shift);

		if (getConnectionDistance(connect) < min_distance) {
			//printf("Unrecording %li\n",
			//       -getNodeID(getConnectionDestination(connect)));
			localConnect =
			    &localScaffold[-getNodeID(getConnectionDestination(connect))
					   + nodeCount(graph)];
			localConnect->frontReference = NULL;
			unmarkNode(getTwinNode(getConnectionDestination(connect)),
				   localConnect);
			destroyConnection(connect, nodeID);
		} else if (getTwinConnection(connect) != NULL)
			incrementConnectionDistance(getTwinConnection(connect), -distance_shift);
	}

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = next) {
		next = getNextConnection(connect);
		incrementConnectionDistance(connect, distance_shift);

		if (getTwinConnection(connect) != NULL)
			incrementConnectionDistance(getTwinConnection(connect), distance_shift);
	}
}

static void recenterLocalScaffold(Node * node, Coordinate oldLength)
{
	MiniConnection *localConnect;
	Coordinate distance_shift = (getNodeLength(node) - oldLength) / 2;
	Coordinate min_distance =
	    getNodeLength(node) / 2 - BACKTRACK_CUTOFF;
	NodeList *nodeList, *next;
	IDnum node2ID;
	Node *node2;

	for (nodeList = markedNodes; nodeList != NULL; nodeList = next) {
		next = nodeList->next;

		node2 = nodeList->node;

		if (node2 == node) {
			setSingleNodeStatus(node2, 1);
			continue;
		}

		node2ID = getNodeID(node2);
		localConnect = &localScaffold[node2ID + nodeCount(graph)];
		localConnect->distance -= distance_shift;

		if (localConnect->distance < min_distance
		    && localConnect->backReference == NULL
		    && localConnect->frontReference == NULL)
			unmarkNode(node2, localConnect);
		else if (getNodeStatus(node2) > 0)
			setSingleNodeStatus(node2, 1);
		else if (getNodeStatus(node2) < 0)
			setSingleNodeStatus(node2, -1);
	}
}

static void adjustShortReads(Node * target, Node * source)
{
	ShortReadMarker *targetArray, *marker;
	IDnum targetLength, index;
	Coordinate position, nodeLength;
	
	if (!readStartsAreActivated(graph))
		return;

	targetArray = getNodeReads(getTwinNode(target), graph);
	targetLength = getNodeReadCount(getTwinNode(target), graph);

	nodeLength = getNodeLength(source);

	for (index = 0; index < targetLength; index++) {
		marker = getShortReadMarkerAtIndex(targetArray, index);
		position = getShortReadMarkerPosition(marker);
		position += nodeLength;
		setShortReadMarkerPosition(marker, position);
	}
}

static void adjustLongReads(Node * target, Node * source)
{
	PassageMarker *marker;
	Coordinate nodeLength = getNodeLength(source);

	for (marker = getMarker(source); marker != NULL;
	     marker = getNextInNode(marker))
		incrementFinishOffset(marker, nodeLength);
}

static boolean goesToNode(PassageMarker * marker, Node * node)
{
	PassageMarker *current;

	for (current = marker; current != NULL;
	     current = getNextInSequence(current))
		if (getNode(current) == node)
			return true;

	return false;
}

static boolean comesFromNode(PassageMarker * marker, Node * node)
{
	Node *target = getTwinNode(node);
	PassageMarker *current;

	for (current = getTwinMarker(marker); current != NULL;
	     current = getNextInSequence(current))
		if (getNode(current) == target)
			return true;

	return false;
}

static void reconnectPassageMarker(PassageMarker * marker, Node * node,
				   PassageMarker ** ptr)
{
	PassageMarker *current;
	PassageMarker *next = getNextInSequence(marker);
	PassageMarker *tmpMarker;

	for (current = marker; getNode(current) != node;
	     current = getPreviousInSequence(current));

	setPreviousInSequence(current, next);
	concatenatePassageMarkers(current, marker);

	for (; marker != current; marker = tmpMarker) {
		tmpMarker = getPreviousInSequence(marker);
		if (*ptr == marker || *ptr == getTwinMarker(marker))
			*ptr = getNextInNode(*ptr);
		setNextInSequence(marker, NULL);
		setPreviousInSequence(NULL, marker);
		destroyPassageMarker(marker);
	}
}

static void concatenateLongReads(Node * node, Node * candidate,
				 Graph * graph)
{
	PassageMarker *marker, *tmpMarker;

	// Passage marker management in node:
	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		if (!goesToNode(marker, candidate))
			incrementFinishOffset(marker,
					      getNodeLength(candidate));
	}

	// Swapping new born passageMarkers from candidate to node
	for (marker = getMarker(candidate); marker != NULL;
	     marker = tmpMarker) {
		tmpMarker = getNextInNode(marker);

		if (!comesFromNode(marker, node)) {
			extractPassageMarker(marker);
			transposePassageMarker(marker, node);
			incrementFinishOffset(getTwinMarker(marker),
					      getNodeLength(node));
		} else {
			reconnectPassageMarker(marker, node, &tmpMarker);
		}
	}
}

static void adjustShortReadsByLength(Node * target, Coordinate nodeLength)
{
	ShortReadMarker *targetArray, *marker;
	IDnum targetLength, index;
	Coordinate position;

	if (!readStartsAreActivated(graph))
		return;

	targetArray = getNodeReads(getTwinNode(target), graph);
	targetLength = getNodeReadCount(getTwinNode(target), graph);

	for (index = 0; index < targetLength; index++) {
		marker = getShortReadMarkerAtIndex(targetArray, index);
		position = getShortReadMarkerPosition(marker);
		position += nodeLength;
		setShortReadMarkerPosition(marker, position);
	}
}

static boolean abs_bool(boolean val)
{
	return val >= 0 ? val : -val;
}

static IDnum abs_ID(IDnum val)
{
	return val >= 0 ? val : -val;
}

static NodeList *pathIsClear(Node * node, Node * oppositeNode,
			     Coordinate distance)
{
	Arc *arc;
	Node *candidate, *dest, *current;
	Coordinate extension_distance = 0;
	boolean maxRepeat = 1;
	Node *repeatEntrance = NULL;
	IDnum counter = 0;
	NodeList *path = NULL;
	NodeList *tail = path;

	setSingleNodeStatus(node, 2);

	current = node;
	while (true) {

		//////////////////////////////////
		//  Selecting destination       //
		//////////////////////////////////
		candidate = NULL;

		// First round for priority nodes
		for (arc = getArc(current); arc != NULL;
		     arc = getNextArc(arc)) {
			dest = getDestination(arc);

			if (dest == node || dest == getTwinNode(node))
				continue;

			if (getNodeStatus(dest) <= 0)
				continue;

			if (candidate == NULL
			    || getNodeStatus(candidate) >
			    getNodeStatus(dest)
			    || (getNodeStatus(candidate) ==
				getNodeStatus(dest)
				&& extension_distance >
				localScaffold[getNodeID(dest) +
					      nodeCount(graph)].
				distance - getNodeLength(dest) / 2)) {
				extension_distance =
				    localScaffold[getNodeID(dest) +
						  nodeCount(graph)].
				    distance - getNodeLength(dest) / 2;
				candidate = dest;
			}
		}

		if (candidate != NULL && repeatEntrance) {
			for (arc = getArc(node); arc != NULL;
			     arc = getNextArc(arc)) {
				dest = getDestination(arc);
				if (dest != candidate
				    && getNodeStatus(dest)) {
					break;
				}
			}
		}
		// In case of failure   
		if (candidate == NULL) {
			for (arc = getArc(current); arc != NULL;
			     arc = getNextArc(arc)) {
				dest = getDestination(arc);

				if (getNodeStatus(dest) == 0)
					continue;

				if (dest == node
				    || dest == getTwinNode(node))
					continue;

				if (candidate == NULL
				    || getNodeStatus(candidate) <
				    getNodeStatus(dest)
				    || (getNodeStatus(candidate) ==
					getNodeStatus(dest)
					&& extension_distance <
					localScaffold[getNodeID(dest) +
						      nodeCount(graph)].
					distance -
					getNodeLength(dest) / 2)) {
					extension_distance =
					    localScaffold[getNodeID(dest) +
							  nodeCount
							  (graph)].
					    distance -
					    getNodeLength(dest) / 2;
					candidate = dest;
				}
			}
		}
		if (candidate == NULL) {
			while (path) {
				tail = path->next;
				deallocateNodeList(path);
				path = tail;
			}
			return false;
		}
		// Loop detection
		if (candidate == repeatEntrance
		    && abs_bool(getNodeStatus(candidate)) ==
		    maxRepeat + 1) {
			while (path) {
				tail = path->next;
				deallocateNodeList(path);
				path = tail;
			}
			return false;
		} else if (abs_bool(getNodeStatus(candidate)) > maxRepeat) {
			maxRepeat = abs_bool(getNodeStatus(candidate));
			repeatEntrance = candidate;
		} else if (abs_bool(getNodeStatus(candidate)) == 1) {
			maxRepeat = 1;
			repeatEntrance = NULL;
		}

		if (getNodeStatus(candidate) > 0)
			setSingleNodeStatus(candidate,
					    getNodeStatus(candidate) + 1);
		else
			setSingleNodeStatus(candidate,
					    getNodeStatus(candidate) - 1);


		// DEBUG 
		if (abs_bool(getNodeStatus(candidate)) > 100
		    || counter > nodeCount(graph)) {
			while (path) {
				tail = path->next;
				deallocateNodeList(path);
				path = tail;
			}
			return false;
		}

		// Missassembly detection
		if (getUniqueness(candidate) && oppositeNode
		    && candidate != oppositeNode
		    && extension_distance > distance) {
			while (path) {
				tail = path->next;
				deallocateNodeList(path);
				path = tail;
			}
			return false;
		}

		if (path == NULL) {
			path = allocateNodeList();
			path->next = NULL;
			path->node = candidate;
			tail = path;
		} else {
			tail->next = allocateNodeList();
			tail = tail->next;
			tail->node = candidate;
			tail->next = NULL;
		}

		if (getUniqueness(candidate))
			return path;

		current = candidate;
	}
}

static boolean pushNeighbours(Node * node, Node * oppositeNode,
			      Coordinate distance, boolean force_jumps)
{
	Node *candidate;
	Node *lastCandidate = NULL;
	Coordinate oldLength = getNodeLength(node);
	Category cat;
	MiniConnection *localConnect;
	NodeList *path, *tmp;

	if ((path = pathIsClear(node, oppositeNode, distance))) {
		while (path) {
			candidate = path->node;
			tmp = path->next;
			deallocateNodeList(path);
			path = tmp;

			///////////////////////////////////////
			//  Stepping forward to destination  //
			///////////////////////////////////////

			if (getUniqueness(candidate)) {
				concatenateReadStarts(node, candidate, graph);
				concatenateLongReads(node, candidate, graph);
				absorbExtension(node, candidate);

				// Scaffold changes
				recenterNode(node, oldLength);
				recenterLocalScaffold(node, oldLength);
				absorbExtensionInScaffold(node, candidate);

				// Read coverage
				for (cat = 0; cat < CATEGORIES; cat++) {
					incrementVirtualCoverage(node, cat,
								 getVirtualCoverage
								 (candidate,
								  cat));
					incrementOriginalVirtualCoverage
					    (node, cat,
					     getOriginalVirtualCoverage
					     (candidate, cat));
				}

				if (getNodeStatus(candidate)) {
					localConnect =
					    &localScaffold[getNodeID
							   (candidate) +
							   nodeCount
							   (graph)];
					if (localConnect->frontReference) {
						destroyConnection
						    (localConnect->
						     frontReference,
						     getNodeID(node));
						localConnect->
						    frontReference = NULL;
					}
					if (localConnect->backReference) {
						destroyConnection
						    (localConnect->
						     backReference,
						     -getNodeID(node));
						localConnect->
						    backReference = NULL;
					}
					unmarkNode(candidate,
						   localConnect);
				}
				if (getNodeStatus(getTwinNode(candidate))) {
					localConnect =
					    &localScaffold[-getNodeID
							   (candidate) +
							   nodeCount
							   (graph)];
					if (localConnect->frontReference) {
						destroyConnection
						    (localConnect->
						     frontReference,
						     getNodeID(node));
						localConnect->
						    frontReference = NULL;
					}
					if (localConnect->backReference) {
						destroyConnection
						    (localConnect->
						     backReference,
						     -getNodeID(node));
						localConnect->
						    backReference = NULL;
					}
					unmarkNode(getTwinNode(candidate),
						   localConnect);
				}
				// Original
				printf("Pebble Concatenated Node %d -- ", 
				       getNodeID(node));
				printf("Node %d\n", getNodeID(candidate));
				// Original

				destroyNode(candidate, graph);
				return true;
			} else {
				adjustShortReads(node, candidate);
				adjustLongReads(node, candidate);
				absorbExtension(node, candidate);
				lastCandidate = candidate;
			}
		}
	}

	if (force_jumps && oppositeNode
	    && abs_ID(getNodeID(oppositeNode)) < abs_ID(getNodeID(node))) {
		distance -= getNodeLength(node) / 2;
		distance -= getNodeLength(oppositeNode) / 2;
		if (distance > 10) {
			adjustShortReadsByLength(node, distance);
			appendGap(node, distance, graph);
		} else {
			adjustShortReadsByLength(node, 10);
			appendGap(node, 10, graph);
		}

		concatenateReadStarts(node, oppositeNode, graph);
		concatenateLongReads(node, oppositeNode, graph);
		absorbExtension(node, oppositeNode);

		// Scaffold changes
		recenterNode(node, oldLength);
		recenterLocalScaffold(node, oldLength);
		absorbExtensionInScaffold(node, oppositeNode);

		// Read coverage
		for (cat = 0; cat < CATEGORIES; cat++)
			incrementVirtualCoverage(node, cat,
						 getVirtualCoverage
						 (oppositeNode, cat));

		if (getNodeStatus(oppositeNode)) {
			localConnect =
			    &localScaffold[getNodeID(oppositeNode) +
					   nodeCount(graph)];
			if (localConnect->frontReference) {
				destroyConnection(localConnect->
						  frontReference,
						  getNodeID(node));
				localConnect->frontReference = NULL;
			}
			if (localConnect->backReference) {
				destroyConnection(localConnect->
						  backReference,
						  -getNodeID(node));
				localConnect->backReference = NULL;
			}
			unmarkNode(oppositeNode, localConnect);
		}
		if (getNodeStatus(getTwinNode(oppositeNode))) {
			localConnect =
			    &localScaffold[-getNodeID(oppositeNode) +
					   nodeCount(graph)];
			if (localConnect->frontReference) {
				destroyConnection(localConnect->
						  frontReference,
						  getNodeID(node));
				localConnect->frontReference = NULL;
			}
			if (localConnect->backReference) {
				destroyConnection(localConnect->
						  backReference,
						  -getNodeID(node));
				localConnect->backReference = NULL;
			}
			unmarkNode(getTwinNode(oppositeNode),
				   localConnect);
		}
		// Original
		printf("Pebble Scaffolded Node %d -- Node %d\n", 
		       getNodeID(node), getNodeID(oppositeNode));
		// Original

		destroyNode(oppositeNode, graph);
	}

	return false;
}

static void unmarkInterestingNodes()
{
	Node *node;
	MiniConnection *localConnect;

	while ((node = popNodeRecord())) {
		setSingleNodeStatus(node, false);
		localConnect =
		    &localScaffold[getNodeID(node) + nodeCount(graph)];
		localConnect->frontReference = NULL;
		localConnect->backReference = NULL;
		localConnect->nodeList = NULL;
	}
}

static void findOppositeNode(Node * node, Node ** oppositeNode,
			     Coordinate * distance)
{
	NodeList *nodeList;
	MiniConnection *localConnect;
	Node *node2;
	IDnum node2ID;

	*oppositeNode = NULL;
	*distance = 0;

	for (nodeList = markedNodes; nodeList != NULL;
	     nodeList = nodeList->next) {
		node2 = nodeList->node;
		node2ID = getNodeID(node2);
		localConnect = &localScaffold[node2ID + nodeCount(graph)];

		if (node2 == node)
			continue;

		if (!getUniqueness(node2))
			continue;

		if (localConnect->distance < 0)
			continue;

		if (*oppositeNode == NULL
		    || *distance > localConnect->distance) {
			*oppositeNode = node2;
			*distance = localConnect->distance;
		}
	}
}

static boolean expandLongNode(Node * node, boolean force_jumps)
{
	boolean hit = true;
	boolean modified = false;
	Node *oppositeNode;
	Coordinate distance = 0;

	markInterestingNodes(node);

	while (hit) {
		correctGraphLocally(node);
		findOppositeNode(node, &oppositeNode, &distance);
		hit = pushNeighbours(node, oppositeNode, 
				     distance, force_jumps);
		modified = modified || hit;
	}

	unmarkInterestingNodes();

	return modified;
}

static boolean expandLongNodes(boolean force_jumps)
{
	IDnum nodeID;
	Node *node;
	boolean modified = false;

	for (nodeID = 1; nodeID <= nodeCount(graph); nodeID++) {
		node = getNodeInGraph(graph, nodeID);

		if (node != NULL && getUniqueness(node)) {
			modified = expandLongNode(node, force_jumps)
			    || modified;
			modified =
			    expandLongNode(getTwinNode(node), force_jumps)
			    || modified;
		}
	}

	return modified;
}

static void cleanMemory()
{
	// Original
	//puts("Cleaning memory");
	// Original

	cleanScaffoldMemory();

	destroyRecycleBin(nodeListMemory);
	nodeListMemory = NULL;

	free(localScaffold);
}

void exploitShortReadPairs(Graph * argGraph, ReadSet * reads,
			   boolean * dubious, boolean force_jumps)
{
	boolean modified = true;

	graph = argGraph;

	if (!readStartsAreActivated(graph))
		return;
	
	// Original
	/*
	puts("Starting pebble resolution...");
	*/
	// Original

	// Prepare graph
	resetNodeStatus(graph);
	prepareGraphForLocalCorrections(graph);

	// Prepare scaffold
	buildScaffold(graph, reads, dubious);

	// Prepare local scaffold 
	localScaffold =
	    callocOrExit(2 * nodeCount(graph) + 1, MiniConnection);

	// Loop until convergence
	while (modified)
		modified = expandLongNodes(force_jumps);

	// Clean up memory
	cleanMemory();
	deactivateLocalCorrectionSettings();

	sortGapMarkers(graph);

	// Original
	/*
	puts("Pebble done.");
	*/
	// Original
}

// Original
static void adjustShortReadsInterRepeat(Node * target, Node * source,
					Graph * argGraph)
{
	ShortReadMarker *targetArray, *marker;
	IDnum targetLength, index;
	Coordinate position, nodeLength;
	Graph * graph = argGraph;
	
	if (!readStartsAreActivated(graph))
		return;

	targetArray = getNodeReads(getTwinNode(target), graph);
	targetLength = getNodeReadCount(getTwinNode(target), graph);

	nodeLength = getNodeLength(source);

	for (index = 0; index < targetLength; index++) {
		marker = getShortReadMarkerAtIndex(targetArray, index);
		position = getShortReadMarkerPosition(marker);
		position += nodeLength;
		setShortReadMarkerPosition(marker, position);
	}
}

static void absorbExtensionInterRepeat(Node * node, Node * extension,
				       Graph * argGraph)
{
	Arc *arc;
	Graph * graph = argGraph;

	appendNodeGaps(node, extension, graph);
	appendDescriptors(node, extension);

	// Destroy old nodes    
	while (getArc(node) != NULL)
		destroyArc(getArc(node), graph);

	// Create new
	for (arc = getArc(extension); arc != NULL; arc = getNextArc(arc))
		createAnalogousArc(node, getDestination(arc), arc, graph);
}

boolean pushNeighboursInterRepeat(Node * node, Node * nodeInterRepeat, 
				  Node * oppositeNode, Graph * argGraph)
{
	Node *candidate;
	Node *lastCandidate = NULL;
	Category cat;
	NodeList *path, *tmp;
	Graph * graph = argGraph;
	
	// Make path (= NodeList of node and oppositeNode)
	path = allocateNodeList();
	path->node = nodeInterRepeat;
	path->next = allocateNodeList();
	path->next->node = oppositeNode;
	path->next->next = NULL;

	// Original
	while (path) {
		candidate = path->node;
		tmp = path->next;
		deallocateNodeList(path);
		path = tmp;
		
		if (getUniqueness(candidate)) {
			concatenateReadStarts(node, candidate, graph);
			concatenateLongReads(node, candidate, graph);
			absorbExtensionInterRepeat(node, candidate, graph);
			// Read coverage
			for (cat = 0; cat < CATEGORIES; cat++) {
				incrementVirtualCoverage(node, cat,
							 getVirtualCoverage
							 (candidate, cat));
				incrementOriginalVirtualCoverage
					(node, cat,
					 getOriginalVirtualCoverage
					 (candidate, cat));
			}

			// Original
			printf("\tConcatenated InNode %d -- OutNode %d\n", 
			       getNodeID(node), getNodeID(candidate));
			// Original
			
			destroyNode(candidate, graph);
			return true;
		} else {
			adjustShortReadsInterRepeat(node, candidate, graph);
			adjustLongReads(node, candidate);
			absorbExtensionInterRepeat(node, candidate, graph);
			lastCandidate = candidate;
		}
	}

	return false;
}
// Original
