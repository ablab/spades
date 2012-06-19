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

#include  "globals.h"
#include "dfib.h"
#include "utility.h"

// Return number of elements stored in heap
IDnum getDFibHeapSize(DFibHeap * heap)
{
	return dfibheap_getSize(heap);
}

// Constructor
// Memory allocated
DFibHeap *newDFibHeap()
{
	DFibHeap* dheap = dfh_makekeyheap();
	if (dheap == NULL)
		exitErrorf(EXIT_FAILURE, true, "Can't allocate DFibHeap");

	return dheap;
}

// Add new node into heap with a key, and a pointer to the specified node
DFibHeapNode *insertNodeIntoDHeap(DFibHeap * heap, Time key,
				  struct node_st * node)
{
	DFibHeapNode *res;
	res = dfh_insertkey(heap, key, node);

	return res;
}

// Replaces the key for a given node
Time replaceKeyInDHeap(DFibHeap * heap, DFibHeapNode * node, Time newKey)
{
	Time res;
	res = dfh_replacekey(heap, node, newKey);

	return res;
}

// Removes the node with the shortest key, then returns it.
Node *removeNextNodeFromDHeap(DFibHeap * heap)
{
	Node *node;
	node = (Node *) dfh_extractmin(heap);

	return node;
}

// Destructor
void destroyDHeap(DFibHeap * heap)
{
	dfh_deleteheap(heap);
}

// Replace the node pointed to by a heap node
void replaceValueInDHeap(DFibHeapNode * node, Node * newValue)
{
	dfh_replacedata(node, newValue);
}

// Remove unwanted node
void destroyNodeInDHeap(DFibHeapNode * node, DFibHeap * heap)
{
	dfh_delete(heap, node);
}

Time getKey(DFibHeapNode * node)
{
	return dfibheap_el_getKey(node);
}
