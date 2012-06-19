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

#include "fib.h"
#include "utility.h"

// Constructor
// Memory allocated
FibHeap *newFibHeap()
{
	FibHeap * heap = fh_makekeyheap();
	if (heap == NULL)
		exitErrorf(EXIT_FAILURE, true, "Can't allocate FibHeap");
	
	return heap;
}

// Add new node into heap with a key, and a pointer to the specified node
FibHeapNode *insertNodeIntoHeap(FibHeap * heap, Coordinate key,
				struct node_st * node)
{
	return fh_insertkey(heap, key, node);
}

// Returns smallest key in heap
Coordinate minKeyOfHeap(FibHeap * heap)
{
	return fh_minkey(heap);
}

// Replaces the key for a given node
Coordinate replaceKeyInHeap(FibHeap * heap, FibHeapNode * node,
			    Coordinate newKey)
{
	return fh_replacekey(heap, node, newKey);
}

// Removes the node with the shortest key, then returns it.
struct node_st *removeNextNodeFromHeap(FibHeap * heap)
{
	return (struct node_st *) fh_extractmin(heap);
}

// Destructor
void destroyHeap(FibHeap * heap)
{
	fh_deleteheap(heap);
}

// Replace the node pointed to by a heap node
void replaceValueInHeap(FibHeapNode * node, Node * newValue)
{
	fh_replacedata(node, newValue);
}

// Remove unwanted node
void destroyNodeInHeap(FibHeapNode * node, FibHeap * heap)
{
	fh_delete(heap, node);
}
