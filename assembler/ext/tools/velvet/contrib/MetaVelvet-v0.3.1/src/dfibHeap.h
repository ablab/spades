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
#ifndef _DFIBHEAP_H_
#define _DFIBHEAP_H_

DFibHeap *newDFibHeap();

DFibHeapNode *insertNodeIntoDHeap(DFibHeap * heap, Time key, Node * node);

Time replaceKeyInDHeap(DFibHeap * heap, DFibHeapNode * node, Time newKey);

Node *removeNextNodeFromDHeap(DFibHeap * heap);

void destroyDHeap(DFibHeap * heap);

void replaceValueInDHeap(DFibHeapNode * node, Node * newValue);

void *destroyNodeInDHeap(DFibHeapNode * node, DFibHeap * heap);

IDnum getDFibHeapSize(DFibHeap * heap);

Time getKey(DFibHeapNode * node);
#endif
