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
#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#ifndef NULL
#define NULL 0
#endif

#ifndef true
#define true 1
#define false 0
#endif

#define VERSION_NUMBER 0
#define RELEASE_NUMBER 7
#define UPDATE_NUMBER 62 

#define MAXLINE 5000 

#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#define inline __inline
extern struct tString_st;
extern struct readSet_st;
extern struct splayTable_st;
extern struct annotation_st;
extern struct roadmap_st;
extern struct insertionMarker_st;
extern struct arc_st;
extern struct node_st;
extern struct graph_st;
extern struct passage_st;
extern struct passageList_st;
extern struct readStart_st;
extern struct preArc_st;
extern struct preNode_st;
extern struct preGraph_st;
extern struct fibheap;
extern struct fibheap_el;
extern struct dfibheap;
extern struct dfibheap_el;
#endif

// Namespace sizes
#include <stdint.h>
typedef int8_t boolean;
typedef int8_t Nucleotide;
typedef uint8_t Descriptor;
typedef int32_t IDnum;
typedef int64_t Coordinate;
typedef int16_t ShortLength;
typedef double Time;
typedef uint8_t Quality;
typedef double Probability;
typedef int8_t Category;

// Atomic word
typedef struct kmer_st Kmer;
typedef int64_t KmerKey;

// Just a sequence string, but with just two bits per character
typedef struct tString_st TightString;

// A simple container when reading files
typedef struct readSet_st ReadSet;
typedef struct sequenceReader_st SequenceReader;

// Hash table structures
typedef struct splayTable_st SplayTable;

// Graph construction structures
typedef struct annotationList_st AnnotationList;
typedef struct annotation_st Annotation;
typedef struct roadmap_st RoadMap;
typedef struct roadMapArray_st RoadMapArray;
typedef struct insertionMarker_st InsertionMarker;

// Pre-Graph elements
typedef struct preArc_st PreArc;
typedef struct preNode_st PreNode;
typedef struct preGraph_st PreGraph;

// Graph elements
typedef struct arc_st Arc;
typedef struct node_st Node;
typedef struct graph_st Graph;
typedef struct shortReadMarker_st ShortReadMarker;
typedef struct passage_st PassageMarker;
typedef struct passageList_st PassageMarkerList;
typedef struct readStart_st ReadStart;
typedef struct gapMarker_st GapMarker;

// Fibonacci heaps used mainly in Tour Bus
typedef struct fibheap FibHeap;
typedef struct fibheap_el FibHeapNode;
typedef struct dfibheap DFibHeap;
typedef struct dfibheap_el DFibHeapNode;

typedef struct nodeList_st NodeList;

#endif
