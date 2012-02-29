///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SAVE_SCAFFOLD_GRAPH_H
#define SAVE_SCAFFOLD_GRAPH_H

#include "Fastavector.h"
#include "VecUtilities.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/reporting/CLinkBundle.h"

/**
 * SaveScaffoldGraph
 *
 * Save the given scaffold graph, together with various other files:
 *  <HEAD>.superb            supers
 *  <HEAD>.graph             scaffold graph
 *  <HEAD>.graph.txt         gnuplot-parseable scaffold graph
 *  <HEAD>.assembly.fasta    fastavector of supers
 *  <HEAD>.contigs.fasta     fastavector of contigs
 */
void SaveScaffoldGraph( const String &HEAD,
			const vec<superb> &supers,
			const digraphE<CLinkBundle> &graph,
			const vec<fastavector> &contigs,
			ostream *log = 0 );

/**
 * SaveScaffoldAssembly
 *
 * Save scaffolds, contigs.fasta, and assembly.fasta:
 *  <HEAD>.superb               supers
 *  <HEAD>.assembly.fasta       fastavector of supers
 *  <HEAD>.contigs.fasta        fastavector of contigs
 *  <HEAD>.contigs.fast{b,amb}  contigs as fastb, fastamb (if save_fastb = true)
 */
void SaveScaffoldAssembly( const String &HEAD,
			   const vec<superb> &supers,
			   const vec<fastavector> &contigs,
			   ostream *log = 0,
			   bool save_fastb = false );

#endif
