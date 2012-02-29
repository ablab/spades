///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Fastavector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "VecUtilities.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/SaveScaffoldGraph.h"
#include "util/RunCommand.h"
// MakeDepend: dependency ScaffoldGraphToGnuplot

/**
 * SaveScaffoldGraph
 */
void SaveScaffoldGraph( const String &HEAD,
			const vec<superb> &supers,
			const digraphE<CLinkBundle> &graph,
			const vec<fastavector> &contigs,
			ostream *log )
{
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;
  
  SaveScaffoldAssembly( HEAD, supers, contigs, log );

  out << Date( ) << ": saving graph" << endl;
  String graph_file = HEAD + ".graph";
  int fd = OpenForWrite( graph_file );
  BinaryWrite( fd, graph );

  String theCommand = "ScaffoldGraphToGnuplot HEAD=" + HEAD;
  RunCommandWithLog( theCommand, "/dev/null" );
  
  out << Date( ) << ": SaveScaffoldGraph done" << endl;

}

/**
 * SaveScaffoldAssembly
 */
void SaveScaffoldAssembly( const String &HEAD,
			   const vec<superb> &supers,
			   const vec<fastavector> &contigs,
			   ostream *log,
			   bool save_fastb )
{
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;
  
  out << Date( ) << ": saving supers" << endl;
  String supers_file = HEAD + ".superb";
  WriteSuperbs( supers_file, supers );
  
  out << Date( ) << ": saving contigs fasta" << endl;
  String contigs_file = HEAD + ".contigs.fasta";
  Ofstream( cg_out, contigs_file );
  for (size_t ii=0; ii<contigs.size( ); ii++)
    contigs[ii].Print( cg_out, "contig_" + ToString( ii ) );
  cg_out.close();
  
  if ( save_fastb ) {
    out << Date( ) << ": saving contigs fastb" << endl;
    String fastb_file = HEAD + ".contigs.fastb";
    vecbasevector as_vecbvec;
    FetchReads( as_vecbvec, 0, contigs_file );
    as_vecbvec.WriteAll( fastb_file );
    
    out << Date( ) << ": saving contigs fastamb" << endl;
    String fastamb_file = HEAD + ".contigs.fastamb";
    vecbitvector amb;
    FetchReadsAmb( amb, contigs_file );
    amb.WriteAll( fastamb_file );
  }
  
  out << Date( ) << ": saving assembly fasta" << endl;
  String assembly_file = HEAD + ".assembly.fasta";
  WriteScaffoldedFasta( assembly_file, contigs, supers );

  out << Date( ) << ": SaveScaffoldAssembly done" << endl;

}

