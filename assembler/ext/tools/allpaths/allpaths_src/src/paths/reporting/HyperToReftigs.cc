/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// HyperToReftigs.  Align assembly unibases or local assemblies to the reference
// sequence.  Then, treating each alignment as defining an interval on the 
// reference sequence, use the graph to piece these alignments together.  The output
// consists of the maximal 'reftig' intervals that can be obtained in this way.

#include "Basevector.h"
#include "Bitvector.h"
#include "Equiv.h"
#include "FetchReadsAmb.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "ParseSet.h"
#include "STLExtensions.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/KmerPath.h"
#include "paths/ReadHBVs.h"
#include "paths/Unipath.h"
#include "paths/reporting/ReftigUtils.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_Bool_OrDefault(USE_CACHE, False);
     CommandArgument_Bool_OrDefault(DISPLAY_GAPS, True);
     CommandArgument_Bool_OrDefault(PRINT_REFTIGS, True);

     // Reference head (relative to <DATA>)
     CommandArgument_String_OrDefault(REF_HEAD, "genome");

     // Optionally generate a dot file. VERBOSITY can be 0, 1, or 2,
     // and it affects how edges are decorated.
     //  0: show base alpha (as in PrettyDOT's edge_labels_base_alpha)
     //  1: specify wich reftig this edge belongs to
     //  2: use integer ids, and show interval of alignment for this edge
     CommandArgument_Bool_OrDefault( DOT, False );
     CommandArgument_Int_OrDefault( VERBOSITY, 1 );
     
     // Default source mode: align <READS>.unibases.k<K> to reference
     CommandArgument_String_OrDefault(READS, "");
     CommandArgument_Int_OrDefault(ORIGIN, 0);
     
     // Alt source mode: align nhoods to reference (if NHOODS != "").
     // If DOT is true, then DOT_OUT must be specified (as full path
     // name of dot output file). Notice: output will be saved as
     // <DOT_OUT>.dot (ie ".dot" will be appended).
     CommandArgument_String_OrDefault( DOT_OUT, "" );
     CommandArgument_String_OrDefault( SUBDIR, "test" );
     CommandArgument_String_OrDefault_Doc( NHOODS, "", "int set (\"{1,2,3}\") or 'all'" );
     
     EndCommandArguments;

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     String full_ref = data_dir + "/" + REF_HEAD;

     String mode = "";
     vec<int> nhoods;
     size_t n_seeds = 0;
     if ( NHOODS != "" ) {
       mode = "alt";
       n_seeds = FirstLineOfFile( sub_dir + "/seeds.ids" ).Int();
       if ( DOT && DOT_OUT == "" ) {
	 cout << "Fatal error, if you select the alt source mode, then you\n"
	      << "must also specify a DOT_OUT.\n" << endl;
	 return 0;
       }
       if ( NHOODS == "all" ) nhoods = vec<int>( n_seeds, vec<int>::IDENTITY );
       else ParseIntSet( NHOODS, nhoods );
     }
     else if ( READS != "" )
       mode = "default";
     else {
       cout << "Neither READS nor NHOODS were selected. Leaving." << endl;
       return 0;
     }
     
     vec< pair<int,ho_interval> > reftigs;

     // Default mode.

     if ( mode == "default" ) 
     {    String unibases = READS + ".unibases.k" + ToString(K);
          String dot_base = run_dir + "/" + READS + ".unibases";
          vec<look_align> ALIGNS;
          GetAlignsFast( K, run_dir + "/" + unibases, full_ref + ".lookup",
			 run_dir + "/" + unibases + ".aligns", ALIGNS, USE_CACHE, run_dir + "/tmp" );
          String strK = "k" + ToString(K);
          String QUERY_HEAD = run_dir + "/" + READS;
          String query_paths_file = QUERY_HEAD + ".paths." + strK;
          String query_paths_rc_file = QUERY_HEAD + ".paths_rc." + strK;
          String query_pathsdb_file = QUERY_HEAD + ".pathsdb." + strK;
          String query_unipaths_file = QUERY_HEAD + ".unipaths." + strK;
          String query_unipathsdb_file = QUERY_HEAD + ".unipathsdb." + strK;
          vecKmerPath paths( query_paths_file ), paths_rc( query_paths_rc_file );
          vecKmerPath unipaths( query_unipaths_file );
          BREAD2( query_pathsdb_file, vec<tagged_rpint>, pathsdb );
          BREAD2( query_unipathsdb_file, vec<tagged_rpint>, unipathsdb );
          digraph G;
          BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, 
               unipathsdb, G );
          HyperKmerPath hkp;
          BuildUnipathAdjacencyHyperKmerPath( K, G, unipaths, hkp );
	  digraph agraph;
          HyperToReftigsCore( K, hkp, ALIGNS, reftigs, &agraph );
	  if ( DOT )
	    GenerateDot( ORIGIN, dot_base, hkp, agraph, ALIGNS,
			 reftigs, VERBOSITY );    }

     // Alt mode.

     else
     {    vec<look_align> all_aligns;
	  digraph agraph;
          int tig_count = 0;
	  
	  // Load selected hypers.
	  vec<HyperBasevector> hypers;
	  vec<bool> select;
	  if ( NHOODS == "all" ) select.resize( n_seeds, True );
	  else {
	    select.resize( n_seeds, False );
	    for (int ii=0; ii<nhoods.isize( ); ii++)
	      select[ nhoods[ii] ] = True;
	  }
	  readHBVs( sub_dir, select, &hypers );

	  if ( NHOODS != "all" ) {
	    vec<HyperBasevector> compacted_hypers;
	    compacted_hypers.reserve( nhoods.size( ) );
	    for (int ii=0; ii<nhoods.isize( ); ii++)
	      compacted_hypers.push_back( hypers[ nhoods[ii] ] );
	    swap( compacted_hypers, hypers );
	  }
		
	  String dot_base = DOT_OUT;
	  
          for (size_t nid=0; nid<nhoods.size( ); nid++) 
          {    int nhood_id = nhoods[nid];

	       // Dir and file names (detect seed directory tree format).
               int nhoodbase = nhood_id / 1000;
               String new_relhead 
                    = ToString( nhoodbase ) + "/" + ToString( nhood_id );
               String old_relhead = GetNestedDirsFromKey( nhood_id );
	       String base_head = sub_dir + "/seed/";
               String new_logfile = base_head + "/" + new_relhead + ".log";
               if ( IsRegularFile( new_logfile ) )
		 base_head += new_relhead + ".hbv";
               else base_head += old_relhead + "/hbv";
               vec<look_align> ALIGNS;
               GetAlignsFast( K, base_head + ".fastb", full_ref + ".lookup",
			      base_head + ".aligns", ALIGNS, USE_CACHE, run_dir + "/tmp" );
	       ForceAssertEq( (size_t)hypers[nid].EdgeObjectCount(),
			      MastervecFileObjectCount( base_head + ".fastb" ) );
               for ( int i = 0; i < ALIGNS.isize( ); i++ )
               {    ALIGNS[i].query_id += tig_count;
                    all_aligns.push_back( ALIGNS[i] );    }
               tig_count += hypers[nid].EdgeObjectCount( );    }
          HyperBasevector h( K, hypers );
          HyperToReftigsCore( K, h, all_aligns, reftigs, &agraph );
	  if ( DOT )
	    GenerateDot( 0, dot_base, h, agraph, all_aligns, reftigs, VERBOSITY );
     }

     // Print reftigs.

     if (PRINT_REFTIGS) {
       vecbitvector amb;
       if ( DISPLAY_GAPS )
	 FetchReadsAmb( amb, full_ref + ".fasta", AMB_EQ_Nn );
       vecbitvector *p_amb = DISPLAY_GAPS ? &amb : 0;
       PrintReftigs( cout, K, ORIGIN, reftigs, p_amb );
     }

     vecbasevector genome( full_ref + ".fastb" );
     vecbitvector cov;
     Mimic( genome, cov );
     for ( int i = 0; i < reftigs.isize( ); i++ )
     {    int t = reftigs[i].first;
          for ( int j = reftigs[i].second.Start( ); 
               j < reftigs[i].second.Stop( ); j++ )
          {    cov[t].Set( j, True );    }    }
     vec<int> reftig_sizes;
     for ( int i = 0; i < reftigs.isize( ); i++ )
          reftig_sizes.push_back( reftigs[i].second.Length( ) );
     Sort(reftig_sizes);
     cout << "total reftigs = " << reftigs.size( ) << "\n";
     cout << "N50 reftig = " << ToStringAddCommas( N50(reftig_sizes) ) << endl;
     cout << "coverage = " << setprecision(3) << 100.0 * Coverage(cov) 
          << "%" << endl;    }
