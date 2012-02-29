/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// ScaffoldLayout.  Figure out how assembly scaffolds lie on a reference.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Bitvector.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "system/Parallel.h"


int main(int argc, char **argv)
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(ASSEMBLY, 
          "fasta file for assembly, with Ns for gaps");
     CommandArgument_String_Doc(REF_LOOKUP, "lookup file for reference" );
     CommandArgument_String_Doc(OUT_HEAD, 
          "generates OUTHEAD.report, plus temp files that are deleted" );
     CommandArgument_String(TMP_DIR);
     // Minimum unambiguous sequence length used for lookup against reference
     CommandArgument_Int_OrDefault(END_SIZE_MIN, 1000);
     // Maximum sequence length used for lookup against reference
     CommandArgument_Int_OrDefault(END_SIZE_MAX, 5000);
     CommandArgument_Int_OrDefault(VERBOSITY, 0);
     // Reject alignments to reference which imply difference in scaffold
     // length by this fraction (e.g., 0.05 == 5%)
     CommandArgument_Double_OrDefault(SIZE_ERROR_LIMIT, 0.25);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");
     EndCommandArguments;

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Load the reference and assembly.

     cout << Date( ) << ": load assembly" << endl;
     vecbasevector assembly;
     vecbitvector assembly_amb;
     vecString names;
     FetchReads( assembly, names, ASSEMBLY );
     FetchReadsAmb( assembly_amb, ASSEMBLY );

     Ofstream( report, OUT_HEAD + ".report" );

     // Create temp dir for alignments
     Mkdir777(TMP_DIR);

     // For each scaffold, find first and last contig.

     cout << Date( ) << ": get ends" << endl;
     vecbasevector ENDS;
     vec<int> end_offsets(2*assembly.size());
     for ( size_t i = 0, j = 0; i < assembly.size( ); i++ ) {
       ostrstream out;
       basevector e;
       int scaffold_size = assembly[i].isize();
       int start, size;

       // This version finds the first and last non-ambiguous sequences of
       // END_SIZE_MIN <= length <= end_size_max
       // --bruce
       int j, k;
       int end_size_max = min(END_SIZE_MAX, scaffold_size/2);
       for (j = k = 0; k < scaffold_size; k++) {
	 if (assembly_amb[i][k]) {
	   if (k-j >= END_SIZE_MIN) break;
	   else j = k+1;
	 } else if (k-j >= end_size_max) break;
       }
       e.SetToSubOf( assembly[i], j, k-j );
       ENDS.push_back_reserve(e);
       end_offsets[2*i] = j;
       if (VERBOSITY >= 5)
	 report << "scaffold " << names[i] << " size " << scaffold_size << " j0=" << j << " k0=" << k << " o=" << j;
       // find the last non-ambiguous sequence of
       // END_SIZE_MIN <= length <= end_size_max
       for ( j = scaffold_size - 1, k = scaffold_size; j >= 0; j-- ) {
	 if (assembly_amb[i][j]) {
	   if (k-j >= END_SIZE_MIN) {
	     j++;
	     break;
	   } else k = j;
	 } else if (j == 0 || k-j >= end_size_max) break;
       }
       if ( j < 0 ) j = k = 0;
       e.SetToSubOf( assembly[i], j, k-j );
       end_offsets[2*i+1] =  k - scaffold_size;
       ENDS.push_back_reserve(e);
       if (VERBOSITY >= 5 )
	 report << " j1=" << j << " k1=" << k  << " o=" << k-scaffold_size << endl;
     }
     String ends_file = OUT_HEAD + ".ends.fastb";
     ENDS.WriteAll(ends_file);

     // Align the ends.

     cout << Date( ) << ": aligning across " << NUM_THREADS << " processors" << endl;
     int batch = ( ENDS.size( ) + NUM_THREADS - 1 ) / NUM_THREADS;
     vec<String> commands;
     for ( size_t i = 0, start = 0; i < NUM_THREADS; i++ )
     {    int stop = Min( static_cast<size_t>(start + batch), ENDS.size( ) );
          String aligns_file = OUT_HEAD + ".ends.aligns." + ToString(i);
          commands.push_back( "QueryLookupTable" + ARG(K, 12) + ARG(MM, 12) 
			      + ARG(MC, 0.15) + ARG(SEQS, ends_file)
			      + ARG(L, REF_LOOKUP) + ARG(VISUAL, True) 
			      + ARG(PARSEABLE, True) + ARG(TMP_DIR, TMP_DIR)
			      + ARG(SEQS_TO_PROCESS, 
				    "\"[" + ToString(start) + "," + ToString(stop) + ")\"")
			      + " > " + aligns_file );
          start = stop;    }
     RunCommandsInParallel( commands, NUM_THREADS );
     Remove(ends_file);
     cout << Date( ) << ": loading alignments" << endl;
     vec<look_align> aligns;
     for ( size_t  i = 0; i < NUM_THREADS; i++ )
     {    vec<look_align> aligns0;
          String aligns_file = OUT_HEAD + ".ends.aligns." + ToString(i);
          LoadLookAligns( aligns_file, aligns0 );
          aligns.append(aligns0);
          Remove(aligns_file);
     }
     vec< vec<int> > aligns_index( ENDS.size( ) );
     for ( int i = 0; i < aligns.isize( ); i++ )
          aligns_index[ aligns[i].query_id ].push_back(i);

     // Analyze end alignments.

     cout << Date( ) << ": analyze alignments" << endl;
     vec< triple<int,pair<int,int>,String> > answer;
     int good_placements = 0;
     const int infinity = 1000000000;
     for ( size_t i = 0, j = 0; i < assembly.size( ); i++ )
     {    ostrstream out;
          int t = infinity, start = -1, stop = -1;
          int scaffold_size = assembly[i].isize();
	  int j1, j2;
	  
	  j1 = j++;
	  j2 = j++;
          out << "[" << i << ": " << names[i] << ", l=" << scaffold_size << "] ";
          if ( !aligns_index[j1].solo( ) || !aligns_index[j2].solo( ) )
               out << "Not solo.";
          else
          {    const look_align& la1 = aligns[ aligns_index[j1][0] ];
               const look_align& la2 = aligns[ aligns_index[j2][0] ];
               Bool bad = False;
               if ( !la1.IsProper( ) || !la2.IsProper( ) )
               {    out << "Not proper.";
                    bad = True;    }
               else if ( la1.target_id != la2.target_id )
               {    out << "Different targets.";
                    bad = True;    }
               else if ( la1.Fw1( ) != la2.Fw1( ) )
               {    out << "Inconsistent orientations.";
                    bad = True;    }
               if ( !bad )
	       {    float size_ratio;
		    int offset1 = end_offsets[j1];
		    int offset2 = end_offsets[j2];
		    start = ( la1.Fw1( ) ? la1.pos2( ) - offset1 : la2.pos2( ) + offset2 );
                    stop = ( la1.Fw1( ) ? la2.Pos2( ) - offset2 : la1.Pos2( ) + offset1 );
		    out << (la1.Fw1() ? "fw " : "rc ");
		    // Check that the implied length on the reference reasonably matches the assembly
		    size_ratio = float(assembly[i].size()) /  float(stop - start);
		    if (abs(size_ratio-1.0) < SIZE_ERROR_LIMIT) {
		      t = la1.target_id;
		      good_placements++;
		    } else {
		      bad = True;
		      out << "Size mismatch: ";
		    }
		    out << la1.target_id << "." << start << "-" << stop << " (sr=" << round(100*size_ratio)/100.0 << ")";
	       }    }
          out << ends;
          answer.push( t, make_pair(start,stop), out.str( ) );    }

     // Generate report.

     cout << Date( ) << ": generate report" << endl;
     Sort(answer);
     vec<int> gap_sizes;
     gap_sizes.reserve( answer.isize( ) );
     for ( int i = 0, maxbase = answer[0].second.second; i < answer.isize( ); i++ ) {
       if ( i > 0 && answer[i-1].first == answer[i].first 
	    && answer[i-1].first != infinity ) {
	 if (answer[i].second.second > maxbase) { // if we are extending
	   int gap_size = answer[i].second.first - maxbase;
	   report << endl << "gap = " << gap_size << endl;
	   maxbase = answer[i].second.second;
	   if ( gap_size > 0 ) gap_sizes.push_back( gap_size );
	 } else report << "\n    ";
       } else report << endl;
       report << answer[i].third;
#ifdef notdef
       if (answer[i].first == infinity) {
	 int scaffold_size = assembly[i].isize();
	 basevector e;
	 vecbasevector contigs;
	 vec<int> offsets(2*assembly.size());

	 for (int j = 0, k = 0; k <= scaffold_size; k++) {
	   if (k == scaffold_size || assembly_amb[i][k]) {
	     if (k-j >= END_SIZE_MIN) {
	       e.SetToSubOf(assembly[i], j, k-j);
	       contigs.push_back_reserve(e);
	       offsets.push_back(j);
	       report << "contig [" << j << "," << k << ") l=" << k-j << endl;
	     }
	     j = k+1;
	   }
	 }
       }
#endif
     }
     
     
     // Write report summary.
     double pct_good = 100.0 * ( (double)good_placements / assembly.size( ) );
     Sort( gap_sizes );
     
     
     report << "\n\nLAYOUT SUMMARY:" << endl;
     report << "Well-placed scaffolds: " << pct_good << "%\t("
	    << good_placements << "/" << assembly.size( ) << ")\n";
     report << "N50 gap size: " << 
       (gap_sizes.empty() ? "No Gaps" : ToString(N50( gap_sizes ))) << endl;

     cout << Date( ) << ": Done!" << endl;
     return 0;
}
