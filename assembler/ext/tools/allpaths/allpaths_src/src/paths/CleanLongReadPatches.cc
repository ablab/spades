///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// CleanLongReadPatches.  Read in long read patches and try to improve them using
// the unipaths.  In more detail,

/*
The idea is to cover the patch with unipaths, then infer the patch sequence from the 
unipaths.  However, we don't want to be restricted to the case where the unipaths 
provide 'continuous' coverage of the region (i.e. cover with successive overlaps of 
K-1), so it's a bit tricky.

Here's the picture:

      left contig bases             patch              right contig bases
----------------------------xxxxxxxxxxxxxxxxxxxxx--------------------------------

Here's the algorithm:

1. Find the last K-mer amongst the left contig bases that is in a unipath.  This 
unipath is the "left" unipath.  Find the first K-mer amongst the right contig bases 
that is in a unipath.  This unipath is the "right" unipath.  Then truncate the 
picture, trimming bases not between the far ends of the two bounding K-mers.

                   ---------xxxxxxxxxxxxxxxxxxxxx-----------

This is the "target", denoted t.

2. By Smith-Waterman, we shall mean the Smith-Waterman alignment algorithm with 
mismatch penalty = 2, gap base penalty = 3.

3. Let L = 12.  Seeding on L-mers, align the unipaths to t.  This is complicated 
because we need to allow the same unipath to go in more than one place.  To take the 
most extreme case, the unipath might be a homopolymer, with a long series of
staggered placements.  At the same time it would be better not to have multiple 
instantiations of what is effectively the same placement.  We proceed as follows:
(a) Make a list of all the L-mer matches between unipath u and the target t.
(b) Choose a match x from the list.  Divide u into three nonoverlapping pieces 
    l, m, r: the middle piece m is L bases long and is defined by x.  The left 
    [resp. right] piece l [resp. r] consists of the bases in u to the left 
    [resp. right] of m.
(c) Now separately Smith-Waterman l and r to the target.  These alignments are 
    completely free, except that in each case the end adjacent to m is fixed 
    (meaning that we penalize for indels there).
(d) Discard alignments having mismatch rate > 10%.
(e) The alignment is now accepted.  Remove all matches in the list that are subsumed
    by the alignment and whose offset is within 5 of x.
(f) Repeat starting at (b), so long as matches remain.

4. Screen the alignments, as follows.  Sort the alignments by mismatch rate.  Now
walk through the sorted alignments and progressively cover the target with them.
If we encounter an alignment that adds nothing to coverage, delete it.

5. For each alignment, record its matching L-mers.  We now think of each alignment 
as being represented by this sequence of L-mer matches.  Discard duplicate 
alignments.

6. The left and right unipaths are aligned to the target at the loci defined by
the K-mers they share with the contig.

7. If alignments of the same unipath share an L-mer match, merge them together into
a single alignment.

8. Form a digraph G.  The vertices of this digraph are the aligned unipaths.  An 
edge from aligned unipath u1 to aligned unipath u2 is defined by an integer (the 
"offset").  An offset is only allowed if it arises from an L-mer that is shared
between u1, u2, and the target, and if it defines a perfect overlap between u1 
and u2, and if in the alignment, u2 extends beyond u1.

9. Find all paths in the graph from the left unipath to the right unipath.  Fail if 
too many paths (> 10^5).  If there is no path, the method fails.

10. Each path defines a possible sequence that replaces the target.  Smith-Waterman 
the path sequences with the target, holding the ends fixed.  Choose the path 
sequence that matches best (and if there is a tie, pick a winner at random).  This 
is the answer. 
*/

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "Fastavector.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "ParallelVecUtilities.h"
#include "Set.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "kmers/KmerRecord.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/AssemblyEdit.h"
#include "paths/LongReadTools.h"

template<int K> void Main( vec< vec<assembly_edit> >& edits, const vecbasevector& U,
     const int L, const vec< vec< pair<int,int> > >& Ulocs, 
     const vec<superb>& scaffolds, const vec<fastavector>& tigsa, 
     const int verbosity, const Bool DIRECT, const String& data_dir, 
     const String& run_dir )
{
     // Hash U.

     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < U.size( ); i++ )
     {    const basevector& u = U[i];
          starts.push_back( starts.back( ) + u.isize( ) - K + 1 );    }
     vec< triple<kmer<K>,int,int> > kmers_plus( starts.back( ) );
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( size_t i = 0; i < U.size( ); i++ )
     {    const basevector& u = U[i];
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               kmers_plus[r].first.SetToSubOf( u, j );
               kmers_plus[r].second = i, kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);

     // Go through the patches.

     vec< pair<int,int> > gaps;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
     {    const superb& S = scaffolds[s];
          for ( int i = 0; i < S.Ngaps( ); i++ )
               if ( edits[ S.Tig(i) ].nonempty( ) ) gaps.push( s, i );    }
     vec<String> report( gaps.size( ) );
     cout << Date( ) << ": start editing " << gaps.size( ) << " patches" 
          << " (100 dots to follow)" << endl;
     uint dots_printed = 0, gaps_processed = 0;
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int gi = 0; gi < gaps.isize( ); gi++ )
     {    
          // Log progress.
          
          uint n = gaps.size( );
          uint& i = gaps_processed;
          #pragma omp critical
          {    i++;
               uint newdots = ((100 * i) / n) - (100 * (i-1)) / n;
               if (newdots > 0) 
               {    for (uint j = 0; j < newdots; j++) 
                    {    cout << ".";
                         dots_printed++;
                         if (dots_printed % 50 == 0) cout << "\n";
                         else if (dots_printed % 10 == 0) cout << " ";    }
                    flush(cout);    }    }

          // Examine patch.

          int s = gaps[gi].first, it = gaps[gi].second;
          const superb& S = scaffolds[s];
          int m1 = S.Tig(it), m2 = S.Tig(it+1);
          ostringstream outx;
          ostream& rout = ( DIRECT ? cout : outx );
          if ( verbosity >= 1 )
          {    rout << "\n" << Date( ) 
                    << ": examining patch from " << m1 << " to " << m2 << endl;    }
          ForceAssert( edits[m1].solo( ) );
          CleanPatch( L, U, Ulocs, kmers_plus, tigsa[m1], tigsa[m2], edits[m1][0],
               outx, DIRECT, verbosity, data_dir, run_dir );
          if ( !DIRECT ) report[gi] = outx.str( );    }

     // Print reports.

     if ( verbosity >= 1 && !DIRECT )
     {    for ( int gi = 0; gi < gaps.isize( ); gi++ )
               if ( report[gi].size( ) > 0 ) cout << report[gi];    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds0.clean");
     CommandArgument_Int(K);
     CommandArgument_String_OrDefault(IN_HEAD, "all_reads");
     CommandArgument_String_OrDefault(EDITS_IN, SCAFFOLDS_IN + ".longread.edits");
     CommandArgument_String_OrDefault(EDITS_OUT,
          SCAFFOLDS_IN + ".longread.fixed.edits");
     CommandArgument_Int_OrDefault(VERBOSITY, 0);
     CommandArgument_Int_OrDefault_Doc(TIG1, -1,
          "look only at gaps whose left contig is this one");
     EndCommandArguments;

     // Start.

     double clock = WallClockTime( );

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Define directories, etc.
     
     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String tigsa_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     Mkdir777( sub_dir + "/tmp" );

     // Load assembly.

     vec<fastavector> tigsa;
     LoadFromFastaFile( tigsa_file, tigsa );
     int ntigs = tigsa.size( );
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds ); 

     // Read the edits.

     vec< vec<assembly_edit> > edits(ntigs);
     String filename =  sub_dir + "/" + EDITS_IN;
     vec<assembly_edit> editsx;
     BinaryReader::readFile( filename.c_str( ), &editsx );
     Bool DIRECT = ( TIG1 >= 0 );
     for ( size_t j = 0; j < editsx.size( ); j++ )
     {    int m1 = editsx[j].Tig1( );
          ForceAssertLt( m1, ntigs );
          if ( TIG1 < 0 || m1 == TIG1 ) edits[m1].push_back( editsx[j] );    }

     // Load the unipaths and hash them.
     
     String unibasefile = run_dir + "/" + IN_HEAD + ".unibases.k" + ToString(K);
     vecbasevector U(unibasefile);
     const int L = 12;
     vec< vec< pair<int,int> > > Ulocs( IPow( 4, L ) );
     for ( size_t i = 0; i < U.size( ); i++ ) 
     {    for ( int j = 0; j <= U[i].isize( ) - L; j++ )
          {    int n = KmerId( U[i], L, j );
               Ulocs[n].push( i, j );    }    }

     // Clean the patches.

     if ( K == 96 )
     {    Main<96>( edits, U, L, Ulocs, scaffolds, tigsa, 
               VERBOSITY, DIRECT, data_dir, run_dir );    }
     else
     {    cout << "K value unsupported" << endl;
          exit(1);    }

     // Insert the new edits and write.

     if ( TIG1 < 0 )
     {    cout << Date( ) << ": writing new edits" << endl;
          for ( size_t j = 0; j < editsx.size( ); j++ )
               editsx[j] = edits[ editsx[j].Tig1( ) ][0];
          String outfile =  sub_dir + "/" + EDITS_OUT;
          BinaryWriter::writeFile( outfile.c_str( ), editsx );    }
     cout << "time used = " << TimeSince(clock) << endl;    
     _exit(0);    }
