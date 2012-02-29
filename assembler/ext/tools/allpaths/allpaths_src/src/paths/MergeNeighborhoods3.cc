///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MergeNeighborhoods3.

#include "MainTools.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/HyperKmerPath.h"

int main( int argc, char *argv[] )
{
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String( PRE );
     CommandArgument_String( DATA );
     CommandArgument_String( RUN );
     CommandArgument_String_OrDefault( SUBDIR, "test" );
     CommandArgument_Int( K );
     CommandArgument_Bool_OrDefault( WRITE, True );
     CommandArgument_String_OrDefault( OUT_SUFFIX, "" );
     EndCommandArguments;

     // Define heuristic constants.

     const int min_unique_in_component = 1000;
  
     // Define filenames.

     cout << Date( ) << ": beginning MergeNeighborhoods3" << endl;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String kK = ".k" + ToString( K );
     String data_dir = PRE + "/" + DATA;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     // Read HyperKmerPath.

     HyperKmerPath h;
     BinaryRead( sub_dir + "/hyper.prelim" + OUT_SUFFIX, h );
  
     // Delete components having less than 1000 unique kmers.  Go from smallest 
     // to largest.  The second sort uses a lot of memory.  The amount of memory
     // is increased a lot by making it a parallel sort.

     cout << Date( ) << ": set up for cleanup" << endl;
     vec< vec<int> > comps;
     h.Components(comps);
     vec<int64_t> kmers;
     vec<int> o, keep, kmers_orig, nkmers( comps.size( ), 0 );
     for ( int i = 0; i < comps.isize( ); i++ )
     {    const vec<int>& o = comps[i];
          for ( int j = 0; j < o.isize( ); j++ )
          {    int v = o[j];
               for ( int t = 0; t < h.From(v).isize( ); t++ )
               {    nkmers[i] += h.EdgeObjectByIndexFrom( v, t )
                         .KmerCount( );    }    }    }
     cout << Date( ) << ": first sort" << endl;
     SortSync( nkmers, comps );
     for ( int i = 0; i < comps.isize( ); i++ )
     {    const vec<int>& o = comps[i];
          int nkmers = 0;
          for ( int j = 0; j < o.isize( ); j++ )
          {    int v = o[j];
               for ( int t = 0; t < h.From(v).isize( ); t++ )
               {    const KmerPath& p = h.EdgeObjectByIndexFrom( v, t );
                    int n = p.KmerCount( );
                    for ( int l = 0; l < n; l++ )
                    {    kmer_id_t x = p.GetKmer(l);
                         kmer_id_t xrc = reverse_kmer(x);
                         kmers.push_back( Min( x, xrc ) );
                         kmers_orig.push_back(i);    }    }    }    }
     cout << Date( ) << ": second sort" << endl;
     SortSync( kmers, kmers_orig );
     cout << Date( ) << ": find stuff to remove" << endl;
     vec<Bool> to_remove( comps.size( ), False );
     for ( int i = 0; i < comps.isize( ); i++ )
     {    const vec<int>& o = comps[i];
          int n_unique = 0;
          for ( int j = 0; j < o.isize( ); j++ )
          {    int v = o[j];
               for ( int t = 0; t < h.From(v).isize( ); t++ )
               {    const KmerPath& p = h.EdgeObjectByIndexFrom( v, t );
                    int n = p.KmerCount( );
                    for ( int l = 0; l < n; l++ )
                    {    kmer_id_t x = p.GetKmer(l);
                         kmer_id_t xrc = reverse_kmer(x);
                         x = Min( x, xrc );
                         size_t start = lower_bound( 
                              kmers.begin( ), kmers.end( ), x ) - kmers.begin( );
                         size_t stop = upper_bound( 
                              kmers.begin( ), kmers.end( ), x ) - kmers.begin( );
                         int mult = 0;
                         for ( size_t u = start; u < stop; u++ )
                              if ( !to_remove[ kmers_orig[u] ] ) ++mult;
                         if ( mult == 1 ) ++n_unique;    }    }    }
          if ( n_unique < min_unique_in_component ) to_remove[i] = True;
          else keep.append(o);    }
     h = HyperKmerPath( h, keep );
     h.RemoveDeadEdgeObjects( );

     // Report statistics on components.

     vec< vec<int> > compe;
     h.ComponentsE(compe);
     vec<int64_t> compsizes;
     for ( int i = 0; i < compe.isize( ); i++ )
     {    int64_t total = 0;
          for ( int j = 0; j < compe[i].isize( ); j++ )
               total += h.EdgeObject( compe[i][j] ).KmerCount( );
          compsizes.push_back(total);    }
     cout << compe.size( ) << " components of total size "
          << ToStringAddCommas( Sum(compsizes) );
     if ( compsizes.nonempty( ) ) cout << " and N50 " << N50(compsizes);
     cout << endl;

     // Write.

     if (WRITE)
     {    cout << Date( ) << ": write assembly to file <sub_dir>/hyper" + 
               OUT_SUFFIX << endl;
          BinaryOverwrite( sub_dir + "/hyper" + OUT_SUFFIX, h );
          Ofstream( dot, sub_dir + "/hyper" + OUT_SUFFIX + ".dot" );
          h.PrintSummaryDOT0w( dot, True, False, True, 0, True );    }

     cout << Date( ) << ": done" << endl;
     _exit(0);    }
