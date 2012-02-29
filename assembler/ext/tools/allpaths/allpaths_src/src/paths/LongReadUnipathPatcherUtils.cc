//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/LongReadTools.h"
#include "paths/LongReadUnipathPatcherUtils.h"

void BuildGLocs( const int LG,
		 const vecbasevector &genome,
		 vec< vec< pair<int,int> > > &Glocs )
{    Glocs.clear( );
     Glocs.resize( IPow( 4, LG ) );
     for ( size_t i = 0; i < genome.size( ); i++ )
     {    for ( int j = 0; j <= genome[i].isize( ) - LG; j++ )
	  {    int n = KmerId( genome[i], LG, j );
	       Glocs[n].push( i, j );    }    }    }

void AlignToGenome( const int LG,
		    const bool COMPUTE_TRUE_READ_LENGTH,
		    const basevector &r,
		    const vecbasevector &genome,
		    const vec< vec< pair<int,int> > > &Glocs,
		    vec<align_data> &adata,
		    ostream &out )
{
     // Define heuristic constants for validation.
  
     const int G_initial = 100;
     const int G_initial_bandwidth_div = 10;
     const double G_initial_error_rate_max = 0.25;
     const double G_bandwidth_frac = 0.2;

     for ( int pass = 1; pass <= 2; pass++ )
     {    basevector b = r;
          if ( pass == 2 ) b.ReverseComplement( );
          vec< triple<int,int,int> > places;
          for ( int s = 0; s <= b.isize( ) - LG; s++ )
          {    int n = KmerId( b, LG, s );
               for ( int z = 0; z < Glocs[n].isize( ); z++ )
               {    int gid = Glocs[n][z].first;
                    places.push( gid, s, Glocs[n][z].second );    }    }
          align a;
          int errors;
          vec<Bool> checked( places.size( ), False );
          for ( int q = 0; q < places.isize( ); q++ )
          {    if ( checked[q] ) continue;
               int gid = places[q].first;
               int rpos = places[q].second, gpos = places[q].third;
               const basevector& g = genome[gid];
	       // Test with an initial, cheaper alignment.
	       int flank = ( G_initial - LG ) / 2;
               int b_start = Max( 0, rpos - flank );
               int b_stop = Min( rpos + LG + flank, b.isize( ) );
               basevector b0( b, b_start, b_stop - b_start );
               int bw = G_initial / G_initial_bandwidth_div;
               int offset = ( rpos - b_start ) - gpos;
               SmithWatBandedA( b0, g, offset, bw, a, errors, 0, 1, 1 );
               if ( a.pos1( ) > 0 || a.Pos1( ) < b0.isize( ) ) continue;
               if ( errors > double( b0.size( ) )
                    * G_initial_error_rate_max )
               {    continue;    }
                // Now do the full alignment.
                     offset = rpos - gpos;
               int bandwidth = int(ceil( 
                    G_bandwidth_frac * double( r.size( ) ) ));
               int sub = 1;
               int ins = 1;
               int del = 1;
               int errors2;
               SmithWatBandedA2<unsigned short>( b, g, offset, 
                    bandwidth, a, errors2, (ostream*) 0, sub, ins, del );
               if ( a.pos1( ) > 0 || a.Pos1( ) < r.isize( ) ) continue;
               adata.push( a, gid, a.pos2( ), a.Pos2( ), errors2, 
                    pass == 1 );    
                // Mark places that are effectively checked by this align.
               // It might save time to switch to a binary search of 
               // "places".
                vec<ho_interval> perfs1, perfs2;
               a.PerfectIntervals1( b, g, perfs1 );
               a.PerfectIntervals2( b, g, perfs2 );
               for ( int w = 0; w < places.isize( ); w++ )
               {    if ( checked[w] ) continue;
                    if ( places[w].first != gid ) continue;
                    int rpos = places[w].second, gpos = places[w].third;
                    for ( int y = 0; y < perfs1.isize( ); y++ )
                    {    if ( !Subset( ho_interval( rpos, rpos + LG ),
                              perfs1[y] ) )
                         {    continue;    }
                         if ( !Subset( ho_interval( gpos, gpos + LG ),
                              perfs2[y] ) )
                         {    continue;    }
                         if ( perfs1[y].Start( ) - perfs2[y].Start( )
                              != rpos - gpos )
                         {    continue;    }
                         checked[w] = True;
                         break;    }    }    }    }
      UniqueSort(adata);
     int min_errors = ( adata.empty( ) ? -1 : adata[0].errors );
     const double max_excess_errors = 0.15;
     vec<Bool> to_delete_a( adata.size( ), False );
     for ( int i = 0; i < adata.isize( ); i++ )
     {    if ( !( adata[i].errors
                    <= double(min_errors) * ( 1.0 + max_excess_errors ) ) )
          {    to_delete_a[i] = True;    }    }
     EraseIf( adata, to_delete_a );
     out << "\n" << adata.size( ) << " genomic placements:\n\n";
     for ( int i = 0; i < adata.isize( ); i++ )
     {    int gid = adata[i].gid;
          out << "[" << i+1 << "] " << gid << "." << adata[i].pos2 << "-" 
               << adata[i].Pos2 << " " << ( adata[i].fw ? "fw" : "rc" )
               << ", " << adata[i].errors << " errors";
          basevector b = r;
          if ( !adata[i].fw ) b.ReverseComplement( );
	  // Compute "true" read length.
	  if (COMPUTE_TRUE_READ_LENGTH)
          {    vec<int> errs_at_pos( b.isize( ), 0 );
               const align& a = adata[i].a;
               int p1 = a.pos1( ), p2 = a.pos2( );
               for ( int j = 0; j < a.Nblocks( ); j++ ) 
               {    if ( a.Gaps(j) > 0 )  
                    {    errs_at_pos[p1] += a.Gaps(j);
                         p2 += a.Gaps(j);    }
                    if ( a.Gaps(j) < 0 ) 
                    {    errs_at_pos[p1] -= a.Gaps(j);
                         p1 -= a.Gaps(j);    }
                    for ( int x = 0; x < a.Lengths(j); x++ ) 
                    {    if ( b[p1] != genome[gid][p2] ) ++errs_at_pos[p1];
                         ++p1; ++p2;    }     }
               int max20 = 0;
               for ( int i1 = 0; i1 < b.isize( ); i1++ )
               {    int e = 0;
                    for ( int i2 = i1; i2 < b.isize( ); i2++ )
                    {    e += errs_at_pos[i2];
                         int n = i2 - i1 + 1;
                         if ( n >= max20 )
                         {    double err_rate = double(e)/double(n);
                              if ( err_rate <= 0.2 )
                                   max20 = Max( max20, n );    }    }    }
               out << ", true read length = " << max20 << "\n";    }
          out << "\n";
	  // Print alignment.
	  PrintVisualAlignment( False, out, b, genome[gid], adata[i].a );

	  // [OBSOLETE CODE]: these lines were in the original
	  //  LongReadUnipathPatcher: both ptl and ptl_loc were
	  //  defined but not filled in the original code (hence they
	  //  were never really used).
	  
          // int low = lower_bound( ptl_loc.begin( ), ptl_loc.end( ),
          //      make_triple( gid, adata[i].pos2, 0 ) ) - ptl_loc.begin( );
          // while( low > 0 && ptl_loc[low-1].first == gid
          //      && adata[i].pos2 < ptl_loc[low-1].third )
          // {    low--;    }
          // int high = upper_bound( ptl_loc.begin( ), ptl_loc.end( ),
          //      make_triple( gid, adata[i].Pos2, 0 ) ) - ptl_loc.begin( );
          // for ( int l = low; l < high; l++ )
          // {    if ( l > low && ptl_loc[l-1].third -
          //           ptl_loc[l].second != K - 1 )
          //      {    out << "----------------------------------------------"
          //                << "--------------------------------------\n";   }
          //      out << ptl[l] << "\n";    }

          out << "\n";    }    }
