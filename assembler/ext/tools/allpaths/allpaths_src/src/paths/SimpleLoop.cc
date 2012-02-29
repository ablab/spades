/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/SimpleLoop.h"
#include "random/NormalRandom.h"
#include "random/Random.h"

void GetSimpleLoops( const HyperKmerPath& h, vec<simple_loop>& loops )
{    loops.clear( );
     for ( int v = 0; v < h.N( ); v++ )
     {    if ( h.From(v).size( ) != 2 || h.To(v).size( ) != 2 ) continue;
          int fself = -1, fother = -1, tself = -1, tother = -1;
          for ( int j = 0; j < h.From(v).isize( ); j++ )
          {    int w = h.From(v)[j];
               if ( w == v ) fself = j;
               else fother = j;    }
          for ( int j = 0; j < h.To(v).isize( ); j++ )
          {    int u = h.To(v)[j];
               if ( u == v ) tself = j;
               else tother = j;    }
          if ( fself < 0 || fother < 0 || tself < 0 || tother < 0 ) continue;
          int u = h.To(v)[tother], w = h.From(v)[fother];
          if ( u == w ) continue;
          int uv = h.EdgeObjectIndexByIndexTo( v, tother );
          int vv = h.EdgeObjectIndexByIndexFrom( v, fself );
          int vw = h.EdgeObjectIndexByIndexFrom( v, fother );
          loops.push_back( simple_loop( u, v, w, uv, vv, vw ) );    }    }

void DisambiguateSimpleLoops( HyperKmerPath& h, const vec<simple_loop>& loops,
     const vec< vec< pair<int,int> > >& sep_dev )
{    for ( int i = 0; i < sep_dev.isize( ); i++ )
     {    static vec< pair<int,int> > sepdev;
          sepdev = sep_dev[i];
          Sort(sepdev);
          if ( sepdev.size( ) >= 6 )
          {    
               // The elements of sepdev are the predicted separations in kmers
               // between the two edges that bracket the loop.

               double sep = 0.0, dev = 0.0;
               int start = sepdev.isize( ) / 3;
               int stop = sepdev.isize( ) - start;
               for ( int j = start; j < stop; j++ )
               {    sep += sepdev[j].first;
                    dev += sepdev[j].second;    }
               sep /= double( stop - start );
               dev /= double( stop - start );
               dev /= sqrt( double( sepdev.size( ) ) );

               // Now sep is the mean of the predicted separations and dev is the
               // mean, divided by the square root of the number of separations.

               int v = loops[i].v, vv = loops[i].vv;
               int vvl = h.EdgeLength(vv);

               // Let m = our best guess as to the number of times the loop is 
               // traversed.  Let md = the error in deviations corresponding to this
               // value of m.

               int m = int( round( sep / double(vvl) ) );
               double md = Abs( sep - double(vvl * m) ) / dev;

               // Let mdalt be the error in deviations corresponding to the second
               // best guess.

               double mdalt = Abs( sep - double(vvl * (m+1)) ) / dev;
               if ( m - 1 >= 0 )
               {    mdalt = Min( mdalt,
                         Abs( sep - double(vvl * (m-1)) ) / dev );    }
               if ( md > 3.0 || mdalt - md < 3.0 ) continue;
               ForceAssert( h.From(v).size( ) == 2 );
               ForceAssert( h.To(v).size( ) == 2 );
               int fself = -1, fother = -1, tself = -1, tother = -1;
               for ( int j = 0; j < h.From(v).isize( ); j++ )   
               {    int w = h.From(v)[j];
                    if ( w == v ) fself = j;
                    else fother = j;    }
               for ( int j = 0; j < h.To(v).isize( ); j++ )
               {    int u = h.To(v)[j];
                    if ( u == v ) tself = j;
                    else tother = j;    }
               ForceAssert( fself >= 0 && fother >= 0 );
               ForceAssert( tself >= 0 && tother >= 0 );
               int u = h.To(v)[tother], w = h.From(v)[fother];
               ForceAssert( u != w );
               int uv = h.EdgeObjectIndexByIndexTo( v, tother );
               int vw = h.EdgeObjectIndexByIndexFrom( v, fother );
               KmerPath p = h.EdgeObject(uv);
               for ( int j = 0; j < m; j++ )
                    p.Append( h.EdgeObject(vv) );
               p.Append( h.EdgeObject(vw) );
               h.DeleteEdgeFrom( v, fself );
               h.JoinEdges( v, p );    }    }
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );    }

double VisibleInsertLength( const int n, const HyperKmerPath& h, 
     const HyperBasevector& hb, const int uv, const int vv, const int vw,
     const int insert_sep, const int insert_dev, NormalRandom& insert, 
     const int readlength, const vec<Bool>& uvcov, 
     const vec<Bool>& vwcov )

{    // cout << "\nComputing insert length correction for multiplicity " 
     //      << n << "\n";
     int K = h.K( );
     int N = 10000;
     int middle = n * h.EdgeLength(vv);
     int left = Min( hb.EdgeLength(uv), insert_sep + 3 * insert_dev - middle );
     int left_shift = hb.EdgeLength(uv) - left;
     // PRINT2( left, left_shift );
     vec<int> insert_sizes;
     for ( int k = 0; k < N; k++ )
     {    int startuv = ( randomx( ) % left ) + left_shift;
          if ( startuv >= uvcov.isize( ) ) continue;
          if ( !uvcov[startuv] ) continue;
          int this_insert = int( round( insert.value( ) ) );
          int startvw = this_insert - ( hb.EdgeLength(uv) - startuv )
               - middle - readlength;
          if ( startvw < 0 || startvw >= vwcov.isize( ) ) continue;
          if ( !vwcov[startvw] ) continue;
          insert_sizes.push_back(this_insert);    }
     // cout << insert_sizes.size( ) << "/" << N << " placed";
     double mean = -1.0;
     if ( insert_sizes.nonempty( ) )
     {    mean = Mean(insert_sizes);
          // cout << ", mean length = " << mean << "\n";    
               }
     return mean;    }

// DisambiguateSimpleLoops2.  This is a fancier version of DisambiguateSimpleLoops.
// For each simple loop u ---> v <---> v ---> w, find the "uniquely coverable" parts
// of u --> v and v --> w.  For each library and each plausible copy number for
// v <---> v, do a simulation to predict the "insert length error term" arising
// from the fact that pairs of uniquely placeable reads have to land on the uniquely
// coverable parts.  Based on this and the placement of the actual pairs, predict
// the copy number of v <---> v if possible, and edit the assembly graph
// accordingly.

void DisambiguateSimpleLoops2( HyperKmerPath& h, HyperBasevector& hb, 
     const vecbasevector& reads, const PairsManager& pairs, 
     const vec<alignment_plus>& Aligns, 
     const Bool verbose )
{   
    // Set up indices.
    vec< vec<int> > Aligns_index( reads.size( ) );
    for ( int i = 0; i < Aligns.isize( ); i++ )
      Aligns_index[ Aligns[i].Id1( ) ].push_back(i);

     srandomx(1694384237);


     // Find simple loops.

     vec<simple_loop> loops;
     GetSimpleLoops( h, loops );

     // Define sep_lib.

     vec<longlong> pairs_index = pairs.getPairsIndex();
     vec< pair<int,int> > rp;
     for ( size_t i = 0; i < reads.size( ); i++ )
          rp.push_back( make_pair( pairs_index[i], i ) );
     Sort(rp);
     vec< vec< pair<int,int> > > lib_sep( loops.size( ) );

     for ( int i = 0; i < rp.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < rp.isize( ); j++ )
               if ( rp[j].first != rp[i].first ) break;
          if ( j - i == 2 )
          {    int id1 = rp[i].second, id2 = rp[i+1].second;
               const int pair_ID = rp[i].first;
               if ( Aligns_index[id1].solo( ) && Aligns_index[id2].solo( ) )
               {    alignment_plus ap1 = Aligns[ Aligns_index[id1][0] ];
                    alignment_plus ap2 = Aligns[ Aligns_index[id2][0] ];
                    if ( ap1.Rc2( ) != ap2.Rc2( ) )
                    {    if ( ap1.Rc2( ) )
                         {    swap( id1, id2 );
                              swap( ap1, ap2 );    }
                         packalign a2(ap2.a);
                         a2.ReverseThis( reads[id2].size( ), 
                              hb.EdgeLength( ap2.Id2( ) ) );
                         int uv = ap1.Id2( ), vw = ap2.Id2( );
                         int li = -1;
                         for ( int l = 0; l < loops.isize( ); l++ )
                         {    if ( loops[l].uv == uv && loops[l].vw == vw )
                              {    li = l;    }    }
                         if ( li < 0 )
                         {    i = j - 1;
                              continue;    }
                         int sep = pairs.sep(pair_ID) - a2.pos2( ) + h.K( ) - 1
                              - ( hb.EdgeLength( ap1.Id2( ) ) - ap1.a.Pos2( ) );
                         int lib = pairs.libraryID( sep, pairs.sd(pair_ID) );
                         lib_sep[li].push_back( 
                              make_pair( lib, sep ) );    }    }    }
          i = j - 1;    }

     // Disambiguate simple loops.

     int K = h.K( );
     int readlength = reads[0].size( ); // Dangerous!

     // Compute coverage of left and right edges of each loop by uniquely placed reads.

     vec< vec<Bool> > uvcovs( loops.size() );
     vec< vec<Bool> > vwcovs( loops.size() );
     vec< vec<int> > edgeToLoopMap( h.EdgeObjectCount() );

     for ( int li = 0; li < loops.isize( ); li++ ) {
       uvcovs[li].resize( h.EdgeLength(loops[li].uv), False );
       vwcovs[li].resize( h.EdgeLength(loops[li].vw), False );
       edgeToLoopMap[ loops[li].uv ].push_back( li );
       edgeToLoopMap[ loops[li].vw ].push_back( li );
     }

     for ( int j = 0; j < Aligns.isize(); ++j ) {
       const alignment_plus& ap = Aligns[j];
       if ( !Aligns_index[ ap.Id1() ].solo() ) continue;
       int edge = ap.Id2();
       if ( edgeToLoopMap[ edge ].empty() ) continue;
       int pos2 = ap.pos2();
       if ( ap.Rc2() ) {
         packalign a(ap.a);
         a.ReverseThis( reads[ ap.Id1( ) ].size( ), 
                        hb.EdgeLength( edge ) );
         pos2 = a.pos2();
       }
       for ( int lii = 0; lii < edgeToLoopMap[edge].isize(); ++lii ) {
         int li = edgeToLoopMap[edge][lii];
         if ( edge == loops[li].uv ) uvcovs[li][pos2] = True;
         if ( edge == loops[li].vw ) vwcovs[li][pos2] = True;
       }
     }

     // Process each loop.

     for ( int li = 0; li < loops.isize( ); li++ )
     {    
          int uv = loops[li].uv, vv = loops[li].vv, vw = loops[li].vw;
          if (verbose)
          {    cout << "\nProcessing loop " << li << "\n";
               int u = loops[li].u, v = loops[li].v, w = loops[li].w;
               PRINT3( u, v, w );
               PRINT3( BaseAlpha(uv), BaseAlpha(vv), BaseAlpha(vw) );
               PRINT3( h.EdgeLength(uv), h.EdgeLength(vv), h.EdgeLength(vw) );    }
          vec<Bool>& uvcov = uvcovs[li];
          vec<Bool>& vwcov = vwcovs[li];

          // Go through the placed pairs.

          Sort( lib_sep[li] );
          for ( int z1 = 0; z1 < lib_sep[li].isize( ); z1++ )
          {    int lib = lib_sep[li][z1].first;
               int z2;
               for ( z2 = z1 + 1; z2 < lib_sep[li].isize( ); z2++ )
                    if ( lib_sep[li][z2].first != lib_sep[li][z1].first ) break;
               const int min_pairs = 6;
               int npairs = z2 - z1;

               int insert_sep = pairs.getLibrarySep(lib); 
               int insert_dev = pairs.getLibrarySD (lib);
               NormalRandom insert( insert_sep, insert_dev );

               if (verbose)
               {    cout << "\nLibrary " << lib << " ("
                         << insert_sep << " +/- " << insert_dev << ")"
                         << ", have " << npairs << " placed pairs.\n";    }

               if ( npairs < min_pairs )
               {    z1 = z2 - 1;
                    if (verbose) cout << "Not enough pairs.\n";
                    continue;    }

               // Now lib_sep[li][z].second, z in [z1,z2), gives the predicted 
               // separation in kmers between the two edges that bracket the loop.
               // However, added to it should be the correction term
               // VisibleInsertLength( n, ... ) - insert_sep, where n is the 
               // true number of copies of the loop.  But of course we don't know n 
               // at this point.

               double sep = 0.0;
               double dev = insert_dev;
               int start = npairs / 3;
               int stop = npairs - start;
               for ( int j = start; j < stop; j++ )
                    sep += lib_sep[li][ z1 + j ].second;
               sep /= double( stop - start );
               dev /= sqrt( double(npairs) );
               if (verbose) PRINT2( sep, dev );

               // Now sep is the mean of the predicted separations and dev is the
               // mean, divided by the square root of the number of separations.
               // But we still have to deal with the correction term.

               int v = loops[li].v;
               int vvl = h.EdgeLength(vv);
               if (verbose) PRINT(vvl);

               // Let m0 = our initial best guess as to the number of times the 
               // loop is traversed.  Let md0 = the error in deviations 
               // corresponding to this guess.

               int m0 = Max( 0, int( round( sep / double(vvl) ) ) );
               double md0 = Abs( sep - double(vvl * m0) ) / dev;
               if (verbose) PRINT2( m0, md0 );

               // Decide on range of values of m that we're going to explore.

               const int m_radius = 8;
               int mlow = Max( 0, m0 - m_radius ), mhigh = m0 + m_radius;

               // Now taking into account correction term, compute md = error in
               // deviations corresponding to the values of m that we're exploring.

               vec<Bool> have_md( mhigh - mlow + 1, True );
               vec<double> md( mhigh - mlow + 1 );
               for ( int m = mlow; m <= mhigh; m++ )
               {    double see = VisibleInsertLength( m, h, hb, uv, vv, vw, 
                         insert_sep, insert_dev, insert, readlength, 
                         uvcov, vwcov );
                    if ( see < 0 ) have_md[m-mlow] = False;
                    else 
                    {    double correction = see - double(insert_sep);
                         md[m-mlow] = Abs( correction + sep - double(vvl*m) ) / dev;
                         if (verbose)
                         {    cout << "m = " << m << ", md = " << md[m-mlow] 
                                   << "\n";    }    }    }

               // See if there is a clear winner.

               vec< pair<double,int> > md_m;
               for ( int m = mlow; m <= mhigh; m++ )
               {    if ( have_md[m-mlow] )
                         md_m.push_back( make_pair( md[m-mlow], m ) );    }
               Sort(md_m);
               if ( md_m.nonempty( ) && md_m[0].first <= 3.0
                    && ( md_m.solo( ) || md_m[1].first - md_m[0].first >= 3.0 ) )
               {    
                    // Edit assembly.

                    int m = md_m[0].second;
                    double md = md_m[0].first;
                    if (verbose)
                         cout << "\nusing m = " << m << ", md = " << md << "\n";
                    ForceAssert( h.From(v).size( ) == 2 );
                    ForceAssert( h.To(v).size( ) == 2 );
                    int fself = -1, fother = -1, tself = -1, tother = -1;
                    for ( int j = 0; j < h.From(v).isize( ); j++ )   
                    {    int w = h.From(v)[j];
                         if ( w == v ) fself = j;
                         else fother = j;    }
                    for ( int j = 0; j < h.To(v).isize( ); j++ )
                    {    int u = h.To(v)[j];
                         if ( u == v ) tself = j;
                         else tother = j;    }
                    ForceAssert( fself >= 0 && fother >= 0 );
                    ForceAssert( tself >= 0 && tother >= 0 );
                    int u = h.To(v)[tother], w = h.From(v)[fother];
                    ForceAssert( u != w );
                    int uv = h.EdgeObjectIndexByIndexTo( v, tother );
                    int vw = h.EdgeObjectIndexByIndexFrom( v, fother );
                    KmerPath p = h.EdgeObject(uv);
                    for ( int j = 0; j < m; j++ )
                         p.Append( h.EdgeObject(vv) );
                    p.Append( h.EdgeObject(vw) );
                    h.DeleteEdgeFrom( v, fself );
                    h.JoinEdges( v, p );
                    break;    }

               // Update placed pair counter.

               z1 = z2 - 1;    }    }

     // Clean up assembly graph.

     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );    }


