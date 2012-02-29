///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// LongReadJoin.  IN PROGRESS.

// TO DO:
// 
// 1. At VALIDATION=2, calls to QLT all use the same temp file.
// 2. Parallelize.
//
// Note that this code will 'elide' long repeats.  We should capture these somehow,
// or at least note their existence for paths that have them.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency QueryLookupTable

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/BigMapTools.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/UnipathScaffold.h"
#include "paths/Uniseq.h"

void ReportExtension( const int l, const vec<int>& x, const int K, 
     const vecbasevector& unibases, const int VALIDATION, const int LG, 
     const vecbasevector& genome2, const vec< vec< pair<int,int> > >& Glocs,
     const String& data_dir, vec<placementy>& PLACES, vecbasevector& all,
     vec< vec<int> >& exts_all, const Bool VERBOSE )
{    exts_all.push_back(x);
     if (VERBOSE)
     {    cout << "[" << l+1 << "]";
          for ( int m = 0; m < x.isize( ); m++ )
               cout << " " << x[m];
          cout << "\n";    }
     vec<basevector> B;
     basevector b;

     // For the special case where start(x) = stop(x), we will align
     // two versions.  This is to accommodate the case were x defines an 
     // overlapping join around a circle.

     if ( x.front( ) != x.back( ) )
     {    for ( int m = 0; m < x.isize( ); m++ )
          {    b = Cat( b, unibases[ x[m] ] );
               if ( m < x.isize( ) - 1 ) b.resize( b.isize( ) - (K-1) );    }
          B.push_back(b);    }
     else
     {    b.resize(0);
          for ( int m = 0; m < x.isize( ) - 1; m++ )
          {    b = Cat( b, unibases[ x[m] ] );
               if ( m < x.isize( ) - 1 ) b.resize( b.isize( ) - (K-1) );    }
          B.push_back(b);
          b.resize(0);
          for ( int m = 1; m < x.isize( ); m++ )
          {    b = Cat( b, unibases[ x[m] ] );
               if ( m < x.isize( ) - 1 ) b.resize( b.isize( ) - (K-1) );    }
          B.push_back(b);    }
     for ( int i = 0; i < (int) B.size( ); i++ )
     {    all.push_back_reserve( B[i] );
          basevector r = B[i];
          r.ReverseComplement( );
          all.push_back_reserve(r);    }

     // Now do the alignment.

     if ( VALIDATION >= 1 )
     {    for ( int bi = 0; bi < (int) B.size( ); bi++ )
          {    const basevector& b = B[bi];
               vec<placementy> places 
                    = FindGenomicPlacementsY( 0, b, LG, genome2, Glocs );
               for ( int m = 0; m < places.isize( ); m++ )
               {    placementy p = places[m];
                    int ng = genome2[p.g].isize( ) / 2;
                    if ( p.pos >= ng ) continue;
                    PLACES.push_back(p);
                    cout << "placed: " << p.g << "." << p.pos << "-" << p.Pos << " "
                         << ( p.fw ? "fw" : "rc" ) << "\n";    }    
               if ( places.empty( ) && VALIDATION >= 2 )
               {    {    Ofstream( out, "slobber.fasta" );
                         b.Print( out, 0 );    }
                    String qlt = AllOfOutput1( "QueryLookupTable K=12 MM=12 MC=0.15 "
                         "SEQS=slobber.fasta L=" + data_dir + "/genome.lookup "
                         "KB=1 VISUAL=True NH=True QUIET=True" );
                    qlt.ReplaceBy( "...", "" );
                    qlt.GlobalReplaceBy( "\n\n\n", "\n\n" );
                    cout << qlt;    }    }    }    }

int main(int argc, char *argv[])
{
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_Int(KBIG);
     CommandArgument_String_OrDefault(IN_HEAD, "extended");
     CommandArgument_String_OrDefault(OUT_HEAD, IN_HEAD + ".long");
     CommandArgument_Int_OrDefault_Doc(VALIDATION, 0,
          "if 1, report perfect alignments; if 2, align everything");
     CommandArgument_Bool_OrDefault(PRINT_INPUT, False);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
          "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(DOT, "");
     CommandArgument_Bool_OrDefault(CIRCO, True);
     CommandArgument_Bool_OrDefault(USE_ALL, True);
     CommandArgument_Bool_OrDefault(USE_UNI, False);
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Int_OrDefault(MIN_SAFE_CN1, 0);
     EndCommandArguments;

     // Heuristics.

     const int min_seed = 1000;
     const int max_overlap_diff = 100;
     const int max_pile = 100;
     const int max_gap = 15000;

     // Define directories, etc.

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String outhead = run_dir + "/" + OUT_HEAD;

     // Thread control.

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads(NUM_THREADS);

     // Load unibases.

     vecbasevector unibases( run_dir + "/" + IN_HEAD + ".unibases.k" + ToString(K) );
     int nuni = unibases.size( );
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc, K );

     // Load right_exts.

     vec< vec<int> > right_exts;
     BinaryReader::readFile( ( run_dir + "/extended.long.right_exts" ).c_str( ),
          &right_exts );
     cout << Date( ) << ": found " << right_exts.size( ) << " right exts" << endl;
     if ( MIN_SAFE_CN1 > 0 )
     {    vec<Bool> to_delete( right_exts.size( ), False );
          for ( int i = 0; i < right_exts.isize( ); i++ )
          {    int u = right_exts[i][0];
               if ( unibases[u].isize( ) < MIN_SAFE_CN1 ) to_delete[i] = True;    }
          EraseIf( right_exts, to_delete );
          cout << Date( ) << ": now have " << right_exts.size( ) << " right exts"
               << endl;    }

     // Load predicted gaps and use them to make a digraph.

     String KS = ToString(K);
     cout << Date( ) << ": loading predicted gaps" << endl;
     digraphE<linklet> H;
     {    vec< vec<int> > from(nuni), to(nuni);
          vec< vec<int> > from_edge_obj(nuni), to_edge_obj(nuni);
          vec<linklet> edges;
          String line;
          fast_ifstream in( run_dir + "/" + IN_HEAD + ".unibases.k" 
               + KS + ".predicted_gaps.txt" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               int u1, u2, sep, dev, nlinks;
               istringstream iline( line.c_str( ) );
               iline >> u1 >> u2 >> sep >> dev >> nlinks;
               if ( unibases[u1].isize( ) < min_seed ) continue;
               if ( unibases[u2].isize( ) < min_seed ) continue;
               if ( u2 == to_rc[u1] ) continue; 
               from[u1].push_back(u2), to[u2].push_back(u1);
               from_edge_obj[u1].push_back( edges.size( ) );
               to_edge_obj[u2].push_back( edges.size( ) );
               edges.push( sep, dev, nlinks, 0 );    }
          for ( int u = 0; u < nuni; u++ )
          {    SortSync( from[u], from_edge_obj[u] );
               SortSync( to[u], to_edge_obj[u] );    }
          H.Initialize( from, to, edges, to_edge_obj, from_edge_obj );    }

     // Get genome.

     vecbasevector genome, genome2;
     const int LG = 12;
     vec< vec< pair<int,int> > > Glocs;
     if ( VALIDATION >= 1 )
     {    genome.ReadAll( data_dir + "/genome.fastb" );
          genome2.resize( genome.size( ) );
          for ( size_t j = 0; j < genome.size( ); j++ )
               genome2[j] = Cat( genome[j], genome[j] );
          Glocs.resize( IPow( 4, LG ) );
          for ( size_t i = 0; i < genome2.size( ); i++ )
          {    for ( int j = 0; j <= genome2[i].isize( ) - LG; j++ )
               {    int n = KmerId( genome2[i], LG, j );
                    Glocs[n].push( i, j );    }    }    }

     // Form word set (all right exts and their reverse complements) and build 
     // index to it.

     vec< vec<int> > seqs = right_exts;
     for ( int i = 0; i < right_exts.isize( ); i++ )
     {    vec<int> x = right_exts[i];
          x.ReverseMe( );
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = to_rc[ x[j] ];
          seqs.push_back(x);    }
     vec< vec< pair<int,int> > > index(nuni);  // (seq id, pos)
     for ( int i = 0; i < seqs.isize( ); i++ )
     {    const vec<int>& x = seqs[i];
          for ( int j = 0; j < x.isize( ); j++ )
               index[ x[j] ].push( i, j );    }

     // Print the right exts.

     if (PRINT_INPUT)
     {    for ( int i = 0; i < right_exts.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < right_exts.isize( ); j++ )
                    if ( right_exts[j][0] != right_exts[i][0] ) break;
               for ( int k = i; k < j; k++ )
               {    const vec<int>& v = right_exts[k];
                    int len = K-1;
                    for ( int m = 0; m < v.isize( ); m++ )
                         len += unibases[ v[m] ].isize( ) - (K-1);
                    cout << "[" << v[0] << "." << k-i+1 << ",l=" << len << "]";
                    for ( int m = 0; m < v.isize( ); m++ )
                         cout << " " << v[m];
                    cout << "\n";    }    
               i = j - 1;   }    }

     // Go through the right exts and try to extend them.

     vec<placementy> PLACES;
     vecbasevector all;
     vec< vec<int> > exts_all;
     for ( int i = 0; i < right_exts.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < right_exts.isize( ); j++ )
               if ( right_exts[j][0] != right_exts[i][0] ) break;
          for ( int zz = i; zz < j; zz++ )
          {    int min_overlap = 1;
               const vec<int>& x = right_exts[zz];
               int u = x[0];

               // Set if there's nothing to do.
     
               if ( unibases[u].isize( ) < min_seed ) continue;
               Bool self = False;
               for ( int m = 1; m < x.isize( ); m++ )
               {    if ( unibases[ x[m] ].isize( ) >= min_seed ) 
                    {    self = True;
                         if (VERBOSE)
                         {    cout << "\nextensions of " << u << "." 
                                   << zz-i+1 << ":\n";    }
                         vec<int> x2(x);
                         x2.resize(m+1);
                         ReportExtension( 0, x2, K, unibases, VALIDATION, LG, 
                              genome2, Glocs, data_dir, PLACES, all, exts_all,
                              VERBOSE );
                         break;    }    }
               if (self) continue;

               // Incrementally builds extensions of x.

               map< vec<int>, int > exts, exts_done; // x' --> minimum overlap
               const int infinity = 1000000000;
               exts[x] = infinity;
               set< vec<int> > elides;
               typedef map< vec<int>, int >::iterator emapit;
               while( exts.size( ) > 0 )
               {    
                    // Test for blowup.
     
                    if ( (int) ( exts.size( ) + exts_done.size( ) ) > max_pile )
                    {    int m = 1000000000;
                         for ( emapit r = exts.begin( ); r != exts.end( ); r++ )
                              m = Min( m, r->second );
                         for ( emapit r = exts_done.begin( ); 
                              r != exts_done.end( ); r++ )
                         {    m = Min( m, r->second );    }
                         min_overlap = m + 1;
                         // cout << "raising min overlap to " << m+1 << endl; // XXX

                         vec<emapit> erase_exts, erase_exts_done;
                         for ( emapit r = exts.begin( ); r != exts.end( ); r++ )
                         {    if ( r->second < min_overlap ) 
                                   erase_exts.push_back(r);    }
                         for ( int j = 0; j < erase_exts.isize( ); j++ )
                              exts.erase( erase_exts[j] );
                         for ( emapit r = exts_done.begin( ); 
                              r != exts_done.end( ); r++ )
                         {    if ( r->second < min_overlap ) 
                                   erase_exts_done.push_back(r);    }
                         for ( int j = 0; j < erase_exts_done.isize( ); j++ )
                              exts_done.erase( erase_exts_done[j] );
                         if ( exts.size( ) == 0 ) break;    }
     
                    // Pop an extension.
     
                    vec<int> y = exts.begin( )->first;
                    int over = exts.begin( )->second;
                    exts.erase( exts.begin( ) );
     
                    /*
                    cout << "\nexamining";
                    for ( int l = 0; l < y.isize( ); l++ )
                         cout << " " << y[l];
                    cout << endl;
                    */

                    // Find the right extensions of y.

                    int ylen = K-1;
                    for ( int r = 1; r < y.isize( ); r++ )
                         ylen += unibases[ y[r] ].isize( ) - (K-1);
                    vec< triple<int,int,int> > poss;
                    int pos1 = y.isize( ) - 1;
                    {    for ( int w = 0; w < index[ y[pos1] ].isize( ); w++ )
                         {    int l = index[ y[pos1] ][w].first;
                              int pos2 = index[ y[pos1] ][w].second;
                              const vec<int>& z = seqs[l];
                              int o = pos1 - pos2;
                              if ( y.isize( ) - o >= z.isize( ) ) continue;
                              Bool mismatch = False;
                              for ( int m = Max( -o, 0 ); 
                                   m < Min( y.isize( ) - o, z.isize( ) ); m++ )
                              {    if ( y[o+m] != z[m] ) 
                                   {    mismatch = True;
                                        break;    }    }
                              if (mismatch) continue;
                              vec<int> ynew(y);
                              Bool terminal = False;
                              for ( int r = y.isize( ) - o; r < z.isize( ); r++ )
                              {    ynew.push_back( z[r] );
                                   if ( unibases[ z[r] ].isize( ) >= min_seed )
                                   {    terminal = True;
                                        break;    }    }

                              // Check for 'illegal' repeat.  First compute the
                              // length rn of the longest terminal sequence in ynew
                              // that appears earlier in ynew.

                              if ( !terminal )
                              {    int rn = 0;
                                   for ( int r = 0; r < ynew.isize( ) - 1; r++ )
                                   {    for ( int s = r; s >= 0; s-- )
                                        {    if ( ynew[ ynew.isize( ) - 1 - (r-s) ]
                                                  != ynew[s] )
                                             {    rn = Max( rn, r - s );    
                                                  break;    }    }    }

                                   // Length of terminal sequence must be >= KBIG.

                                   int rnlen = 0, nyn = ynew.size( );
                                   for ( int t = 0; t < rn; t++ )
                                   {    rnlen += unibases[ nyn - t - 1 ].isize( )
                                             - (K-1);    }
                                   if ( rnlen >= KBIG )
                                   {
                                        // There must not exist any read that aligns
                                        // to the repeat, extending it on both sides.

                                        Bool biextended = False;

                                        int p1 = nyn - 1;
                                        for ( int w = 0; w 
                                             < index[ ynew[p1] ].isize( ); w++ )
                                        {    int l = index[ ynew[p1] ][w].first;
                                             int p2 = index[ ynew[p1] ][w].second;
                                             int o = p1 - p2;
                                             if ( nyn - rn - o <= 0 ) continue;
                                             const vec<int>& z = seqs[l];
                                             if ( nyn - o >= z.isize( ) ) continue;
                                             Bool mismatch = False;
                                             for ( int m = Max( nyn - o - rn, 0 ); 
                                                  m < Min( nyn - o, z.isize( ) ); 
                                                  m++ )
                                             {    if ( ynew[o+m] != z[m] ) 
                                                  {    mismatch = True;
                                                       break;    }    }
                                             if ( !mismatch ) 
                                             {    biextended = True;
                                                  break;    }    }
                                        if ( !biextended ) 
                                        {    vec<int> rep;
                                             for ( int l = rn; l >= 1; l-- )
                                                  rep.push_back( ynew[ nyn - l ] );
                                             elides.insert(rep);
                                             continue;    }    }    }

                              // Check for too long.

                              /*
                              if ( !terminal )
                              {    int len = ylen;
                                   for ( int r = y.isize( ); r < ynew.isize( ); r++ )
                                   {    len += unibases[ ynew[r] ].isize( ) - (K-1);
                                        if ( len > max_gap ) break;    }
                                   if ( len > max_gap )
                                   {    cout << "\nhuge:";
                                        for ( int l = 0; l < ynew.isize( ); l++ )
                                             cout << " " << ynew[l];
                                        cout << "\n";
                                        return 0;    }
                                   if ( len > max_gap ) continue;    }
                              */

                              int over2 = 0;
                              for ( int s = Max( o, 0 ); s < y.isize( ); s++ )
                                   over2 += unibases[ y[s] ].isize( );
                              if ( over2 < min_overlap ) continue;
                              int new_over = Min( over, over2 );
                              map< vec<int>, int >&
                                   target = ( terminal ? exts_done : exts );
                              target[ynew] 
                                   = Max( target[ynew], new_over );    }    }    }

               vec< pair< int, vec<int> > > exts_donev;
               for ( emapit r = exts_done.begin( ); r != exts_done.end( ); r++ )
                    exts_donev.push( r->second, r->first );
               ReverseSort(exts_donev);
               if (VERBOSE)
                    cout << "\nextensions of " << u << "." << zz-i+1 << ":\n";
               for ( int l = 0; l < exts_donev.isize( ); l++ )
               {    if ( exts_donev[l].first 
                         < exts_donev[0].first - max_overlap_diff )
                    {    break;     }
                    const vec<int>& x = exts_donev[l].second;
                    int places0 = PLACES.size( );
                    ReportExtension( l, x, K, unibases, VALIDATION, LG, genome2, 
                         Glocs, data_dir, PLACES, all, exts_all, VERBOSE );    
                    int nplaces = PLACES.isize( ) - places0;
                    int elide_count = 0, nlinks = 0, v = x.back( ), dist = 0;
                    for ( set< vec<int> >::iterator ei = elides.begin( );
                         ei != elides.end( ); ei++ )
                    {    const vec<int>& e = *ei;
                         if ( x.Contains(e) )
                         {    if (VERBOSE)
                              {    cout << "Warning: contains elided repeat";
                                   for ( int m = 0; m < e.isize( ); m++ )
                                        cout << " " << e[m];
                                   cout << "\n";    }
                              elide_count++;    }    }
                    for ( int m = 0; m < H.From(u).isize( ); m++ )
                    {    if ( H.From(u)[m] == v )
                              nlinks = H.EdgeObjectByIndexFrom( u, m ).nlinks;    }
                    for ( int m = 1; m < x.isize( ) - 1; m++ )
                         dist += unibases[ x[m] ].isize( ) - (K-1);
                    if (VERBOSE) PRINT3( elide_count, dist, nlinks );    
                    if ( VALIDATION >= 1 && nlinks == 0 && nplaces > 0 )
                    {    cout << "Ouch!  Placed without links.\n";    }    }    }
          i = j - 1;    }

     // Build a snark.

     uniseq dummy;
     dummy.SetUnibases(unibases);
     int ne = exts_all.size( );
     for ( int i = 0; i < ne; i++ )
     {    vec<int> x = exts_all[i];
          x.ReverseMe( );
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = to_rc[ x[j] ];
          exts_all.push_back(x);    }
     UniqueSort(exts_all);
     vec<int> longs;
     vec<uniseq> seq;
     for ( int u = 0; u < (int) unibases.size( ); u++ )
     {    if ( unibases[u].isize( ) >= min_seed ) 
          {    longs.push_back(u);
               vec<int> us(1), overlap;
               us[0] = u;
               seq.push( us, overlap );    }    }
     int N = seq.size( );
     vec< triple< int, int, vec<int> > > exts;
     for ( int i = 0; i < exts_all.isize( ); i++ )
     {    const vec<int>& x = exts_all[i];
          exts.push( x.front( ), x.back( ), x );    }
     Sort(exts);
     vec< vec<int> > from(N), to(N), from_edge_obj(N), to_edge_obj(N);
     vec<gapster> edges;
     for ( int i = 0; i < exts.isize( ); i++ )
     {    int u1 = exts[i].first, u2 = exts[i].second;
          int pu1 = Position(longs, u1), pu2 = Position(longs, u2);
          int j;
          for ( j = i + 1; j < exts.isize( ); j++ )
               if ( exts[j].first != u1 || exts[j].second != u2 ) break;
          vec<uniseq> c;
          for ( int k = i; k < j; k++ )
          {    const vec<int>& x = exts_all[k];
               vec<int> over( x.isize( ) - 1, K - 1 );
               c.push_back( uniseq( x, over ) );    }
          from[pu1].push_back(pu2), to[pu2].push_back(pu1);
          from_edge_obj[pu1].push_back( edges.size( ) );
          to_edge_obj[pu2].push_back( edges.size( ) );
          edges.push(c);
          i = j - 1;    }
     for ( int i = 0; i < seq.isize( ); i++ )
     {    SortSync( from[i], from_edge_obj[i] );
          SortSync( to[i], to_edge_obj[i] );    }
     digraphE<gapster> G( from, to, edges, to_edge_obj, from_edge_obj );
     snark S( G, seq );
     S.SetUnibases(unibases);

     // Look for vertices that could be pulled apart (to do).

     /*
     cout << "\n" << Date( ) << ": pulling apart vertices" << endl;
     for ( int vp = 0; vp < S.VertN( ); vp++ )
     {    if ( S.From(vp).isize( ) == 2 && S.To(vp).isize( ) == 2 )
          {    int v = S.Vert(vp).U(0);
               int x1p = S.To(vp)[0], x2p = S.To(vp)[1];
               int y1p = S.From(vp)[0], y2p = S.From(vp)[1];
               int x1 = S.Vert(x1p).U(0);
               int x2 = S.Vert(x2p).U(0);
               int y1 = S.Vert(y1p).U(0);
               int y2 = S.Vert(y2p).U(0);
               if ( x1 == v || x2 == v || y1 == v || y2 == v ) continue;
               cout << "\ncout try " << v << endl;
               PRINT4( x1, x2, y1, y2 );
               for ( int j = 0; j < H.From(x1).isize( ); j++ )
               {    int z = H.From(x1)[j];
                    if ( z == y1 )
                    {    int nlinks = H.EdgeObjectByIndexFrom( x1, j ).nlinks;
                         cout << "see link from x1 to y1, count = " << nlinks
                              << endl;    }
                    if ( z == y2 )
                    {    int nlinks = H.EdgeObjectByIndexFrom( x1, j ).nlinks;
                         cout << "see link from x1 to y2, count = " << nlinks
                              << endl;    }    }
               for ( int j = 0; j < H.From(x2).isize( ); j++ )
               {    int z = H.From(x2)[j];
                    if ( z == y1 )
                    {    int nlinks = H.EdgeObjectByIndexFrom( x2, j ).nlinks;
                         cout << "see link from x2 to y1, count = " << nlinks
                              << endl;    }
                    if ( z == y2 )
                    {    int nlinks = H.EdgeObjectByIndexFrom( x2, j ).nlinks;
                         cout << "see link from x2 to y2, count = " << nlinks
                              << endl;    }    }    }    }
     */

     // Clean up.

     S.SwallowSimpleGaps( );
     S.BringOutTheDead( );
     if (VERBOSE) cout << "\n";
     cout << Date( ) << ": assembly has " << S.EdgeN( ) << " edges" << endl;

     // Generate dot file.

     if ( DOT != "" )
     {    Ofstream( dout, DOT );
          vec< vec<String> > legends;
          vec<String> vertex_labels, legend1;
          legend1.push_back( "<FONT POINT-SIZE=\"12\">Legend</FONT>" );
          int count = 0;
          for ( int x = 0; x < S.VertN( ); x++ )
          {    if ( S.Vert(x).N( ) == 1 )
                    vertex_labels.push_back( ToString( S.Vert(x).U(0) ) );
               else
               {    String s = BaseAlpha(count) + " = ";
                    int page_width = 120;
                    for ( int j = 0; j < S.Vert(x).U( ).isize( ); j++ )
                    {    if ( j > 0 ) 
                         {    int over = S.Vert(x).Over(j-1);
                              if ( over == K - 1 ) s += "+";
                              else s += "--(" + ToString(over) + ")--";     }
                         s += ToString( S.Vert(x).U(j) );
                         if ( j < S.Vert(x).N( ) - 1 && s.isize( ) >= page_width )
                         {    legend1.push_back(s);
                              s = "    ";    }    }
                    vertex_labels.push_back( BaseAlpha(count++) );
                    legend1.push_back(s);    }    }
          legends.push_back(legend1);
          vec<String> colors;
          vec< vec<String> > edge_labels( S.VertN( ) );
          for ( int v = 0; v < S.VertN( ); v++ )
          {    for ( int j = 0; j < S.G( ).From(v).isize( ); j++ )
               {   int n = S.G( ).EdgeObjectByIndexFrom( v, j ).ClosureCount( );
                    edge_labels[v].push_back( ToString(n) );    }    }
          S.G( ).DOT_vl( dout, vertex_labels, ( CIRCO ? "circo" : "" ), 
               legends, colors, edge_labels );    }

     // Generate report.

     if ( VALIDATION >= 1 )
     {    Sort(PLACES);
          cout << "\n==========================================================="
               << "=========================\n\n";
          cout << "covered parts of genome:\n";
          for ( int i = 0; i < PLACES.isize( ); i++ )
          {    int j = i + 1;
               placementy& p = PLACES[i];
               while(1)
               {    if ( PLACES[j].g != p.g ) break;
                    if ( PLACES[j-1].Pos - PLACES[j].pos < min_seed ) break;
                    if ( ++j == PLACES.isize( ) ) break;    }
               cout << p.g << ".";
               p.Pos = PLACES[j-1].Pos;
               int ng = genome2[p.g].isize( ) / 2;
               if ( p.Pos >= ng && p.Pos - ng >= p.pos + min_seed ) 
               {    cout << "all (circ)\n";    }
               else if ( p.Pos >= ng )
               {    cout << p.pos << "-" << p.Pos - ng << " (circ)\n";    }
               else cout << p.pos << "-" << p.Pos << "\n";
               i = j - 1;    }    }

     // Write output.
               
     if (WRITE) 
     {    cout << Date( ) << ": building output" << endl;
          if ( !USE_ALL ) all.resize(0);
          if (USE_UNI) all.Append(unibases);
          for ( int i = 0; i < seqs.isize( ); i++ )
          {    const vec<int>& x = seqs[i];
               basevector b;
               for ( int m = 0; m < x.isize( ); m++ )
               {    b = Cat( b, unibases[ x[m] ] );
                    if ( m < x.isize( ) - 1 ) b.resize( b.isize( ) - (K-1) );    }
               all.push_back_reserve(b);    }
          // all.WriteAll( "temp.fastb" );
          vecKmerPath newpaths, newpathsrc, newunipaths;
          vec<tagged_rpint> newpathsdb, newunipathsdb;
          ReadsToPathsCoreY( all, KBIG, newpaths, newpathsrc, newpathsdb,
               run_dir + "/LRJ", NUM_THREADS );
          Unipath( newpaths, newpathsrc, newpathsdb, newunipaths, newunipathsdb );
	  digraph A;
	  BuildUnipathAdjacencyGraph( newpaths, newpathsrc, newpathsdb, 
               newunipaths, newunipathsdb, A);
          KmerBaseBroker newkbb( KBIG, newpaths, newpathsrc, newpathsdb, all );
          vecbasevector newunibases;
          for ( size_t i = 0; i < newunipaths.size( ); i++ )
               newunibases.push_back_reserve( newkbb.Seq( newunipaths[i] ) );
          cout << Date( ) << ": writing output files" << endl;
          String KBIGS = ToString(KBIG);
          newunipaths.WriteAll( outhead + ".unipaths.k" + KBIGS );
          BinaryWrite3( outhead + ".unipathsdb.k" + KBIGS, newunipathsdb );
	  BinaryWrite( outhead + ".unipath_adjgraph.k" + KBIGS, A );
          newunibases.WriteAll( outhead + ".unibases.k" + KBIGS );    }
     cout << Date( ) << ": done" << endl;    }
