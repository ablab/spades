///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ConnectScaffolds.

#include "Equiv.h"
#include "MainTools.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/ReadLoc.h"
#include "system/ParsedArgs.h"

class olink {
     public:

     olink( ) { }
     olink( const int s1, const int s2, const Bool rc1, const Bool rc2,
          const int overlap, const double error_rate )
          : s1(s1), s2(s2), rc1(rc1), rc2(rc2), overlap(overlap),
          error_rate(error_rate) { }

     friend Bool operator<( const olink& o1, const olink& o2 )
     {    if ( o1.s1 < o2.s1 ) return True;
          if ( o1.s1 > o2.s1 ) return False;
          if ( o1.s2 < o2.s2 ) return True;
          if ( o1.s2 > o2.s2 ) return False;
          if ( o1.rc1 < o2.rc1 ) return True;
          if ( o1.rc1 > o2.rc1 ) return False;
          if ( o1.rc2 < o2.rc2 ) return True;
          if ( o1.rc2 > o2.rc2 ) return False;
          if ( o1.overlap < o2.overlap ) return True;
          if ( o1.overlap > o2.overlap ) return False;
          return o1.error_rate < o2.error_rate;    }

     friend Bool operator==( const olink& o1, const olink& o2 )
     {    return o1.s1 == o2.s1 && o1.s2 == o2.s2 && o1.rc1 == o2.rc1
               && o1.rc2 == o2.rc2 && o1.overlap == o2.overlap
               && o1.error_rate == o2.error_rate;    }

     int s1;
     int s2;
     Bool rc1;
     Bool rc2;
     int overlap;
     double error_rate;

};

class Link {
     public:

     Link( ) { }
     Link( const int m1, const int s1, const int m2, const int s2, const int sep, 
          const int dev, const int read_class, const int dist_to_end1, 
          const int dist_to_end2, const Bool rc1, const Bool rc2 ) : m1(m1), s1(s1),
          m2(m2), s2(s2), sep(sep), dev(dev), read_class(read_class), 
          dist_to_end1(dist_to_end1), dist_to_end2(dist_to_end2), rc1(rc1), 
          rc2(rc2) { }

     friend Bool operator<( const Link& l1, const Link& l2 )
     {    return l1.s2 < l2.s2;     }

     int m1;
     int s1;
     int m2;
     int s2;
     int sep;
     int dev;
     int read_class;
     int dist_to_end1;
     int dist_to_end2;
     Bool rc1;
     Bool rc2;
     int overlap;
     double error_rate;

};

int main(int argc, char *argv[])
{
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
     CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
          "if specified, file extension is .READLOCS_PREFIX.readlocs "
          "instead of .readlocs");
     CommandArgument_Int_OrDefault(MAX_OVERLAP, 20000);
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Bool_OrDefault(MERGE_SCAFFOLDS, True);
     CommandArgument_Bool_OrDefault(MERGE_CONTIGS, True);
     EndCommandArguments;

     // Define directories.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

     // Load fastb.

     vecbasevector tigs( sub_dir + "/" + ASSEMBLY + ".contigs.fastb" );
     vec<fastavector> tigsa;
     LoadFromFastaFile( sub_dir + "/" + ASSEMBLY + ".contigs.fasta", tigsa );
     int ntigs = tigs.size( );

     // Load scaffolds.

     vec<superb> scaffolds;
     ReadSuperbs( sub_dir + "/" + ASSEMBLY + ".superb", scaffolds );
     vec<int> to_super( ntigs, -1 ), to_super_pos( ntigs, -1 );
     vec<int> to_super_posr( ntigs, - 1 );
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    int n = scaffolds[i].Ntigs( );
          for ( int j = 0; j < n; j++ )
          {    to_super[ scaffolds[i].Tig(j) ] = i;
               to_super_pos[ scaffolds[i].Tig(j) ] = j;   
               to_super_posr[ scaffolds[i].Tig(j) ] = n - j - 1;    }    }

     // Set up to fetch read locations.

     String head = sub_dir + "/" + ASSEMBLY;
     if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
     read_locs_on_disk locs_file( head, run_dir );

     // Create data structures to record what is to be deleted.

     vec<Bool> tigs_to_delete( tigs.size( ), False );
     vec<Bool> scaffolds_to_delete( scaffolds.size( ), False );

     // Go through the scaffolds.

     const double max_error_rate = 0.05;
     if (MERGE_SCAFFOLDS)
     {    const int bandwidth = 100;
          vec<olink> olinks;
          for ( int s1 = 0; s1 < scaffolds.isize( ); s1++ )
          {    const int min_length = 20000;
               const superb& S1 = scaffolds[s1];
               if ( S1.FullLength( ) < min_length ) continue;
               if (VERBOSE) cout << "\nlooking at scaffold " << s1 << endl;
     
               vec<Link> Links;
               for ( int p1 = 0; p1 < S1.Ntigs( ); p1++ )
               {    int m1 = S1.Tig(p1);
                    vec<read_loc> locs;
                    locs_file.LoadContig( m1, locs );
                    for ( int j = 0; j < locs.isize( ); j++ )
                    {    const read_loc& rl = locs[j];
                         if ( !rl.PartnerPlaced( ) ) continue;
                         int x1 = rl.Start( );
                         int y1 = scaffolds[s1].SubSuperLength( 0, p1 - 1 ) + x1;
                         int m2 = rl.PartnerContigId( );
                         int s2 = to_super[m2], p2 = to_super_pos[m2];
                         const superb& S2 = scaffolds[s2];
                         if ( s2 == s1 ) continue;
                         int x2 = rl.PartnerStart( );
                         int y2 = scaffolds[s2].SubSuperLength( 0, p2 - 1 ) + x2;
          
                         // For now just consider 1/2 of the cases:
     
                         if ( rl.Fw( ) && rl.PartnerRc( ) )
                         {    int dist_to_end1 = S1.Len(p1) - rl.Stop( ) 
                                   + S1.SubSuperLength( p1 + 1, S1.Ntigs( ) - 1 );
                              int dist_to_end2 = rl.PartnerStart( ) 
                                   + S2.SubSuperLength( 0, p2 - 1 );
                              int sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
                              int dev = rl.Dev( );
                              if ( sep < -MAX_OVERLAP ) continue;
                              Links.push( m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
                                   dist_to_end1, dist_to_end2, False, False );    }
                         if ( rl.Fw( ) && rl.PartnerFw( ) )
                         {    int dist_to_end1 = S1.Len(p1) - rl.Stop( ) 
                                   + S1.SubSuperLength( p1 + 1, S1.Ntigs( ) - 1 );
                              int dist_to_end2 = S2.Len(p2) - rl.PartnerStop( ) 
                                   + S2.SubSuperLength( p2 + 1, S2.Ntigs( ) - 1 );
                              int sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
                              int dev = rl.Dev( );
                              if ( sep < -MAX_OVERLAP ) continue;
                              Links.push( m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
                                   dist_to_end1, dist_to_end2, False, True );    }
                         if ( rl.Rc( ) && rl.PartnerRc( ) )
                         {    int dist_to_end1 
                                   = rl.Start( ) + S1.SubSuperLength( 0, p1 - 1 );
                              int dist_to_end2 = rl.PartnerStart( ) 
                                   + S2.SubSuperLength( 0, p2 - 1 );
                              int sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
                              int dev = rl.Dev( );
                              if ( sep < -MAX_OVERLAP ) continue;
                              Links.push( m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
                                   dist_to_end1, dist_to_end2, 
                                   True, False );   }   }   }

               Sort(Links);
               for ( int i = 0; i < Links.isize( ); i++ )
               {    const Link& L = Links[i];
                    int j;
                    for ( j = i + 1; j < Links.isize( ); j++ )
                         if ( Links[j].s2 != Links[i].s2 ) break;
                    equiv_rel e( j - i );
                    for ( int k1 = i; k1 < j; k1++ )
                    {    for ( int k2 = i + 1; k2 < j; k2++ )
                         {    int delta = Abs( Links[k1].dist_to_end1 
                                   - Links[k2].dist_to_end1 )
                                   + Abs( Links[k1].dist_to_end2 
                                   - Links[k2].dist_to_end2 );
                              int min_delta = 10;
                              if ( delta < min_delta 
                                   && Links[k1].rc1 == Links[k2].rc1 
                                   && Links[k1].rc2 == Links[k2].rc2 ) 
                              {    e.Join( k1 - i, k2 - i );    }    }    }
                    vec<int> reps;
                    e.OrbitReps(reps);
                    if ( reps.size( ) > 1 )
                    {    if (VERBOSE) cout << "\n";
                         for ( int r = 0; r < reps.isize( ); r++ )
                         {    int k = i + reps[r];
                              const Link& L = Links[k];
                              int s2 = L.s2, sep = L.sep, dev = L.dev;
                              int m1 = L.m1, m2 = L.m2;
                              int dist_to_end1 = L.dist_to_end1;
                              int dist_to_end2 = L.dist_to_end2;
                              String dist_to_end = ToString(dist_to_end1) + "/"
                                   + ToString(dist_to_end2);
     
                              // Try aligning.  Doesn't belong here.
     
                              const int min_sep = -1000;
                              int offset = tigs[m1].size( ) + sep;
                              if ( sep < min_sep )
                              {    align a;
                                   int errors;
                                   basevector T1 = tigs[m1], T2 = tigs[m2];
                                   if (L.rc1) T1.ReverseComplement( );
                                   if (L.rc2) T2.ReverseComplement( );
                                   if ( offset > T1.isize( ) 
                                        || offset < -T2.isize( ) )
                                   {    continue;    }
                                   int score = SmithWatBandedA( T1, T2,
                                        offset, bandwidth, a, errors, 0, 1, 1 );
                                   if ( a.pos2( ) > 0 ) continue;
                                   int overlap = a.Pos1( ) - a.pos1( );
               
                                   // Require low error rate.
     
                                   double error_rate = double(score)/double(overlap);
                                   if ( error_rate > max_error_rate ) continue;
     
                                   // Require at ends.
     
                                   int p1 = to_super_pos[m1], p2 = to_super_pos[m2];
                                   const superb& S1 = scaffolds[s1];
                                   const superb& S2 = scaffolds[s2];
                                   if ( !L.rc1 && S1.Ntigs( ) - p1 != 1 ) continue;
                                   if ( L.rc1 && p1 != 0 ) continue;
                                   if ( L.rc2 && S2.Ntigs( ) - p2 != 1 ) continue;
                                   if ( !L.rc2 && p2 != 0 ) continue;
                                        
                                   // Print.
     
                                   if (VERBOSE)
                                   {    cout << "p1 = " << p1+1 << ", "
                                             << scaffolds[s1].Ntigs( ) - p1 << "; ";
                                        cout << "p2 = " << p2+1 << ", "
                                             << scaffolds[s2].Ntigs( ) - p2 << "; ";
                                        if ( L.read_class == 0 ) cout << "frag ";
                                        if ( L.read_class == 1 ) cout << "jump ";
                                        if ( L.read_class == 2 ) cout << "long ";
                                        cout << "s" << s1 
                                             << ( L.rc1 ? ".rc" : ".fw" ) 
                                             << "[" << m1 << "] --> ";
                                        cout << "s" << s2 
                                             << ( L.rc2 ? ".rc" : ".fw" ) 
                                             << "[" << m2 << "], ";
                                        cout << "sep: " << sep << " +/- " << dev 
                                             << ", ->end: " << dist_to_end << "\n";
                                        cout << "score = " << score << ", overlap = "
                                             << overlap << ", error rate = " 
                                             << 100.0 * error_rate << "\n";
                                        PRINT2( a.pos1( ), a.pos2( ) );    }
                                   error_rate = 100.0 * error_rate;
                                   if (VERBOSE)
                                   {    cout << "WOOF ";
                                        cout << "s" << s1 << ( L.rc1 ? "rc" : "fw" );
                                        cout << " --> ";
                                        cout << "s" << s2 << ( L.rc2 ? "rc" : "fw" );
                                        cout << ", ";
                                        PRINT2( overlap, error_rate );    }
                                   olinks.push( s1, s2, L.rc1, L.rc2, 
                                        overlap, error_rate );    }    }    }
                    i = j - 1;    }    }

          // Build graph from links.

          UniqueSort(olinks);
          int N = 2 * scaffolds.size( );
          vec< vec<int> > from(N), to(N), from_edge_obj(N), to_edge_obj(N);
          vec<int> edges;
          for ( int j = 0; j < olinks.isize( ); j++ )
          {    int k;
               for ( k = j+1; k < olinks.isize( ); k++ )
               {    if ( olinks[k].s1 != olinks[j].s1 ) break;
                    if ( olinks[k].s2 != olinks[j].s2 ) break;
                    if ( olinks[k].rc1 != olinks[j].rc1 ) break;
                    if ( olinks[k].rc2 != olinks[j].rc2 ) break;    }
               if (VERBOSE) cout << "\n";
               Bool too_far = False;
               const int max_dist = 100;
               for ( int u1 = j; u1 < k; u1++ )
               {    for ( int u2 = u1+1; u2 < k; u2++ )
                    {    if ( Abs( olinks[u1].overlap - olinks[u2].overlap ) 
                              > max_dist )
                         {    too_far = True;    }    }    }
               if ( !too_far )
               {    double min_error_rate = 100.0;
                    int best = -1;
                    for ( int u = j; u < k; u++ )
                    {    if ( olinks[u].error_rate < min_error_rate )
                         {    min_error_rate = olinks[u].error_rate;
                              best = u;    }    }
                    ForceAssertGe( best, 0 );
                    const olink& o = olinks[best];
                    if (VERBOSE)
                    {    cout << "SNORT ";
                         cout << "s" << o.s1 << ( o.rc1 ? "rc" : "fw" );
                         cout << " --> ";
                         cout << "s" << o.s2 << ( o.rc2 ? "rc" : "fw" );
                         cout << ", ";
                         PRINT2( o.overlap, o.error_rate );    }
                    int id1 = 2*o.s1 + ( o.rc1 ? 1 : 0 );
                    int id2 = 2*o.s2 + ( o.rc2 ? 1 : 0 );
                    if ( !Member( from[id1], id2 ) )
                    {    from[id1].push_back(id2);
                         to[id2].push_back(id1);
                         from_edge_obj[id1].push_back( edges.size( ) );
                         to_edge_obj[id2].push_back( edges.size( ) );
                         edges.push_back( o.overlap );    }
                    int id1p = 2*o.s2 + ( o.rc2 ? 0 : 1 );
                    int id2p = 2*o.s1 + ( o.rc1 ? 0 : 1 );
                    if ( !Member( from[id1p], id2p ) )
                    {    from[id1p].push_back(id2p);
                         to[id2p].push_back(id1p);
                         from_edge_obj[id1p].push_back( edges.size( ) );
                         to_edge_obj[id2p].push_back( edges.size( ) );
                         edges.push_back( o.overlap );    }    }
               j = k - 1;    }
          for ( int i = 0; i < N; i++ )
          {    Sort( from[i] );
               Sort( to[i] );
               Sort( from_edge_obj[i] );
               Sort( to_edge_obj[i] );    }
          digraphE<int> G( from, to, edges, to_edge_obj, from_edge_obj );

          // Find lines in the graph.

          vec<int> used( N, False );
          vec< vec<int> > lines;
          for ( int v = 0; v < N; v++ )
          {    int vp = ( v % 2 == 0 ? v + 1 : v - 1 );
               if ( used[v] || used[vp] ) continue;
               used[v] = True;
               vec<int> line;
               line.push_back(v);
               int w = v;
               while(1)
               {    if ( !G.To(w).solo( ) ) break;
                    int z = G.To(w)[0];
                    int zp = ( z % 2 == 0 ? z + 1 : z - 1 );
                    if ( used[z] || used[zp] ) break;
                    line.push_front(z);
                    used[z] = True;
                    w = z;    }
               w = v;
               while(1)
               {    if ( !G.From(w).solo( ) ) break;
                    int z = G.From(w)[0];
                    int zp = ( z % 2 == 0 ? z + 1 : z - 1 );
                    if ( used[z] || used[zp] ) break;
                    line.push_back(z);
                    used[z] = True;
                    w = z;    }
               if (VERBOSE)
               {    cout << "LINE:";
                    for ( int j = 0; j < line.isize( ); j++ )
                    {    cout << " " << line[j]/2 
                              << ( line[j] % 2 == 0 ? "fw" : "rc" );    }
                    cout << "\n";    }
               lines.push_back(line);    }

          // Merge scaffolds.

          for ( int i = 0; i < lines.isize( ); i++ )
          {    const vec<int>& L = lines[i];
               // cout << "reversing" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               for ( int j = 0; j < L.isize( ); j++ )
               {    int s = L[j]/2;
                    if ( L[j] % 2 == 1 )
                    {    scaffolds[s].Reverse( );
                         for ( int k = 0; k < scaffolds[s].Ntigs( ); k++ )
                         {    tigs[ scaffolds[s].Tig(k) ].ReverseComplement( );
                              tigsa[ scaffolds[s].Tig(k) ].
                                   ReverseComplement( );    }    }    }
               // PRINT(i); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               int s1 = L[0]/2;
               superb& S1 = scaffolds[s1];
               for ( int r = 1; r < L.isize( ); r++ )
               {    int s2 = L[r]/2;
                    // PRINT(r); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    superb& S2 = scaffolds[s2];
                    int m1 = S1.Tig( S1.Ntigs( ) - 1 ), m2 = S2.Tig(0);
                    int o = G.EdgeObjectByIndexFrom( L[r-1], 0 );
     
                    // Merge m1 to m2 along overlap o.
     
                    // cout << "aligning" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    align a;
                    int errors;
                    int offset = tigs[m1].isize( ) - o;
                    // PRINT2( m1, m2 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    // PRINT2( tigs[m1].size( ), tigs[m2].size( ) ); // XXXXXXXXXXXX
                    int score = SmithWatBandedA( tigs[m1], tigs[m2],
                         offset, bandwidth, a, errors, 0, 1, 1 );
                    // PRINT(score); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    // PRINT2( s1, s2 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    // Note that we're not checking that a.pos2( ) == 0, and probably
                    // should.  I believe that it might assert.
                    basevector newtig;
                    fastavector newtiga;
                    int overhang;
                    // PRINT4( a.pos1( ), a.Pos1( ), a.pos2( ), a.Pos2( ) ); // XXXX
                    // PRINT2( a.extent1( ), tigs[m2].size( ) ); // XXXXXXXXXXXXXXXX
                    if ( a.extent1( ) <= tigs[m2].isize( ) )
                    {    overhang = 0;
                         newtig = Cat( basevector( tigs[m1], 0, a.pos1( ) ), 
                              tigs[m2] );
                         fastavector n1;
                         n1.SetToSubOf( tigsa[m1], 0, a.pos1( ) );
                         newtiga = Cat( n1, tigsa[m2] );    }
                    else
                    {    newtig = tigs[m1];
                         newtiga = tigsa[m1];
                         overhang = a.extent1( ) - tigs[m2].isize( );    }
                    // cout << "subtracting overhang" << endl; // XXXXXXXXXXXXXXXXXX
                    if ( S2.Ntigs( ) > 1 ) 
                         S2.SetGap(0, scaffolds[s2].Gap(0) - overhang);
                    // PRINT( newtig.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    tigs[m1] = newtig, tigs[m2].resize(0);
                    tigsa[m1] = newtiga, tigsa[m2].resize(0);
                    tigs_to_delete[m2] = True;
                    S1.SetLen( S1.Ntigs( ) - 1, tigs[m1].size( ) );
                    // cout << "appending tigs" << endl; // XXXXXXXXXXXXXXXXXXXXXXXX
                    for ( int j = 1; j < S2.Ntigs( ); j++ )
                    {    scaffolds[s1].AppendTig( S2.Tig(j), 
                              S2.Len(j), S2.Gap(j-1), S2.Dev(j-1) );    }
                    // cout << "clearing" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    S2.Clear( );
                    scaffolds_to_delete[s2] = True;    }    }    }

     // Merge adjacent overlapping contigs within scaffolds.

     if (MERGE_CONTIGS)
     {    double contig_merge_clock = WallClockTime( );
          int contig_merges = 0;
          const int bandwidth = 1500;
          const int min_overlap_to_close = 500;
          const double max_delta = 0.2;
          for ( int s = 0; s < scaffolds.isize( ); s++ )
          {    if ( scaffolds_to_delete[s] ) continue;
               superb& S = scaffolds[s];
               for ( int j = 0; j < S.Ntigs( ) - 1; j++ )
               {    if ( S.Gap(j) > -min_overlap_to_close ) continue;
                    int m1 = S.Tig(j), m2 = S.Tig(j+1);
                    const basevector &T1 = tigs[m1], &T2 = tigs[m2];
                    int offset = T1.isize( ) + S.Gap(j);
                    align a;
                    int errors;
                    if ( offset > T1.isize( ) || offset < -T2.isize( ) ) continue;
                    int score = SmithWatBandedA( T1, T2, offset, bandwidth, 
                         a, errors, 0, 1, 1 );
                    int overlap = a.Pos1( ) - a.pos1( );
                    double error_rate = double(score)/double(overlap);
                    if ( error_rate > max_error_rate ) continue;
                    int actual_offset = a.Pos2( ) + T1.isize( ) - a.Pos1( );
                    double delta = double( Abs( actual_offset + S.Gap(j) ) ) 
                         / double( -S.Gap(j) );
                    if ( delta > max_delta ) continue;
                    vec<ho_interval> perf1, perf2;
                    a.PerfectIntervals1( T1, T2, perf1 );
                    a.PerfectIntervals2( T1, T2, perf2 );
                    int best = -1, p = -1;
                    for ( int z = 0; z < perf1.isize( ); z++ )
                    {    if ( perf1[z].Length( ) > p )
                         {    p = perf1[z].Length( );
                              best = z;    }    }
                    if ( best < 0 ) continue; // not sure how this can happen
                    const ho_interval &P1 = perf1[best], &P2 = perf2[best];
                    basevector newtig;
                    fastavector newtiga;
                    int del1 = 0, del2 = 0;
                    int n1 = T1.size( ), n2 = T2.size( );
                    basevector left, middle, right;
                    fastavector lefta, middlea, righta;
                    middle.SetToSubOf( T1, P1.Start( ), P1.Length( ) );
                    middlea.SetToSubOf( tigsa[m1], P1.Start( ), P1.Length( ) );
                    if ( P1.Start( ) >= P2.Start( ) )
                    {    left.SetToSubOf( T1, 0, P1.Start( ) );
                         lefta.SetToSubOf( tigsa[m1], 0, P1.Start( ) );    }
                    else
                    {    left.SetToSubOf( T2, 0, P2.Start( ) );
                         lefta.SetToSubOf( tigsa[m2], 0, P2.Start( ) );
                         del1 = P2.Start( ) - P1.Start( );    }
                    if ( n1 - P1.Stop( ) >= n2 - P2.Stop( ) )
                    {    right.SetToSubOf( T1, P1.Stop( ), n1 - P1.Stop( ) );
                         righta.SetToSubOf( tigsa[m1], P1.Stop( ), n1 - P1.Stop( ) );
                         del2 = ( n1 - P1.Stop( ) ) - ( n2 - P2.Stop( ) );    }
                    else
                    {    right.SetToSubOf( T2, P2.Stop( ), n2 - P2.Stop( ) );
                         righta.SetToSubOf( tigsa[m2], P2.Stop( ), n2 - P2.Stop( ) );
                              }
                    newtig = Cat( left, middle, right );
                    newtiga = Cat( lefta, middlea, righta );

                    /*
                    fastavector n1;
                    if ( a.pos2( ) == 0 && a.Pos1( ) == T1.isize( ) )
                    {    newtig = Cat( basevector( tigs[m1], 0, a.pos1( ) ), 
                              tigs[m2] );
                         n1.SetToSubOf( tigsa[m1], 0, a.pos1( ) );
                         newtiga = Cat( n1, tigsa[m2] );    }
                    else if ( a.pos2( ) == 0 && a.Pos2( ) == T2.isize( ) )
                    {    newtig = tigs[m1];
                         newtiga = tigsa[m1];
                         del2 = T1.isize( ) - a.Pos1( );    }
                    else if ( a.pos1( ) == 0 && a.Pos1( ) == T1.isize( ) )
                    {    newtig = tigs[m2];
                         newtiga = tigsa[m2];
                         del1 = a.pos2( );    }
                    else if ( a.pos1( ) == 0 && a.Pos2( ) == T2.isize( ) )
                    {    del2 = tigs[m1].isize( ) - a.Pos1( );
                         newtig = Cat( tigs[m2], basevector( tigs[m1], a.Pos1( ),
                              del2 ) );
                         n1.SetToSubOf( tigsa[m1], a.Pos1( ), del2 );
                         newtiga = Cat( tigsa[m2], n1 );
                         del1 = a.pos2( );    }
                    else continue;
                    */

                    if (VERBOSE)
                    {    cout << "\n";
                         PRINT6( a.pos1( ), a.Pos1( ), T1.size( ),
                                 a.pos2( ), a.Pos2( ), T2.size( ) );
                         cout << "s = " << s << ", m1 = " << m1 << ", m2 = " << m2
                              << ", overlap = " << overlap << ", errs = " 
                              << setprecision(3) << 100.0 * error_rate
                              << "%, delta = " << 100.0 * delta << "%\n";    }
                    contig_merges++;
                    tigs[m1] = newtig, tigs[m2].resize(0);
                    tigsa[m1] = newtiga, tigsa[m2].resize(0);
                    tigs_to_delete[m2] = True;
                    S.SetLen( j, newtig.size( ) );
                    S.RemoveTigByPosAlt( j+1, -del1, -del2 );
                    j--;    }    }
          cout << "\n" << contig_merges << " contig merges, time used = "
               << TimeSince(contig_merge_clock) << endl;    }
     
     // Clean up data structures.

     EraseIf( scaffolds, scaffolds_to_delete );
     vec<int> new_id( tigs.size( ) );
     size_t count = 0;
     for ( size_t j = 0; j < tigs.size( ); j++ )
     {    if ( !tigs_to_delete[j] )
          {    new_id[j] = count;
               if ( j > count ) 
               {    tigs[count] = tigs[j], tigsa[count] = tigsa[j];    }
               count++;    }    }
     tigs.resize(count);
     tigsa.resize(count);
     for ( int j = 0; j < scaffolds.isize( ); j++ )
     {    superb& S = scaffolds[j];
          for ( int l = 0; l < S.Ntigs( ); l++ )
               S.SetTig( l, new_id[ S.Tig(l) ] );    }

     // Write new assembly.

     if ( !WRITE ) return 0;
     vecfastavector tigsav;
     for ( size_t i = 0; i < tigsa.size( ); i++ )
          tigsav.push_back_reserve( tigsa[i] );
     tigsav.WriteAll( sub_dir + "/" + ASSEMBLY + ".connected.contigs.vecfasta" );
     tigs.WriteAll( sub_dir + "/" + ASSEMBLY + ".connected.contigs.fastb" );
     Ofstream( out, sub_dir + "/" + ASSEMBLY + ".connected.contigs.fasta" );
     for ( size_t i = 0; i < tigsa.size( ); i++ )
          tigsa[i].Print( out, ToString(i) );
     WriteSuperbs( sub_dir + "/" + ASSEMBLY + ".connected.superb", scaffolds );
     WriteSummary( sub_dir + "/" + ASSEMBLY + ".connected" , scaffolds );
     WriteScaffoldedFasta( sub_dir + "/" + ASSEMBLY + ".connected.assembly.fasta",
          tigsa, scaffolds );    }
