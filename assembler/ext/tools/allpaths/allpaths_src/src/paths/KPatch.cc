///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// KPatch.  Patch captured gaps in an efasta assembly.  Roughly, for each gap,
// take the last kmer on the left contig, and the first kmer on the right contig,
// and find all the paths between them in the unipath graph, having the 'right'
// size, and return these paths as alternative closures of the gap.  Gaps that
// could be negative require a slightly different treatment.
//
// Currently we filter patches in the following way: if the sum of the lengths of
// the 'variable' paths in a given patch is > 50, then we discard the patch.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Fastavector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "kmers/KmerRecord.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/AssemblyEdit.h"
#include "paths/GetNexts.h"

template<int K> void PatchMe( const vec<basevector>& jbases_sorted,
     const vec<int64_t>& jbases_sorted_id, const PairsManager& jpairs,
     const vec< triple<int64_t,int,int> >& jaligns, const vec<superb>& scaffolds,
     const vec<efasta>& tigse, const vecbasevector& tigs, 
     const vecbasevector& unibases, const vec< vec<int> >& nexts, 
     vec< vec<Bool> >& patched, vec< vec< triple<efasta,int,int> > >& results, 
     const int max_iterations, const int verbosity )
{
     // Sort the kmers in the unibases.

     cout << Date( ) << ": getting unibase " << K << "-mers" << endl;
     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + u.isize( ) - K + 1 );    }
     vec< triple<kmer<K>,int,int> > kmers_plus( starts.back( ) );
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               kmers_plus[r].first.SetToSubOf( u, j );
               kmers_plus[r].second = i, kmers_plus[r].third = j;    }    }

     cout << Date( ) << ": sorting " << kmers_plus.size() << " " << K << "-mers" << endl;

     sortInPlaceParallel(kmers_plus.begin(), kmers_plus.end());

     // Go through the gaps.  For each gap we store the outcome in "patched"
     // (boolean: patched or not) and "results" (patch, left_trim, right_trim),
     // where left_trim is the amount to be trimmed from the right end of the
     // left contig and right_trim is the amount to be trimmed from the left end
     // of the right contig, before inserting the patch.

     cout << Date( ) << ": going through the gaps" << endl;
     vec< pair<int,int> > gaps;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
     {    const superb& S = scaffolds[s];
          for ( int i = 0; i < S.Ngaps( ); i++ )
               gaps.push( s, i );     }
     patched.resize( scaffolds.size( ) ), results.resize( scaffolds.size( ) );
     for ( int s = 0; s < scaffolds.isize( ); s++ )
     {    patched[s].resize( scaffolds[s].Ngaps( ), False );
          results[s].resize( scaffolds[s].Ngaps( ) );    }
     vec<String> report( gaps.size( ) );
     #pragma omp parallel for
     for ( int gi = 0; gi < gaps.isize( ); gi++ )
     {    int s = gaps[gi].first, it = gaps[gi].second;
          const superb& S = scaffolds[s];
          int m1 = S.Tig(it), m2 = S.Tig(it+1), gap = S.Gap(it), dev = S.Dev(it);
          ostringstream out;
          if ( verbosity >= 2 ) 
               out << "\nlooking at gap from " << m1 << " to " << m2 << "\n";

          // Define bounds on walk size.

          const int dev_mult = 3;
          int min_walk = K + 1 + gap - dev_mult * dev;
          int max_walk = K + 1 + gap + dev_mult * dev;

          // Define the "end" kmers.  In the case of a gap that could be negative,
          // we can't take the actual ends, and instead have to back up from there.
          // For the moment we just move the right end forward, but we could move
          // both ends back equally, or something else.  (NOTE: CHANGED.)
          // Note also that we stay away from ambiguities.

          int n1 = tigse[m1].size( ), n2 = tigse[m2].size( );
          basevector b1(K), b2(K);
          int over = Max( 0, 1 - min_walk );

          // Find left.

          int stop1 = n1 - 1;
          find_left:
          int count1 = 0;
          for ( int j = stop1; j >= 0; j-- )
          {    if ( tigse[m1][j] == '}' )
               {    count1 = 0;
                    while( j >= 0 && tigse[m1][j] != '{' ) j--;
                    stop1 = j - 1;    }
               else 
               {    count1++;
                    if ( count1 == K ) break;    }    }
          if ( count1 < K )
          {    if ( verbosity >= 2 ) PRINT_TO( out, count1 );
               report[gi] = out.str( );
               continue;    }
          int over1 = n1 - 1 - stop1;
          for ( int j = 0; j < K; j++ )
               b1.Set( j, as_char( tigse[m1][ n1 - K + j - over1 ] ) );
          kmer<K> x1(b1);
          int64_t p1 = BinPosition1( kmers_plus, x1 );
          if ( p1 < 0 )
          {    stop1--;
               goto find_left;    }

          // Find right.

          int start2 = Max( 0, over - over1 );
          if ( start2 >= n2 )
          {    if ( verbosity >= 2 ) PRINT_TO( out, start2 );
               report[gi] = out.str( );
               continue;    }
          int brack_count = 0;
          for ( int j = 0; j < start2; j++ )
          {    if ( tigse[m2][j] == '{' ) brack_count++;
               if ( tigse[m2][j] == '}' ) brack_count--;    }
          for ( ; start2 < n2; start2++ )
          {    if ( brack_count == 0 ) break;
               if ( tigse[m2][start2] == '{' ) brack_count++;
               if ( tigse[m2][start2] == '}' ) brack_count--;    }
          find_right:
          int count2 = 0;
          for ( int j = start2; j < n2; j++ )
          {    if ( tigse[m2][j] == '{' ) 
               {    count2 = 0;
                    while ( j < n2 && tigse[m2][j] != '}' ) j++;
                    start2 = j + 1;    }
               else
               {    count2++;
                    if ( count2 == K ) break;    }    }
          if ( count2 < K )
          {    if ( verbosity >= 2 ) PRINT_TO( out, count2 );
               report[gi] = out.str( );
               continue;    }
          int over2 = start2;
          for ( int j = 0; j < K; j++ )
               b2.Set( j, as_char( tigse[m2][ over2 + j ] ) );
          kmer<K> x2(b2);
          int64_t p2 = BinPosition1( kmers_plus, x2 );
          if ( p2 < 0 )
          {    start2++;
               goto find_right;    }

          // Next.

          if ( verbosity >= 2 ) PRINT3_TO( out, K, over1, over2 );
          over = over1 + over2;
          min_walk += over, max_walk += over;
          int u1 = kmers_plus[p1].second, u2 = kmers_plus[p2].second;
          int s1 = kmers_plus[p1].third, s2 = kmers_plus[p2].third;
          pair<int,int> loc1 = make_pair(u1,s1), loc2 = make_pair(u2,s2);
          if ( verbosity >= 2 ) PRINT2_TO( out, u1, u2 );

          // Test for a case that we don't deal with.

          if ( u1 == u2 && !( s1 < s2 ) )
          {    if ( verbosity >= 2 ) PRINT4_TO( out, u1, u2, s1, s2 );
               report[gi] = out.str( );
               continue;    }

          // Try to walk through the unipath graph from one end to the other.

          vec< vec< pair<int,int> > > walks0, walks1;
          vec< pair<int,int> > winit;
          winit.push_back(loc1);
          walks0.push(winit);
          int bad = 0, iterations = 0;
          while( walks0.nonempty( ) )
          {    if ( ++iterations > max_iterations )
               {    bad = 1;
                    break;    }
               vec< pair<int,int> > w = walks0.back( );
               walks0.pop_back( );
               int n = 0;
               for ( int j = 1; j < w.isize( ); j++ )
                    n += unibases[ w[j-1].first ].isize( ) - w[j-1].second - K;
               if ( n > max_walk ) continue;
               int np = 
                    n + unibases[ w.back( ).first ].isize( ) - w.back( ).second - K;
               if ( np >= min_walk && n <= max_walk && w.back( ).first == u2 )
               {    if ( w.solo( ) ) w.push_back(loc2);
                    else w.back( ).second = s2;
                    walks1.push_back(w);    }
               int u = w.back( ).first, s = w.back( ).second;
               for ( int j = 0; j < nexts[u].isize( ); j++ )
               {    int v = nexts[u][j];
                    vec< pair<int,int> > wnew(w);
                    wnew.push( v, 0 );
                    walks0.push_back(wnew);    }    }
          if ( verbosity >= 2 ) PRINT3_TO( out, walks1.size( ), iterations, bad );
          if ( walks1.empty( ) || bad == 1 )
          {    report[gi] = out.str( );
               continue;    }

          // Generate the sequences associated to the patches.  These extend from
          // the right end of the left contig to the right end of the "end" kmer on 
          // the right contig.

          int npatches = walks1.size( );
          vec<basevector> patches(npatches);
          for ( int j = 0; j < npatches; j++ )
          {    const vec< pair<int,int> >& w = walks1[j];
               vec< pair<int,int> > wx;
               if ( w.size( ) == 2 && w[0].first == w[1].first )
               {    for ( int m = s1; m <= s2; m++ )
                         wx.push( u1, m );    }
               else
               {    for ( int l = 0; l < w.isize( ); l++ )
                    {    if ( l == 0 ) wx.push_back( w[l] );
                         else
                         {    for ( int m = w[l-1].second + 1;
                                   m <= unibases[ w[l-1].first ].isize( ) - K; m++ )
                              {    wx.push( w[l-1].first, m );    }
                              for ( int m = 0; m <= w[l].second; m++ )
                                   wx.push( w[l].first, m );    }    }    }
               int u1 = wx[0].first, s1 = wx[0].second;
               patches[j].resize( wx.isize( ) - 1 );
               for ( int l = 1; l < wx.isize( ); l++ )
               {    int u = wx[l].first, s = wx[l].second;
                    patches[j].Set( l-1, unibases[u][s+K-1] );    }    }

          // Filter using jump alignments.

          if ( jbases_sorted.nonempty( ) && patches.nonempty( ) )
          {    int L = jbases_sorted[0].size( );
               vec<int> hits( npatches, 0 ), id( npatches, vec<int>::IDENTITY );
               vec< vec<int> > HITS(npatches);
               for ( int i = 0; i < npatches; i++ )
               {    for ( int j = 0; j <= patches[i].isize( ) - L; j++ )
                    {    basevector b( patches[i], j, L );
                         for ( int pass = 1; pass <= 2; pass++ )
                         {    if ( pass == 2 ) b.ReverseComplement( );
                              int64_t low = LowerBound( jbases_sorted, b );
                              int64_t high = UpperBound( jbases_sorted, b );
                              for ( int64_t x = low; x < high; x++ )
                              {    int64_t id = jbases_sorted_id[x];
                                   int64_t idp = jpairs.getPartnerID(id);
                                   int xp = BinPosition1( jaligns, idp );
                                   if ( xp < 0 ) continue;
                                   int mx = jaligns[xp].second;
                                   int p = jaligns[xp].third;
                                   const int max_dist = 10000; // ******************
                                   if ( mx == m1 && tigs[m1].isize( ) - p <= max_dist
                                        && pass == 2 )
                                   {    // hits[i]++;
                                        HITS[i].push_back(id);    }
                                   if ( mx == -m2-1 && p <= max_dist && pass == 1 )
                                   {    // hits[i]++;    
                                        HITS[i].push_back(id);    }    }    }    }
                    if ( verbosity >= 3 ) PRINT2_TO( out, i, hits[i] );    }
               for ( int i = 0; i < npatches; i++ )
                    UniqueSort( HITS[i] );
               vec<Bool> to_delete( npatches, False );
               for ( int i1 = 0; i1 < npatches; i1++ )
               {    if ( to_delete[i1] ) continue;
                    for ( int i2 = 0; i2 < npatches; i2++ )
                    {    if ( !( HITS[i1].size( ) > HITS[i2].size( ) ) ) continue;
                         int better1 = 0, better2 = 0;
                         for ( int j = 0; j < HITS[i1].isize( ); j++ )
                              if ( !BinMember( HITS[i2], HITS[i1][j] ) ) better1++;
                         for ( int j = 0; j < HITS[i2].isize( ); j++ )
                              if ( !BinMember( HITS[i1], HITS[i2][j] ) ) better2++;
                         const int min_mult = 4;
                         if ( better1 >= min_mult * better2 )
                              to_delete[i2] = True;    }    }
               EraseIf( patches, to_delete );
               npatches = patches.size( );
               if ( verbosity >= 3 )
               {    out << "\npatches:\n";
                    for ( int j = 0; j < npatches; j++ )
                         patches[j].Print( out, "patch_" + ToString(j) );    }    }

          // Find the number of bases shared by the patches have on the left.  Then,
          // excluding the shared bases that have already been declared on the left,
          // determine the number of shared bases that the patches have on the
          // right.

          int left_share = 0, right_share = 0;
          GetShares( patches, left_share, right_share );
          basevector Left_share( patches[0], 0, left_share );
          basevector Right_share( patches[0], 
               patches[0].isize( ) - right_share, right_share );

          // Find the variant part "patches_diff" of the patches.  After this, the
          // joined contigs may be represented as
          // tigse[m1] 
          // + Left_share + {patches_diff} + Right_share 
          // + tigse[m2] (after deleting over2+K bases on left).
          
          vec<basevector> patches_diff(npatches);
          for ( int i = 0; i < npatches; i++ )
          {    patches_diff[i].SetToSubOf( patches[i], left_share, 
                    patches[i].isize( ) - left_share - right_share );    }
          vec<int> varsize(npatches);
          for ( int i = 0; i < npatches; i++ )
               varsize[i] = patches_diff[i].size( );

          // Build and filter efasta patches.  The current filter is very stringent.

          const int max_var_sum = 50;
          efasta epatch = Left_share.ToString( ) + ( npatches > 1 ? "{" : "" );
          for ( int i = 0; i < npatches; i++ )
               epatch += ( i > 0 ? "," : "" ) + patches_diff[i].ToString( );
          epatch += ( npatches > 1 ? "}" : "" ) + Right_share.ToString( );
          if ( Sum(varsize) <= max_var_sum )
          {    patched[s][it] = True;
               results[s][it] = make_triple( epatch, over1, over2 + K );    }

          // Generate report.

          if ( verbosity >= 1 )
          {    Sort(varsize);
               String varsizes = "{";
               for ( int i = 0; i < npatches; i++ )
               {    if ( i > 0 ) varsizes += ",";
                    varsizes += ToString( varsize[i] );    
                    int j;
                    for ( j = i + 1; j < npatches; j++ )
                         if ( varsize[j] != varsize[i] ) break;
                    if ( j - i > 1 ) varsizes += "^" + ToString(j-i);
                    i = j - 1;    }
               varsizes += "}";
               String share = "{" + ToString(left_share) + ","
                    + ToString(right_share) + "}";
               out << m1 << " --> " << m2 << ": ";
               String outcome = ( patched[s][it] ? "keep" : "trash" );
               PRINT5_TO( out, over1, over2, share, varsizes, outcome );
               report[gi] = out.str( );    }    }

     // Print reports.

     if ( verbosity >= 1 )
     {    cout << String( 85, '=' ) << "\n";
          for ( int gi = 0; gi < gaps.isize( ); gi++ )
               if ( report[gi].size( ) > 0 ) cout << report[gi];
          cout << String( 85, '=' ) << "\n";    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, 
          "linear_scaffolds0.clean.patched.fixed");
     CommandArgument_String_OrDefault(JUMPS_IN, "jump_reads_ec");
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, SCAFFOLDS_IN + ".kpatch");
     CommandArgument_Int(K);
     CommandArgument_String_OrDefault(HEAD, "all_reads");
     CommandArgument_Int_OrDefault(VERBOSITY, 0);
     CommandArgument_Bool_OrDefault(WRITE_EDITS, True);
     CommandArgument_Bool_OrDefault(WRITE_ASSEMBLY, False);
     CommandArgument_Int_OrDefault(MAX_ITERATIONS, 10000);
     CommandArgument_Int_OrDefault(MAX_GENOME_TO_FILTER_BY_JUMPS, 10 * 1000 * 1000);
     EndCommandArguments;

     // Define directories, etc.

     String data_dir = PRE + "/" + DATA; 
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     // Define input files.

     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     String efasta_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.efasta";
     String unibases_file = run_dir + "/" + HEAD + ".unibases.k" + ToString(K);

     // Load assembly and unibases.

     cout << Date( ) << ": loading assembly and unibases" << endl;
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );
     vec<efasta> tigse;
     LoadEfastaIntoStrings( efasta_file, tigse );
     vecbasevector tigs;
     tigs.ReadAll( sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fastb" );
     vecbasevector unibases(unibases_file);
     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );

     // Load jumps.

     vec<basevector> jbases_sorted;
     vec<int64_t> jbases_sorted_id;
     PairsManager jpairs;
     vec< triple<int64_t,int,int> > jaligns;
     int64_t assembly_size = 0;
     for ( size_t j = 0; j < tigs.size( ); j++ )
          assembly_size += tigs[j].size( );
     if ( assembly_size <= MAX_GENOME_TO_FILTER_BY_JUMPS )
     {    cout << Date( ) << ": loading jump data" << endl;
          vecbasevector jbases( run_dir + "/" + JUMPS_IN + ".fastb" );
          jpairs.Read( run_dir + "/" + JUMPS_IN + ".pairs" );
          jpairs.makeCache( );
          int min_read = 1000000000;
          for ( size_t id = 0; id < jbases.size( ); id++ )
               min_read = Min( min_read, jbases[id].isize( ) );
          DPRINT(min_read);
          cout << Date( ) << ": sorting jump data" << endl;
          jbases_sorted.resize( jbases.size( ) );
          jbases_sorted_id.resize( jbases.size( ) );
          for ( size_t id = 0; id < jbases.size( ); id++ )
          {    jbases_sorted[id] = jbases[id];
               if ( jbases[id].isize( ) > min_read )
               {    jbases_sorted[id].SetToSubOf( 
                         jbases[id], jbases[id].isize( ) - min_read, min_read );    }
               jbases_sorted_id[id] = id;    }
          ParallelSortSync( jbases_sorted, jbases_sorted_id );
          cout << Date( ) << ": looking up contigs" << endl;
          vec< triple<int64_t,int,int> > jaligns0;
          #pragma omp parallel for
          for ( size_t m = 0; m < tigs.size( ); m++ )
          {    vec< triple<int64_t,int,int> > jaligns0m;
               for ( int p = 0; p <= tigs[m].isize( ) - min_read; p++ )
               {    basevector b( tigs[m], p, min_read );
                    for ( int pass = 1; pass <= 2; pass++ )
                    {    if ( pass == 2 ) b.ReverseComplement( );
                         int64_t low = LowerBound( jbases_sorted, b );
                         int64_t high = UpperBound( jbases_sorted, b );
                         for ( int64_t l = low; l < high; l++ )
                         {    int64_t id = jbases_sorted_id[l];
                              int mx = m;
                              if ( pass == 2 ) mx = -m-1;
                              jaligns0m.push( id, mx, p );    }    }    }
               #pragma omp critical
               {    jaligns0.append(jaligns0m);    }    }
          cout << Date( ) << ": sorting aligns" << endl;
          ParallelSort(jaligns0);
          for ( size_t j1 = 0; j1 < jaligns0.size( ); j1++ )
          {    size_t j2;
               for ( j2 = j1 + 1; j2 < jaligns0.size( ); j2++ )
                    if ( jaligns0[j2].first != jaligns0[j1].first ) break;
               if ( j2 - j1 == 1 ) jaligns.push_back( jaligns0[j1] );
               j1 = j2 - 1;    }    }

     // Generate patches.

     vec< vec<Bool> > patched;
     vec< vec< triple<efasta,int,int> > > results;
     if ( K == 96 ) 
     {    PatchMe<96>( jbases_sorted, jbases_sorted_id, jpairs, jaligns, scaffolds, 
               tigse, tigs, unibases, nexts, patched, results, MAX_ITERATIONS, 
               VERBOSITY );    }
     else if ( K == 400 ) 
     {    PatchMe<400>( jbases_sorted, jbases_sorted_id, jpairs, jaligns, scaffolds,
               tigse, tigs, unibases, nexts, patched, results, MAX_ITERATIONS, 
               VERBOSITY );    }
     else if ( K == 640 ) 
     {    PatchMe<640>( jbases_sorted, jbases_sorted_id, jpairs, jaligns, scaffolds,
               tigse, tigs, unibases, nexts, patched, results, MAX_ITERATIONS, 
               VERBOSITY );    }
     else
     {    cout << "K value unsupported" << endl;
          exit(1);    }

     // Save patches.

     if (WRITE_EDITS)
     {    cout << Date( ) << ": saving patches" << endl;
          vec<assembly_edit> edits;
          for ( int s = 0; s < scaffolds.isize( ); s++ ) 
          {    const superb& S = scaffolds[s];
               for ( int p = S.Ngaps( ) - 1; p >= 0; p-- )
               {    if ( !patched[s][p] ) continue;
                    int m1 = S.Tig(p), m2 = S.Tig(p+1);
                    efasta M1 = tigse[m1], M2 = tigse[m2];
                    int trim1 = results[s][p].second, trim2 = results[s][p].third;
                    if ( trim1 > M1.isize( ) || trim2 > M2.isize( ) ) continue;
                    int start1 = M1.isize( ) - trim1;
                    int stop2 = trim2;
                    vec<basevector> reps;
                    results[s][p].first.ExpandTo(reps);
                    edits.push( assembly_edit::GAP_CLOSER, m1, 
                         start1, m2, stop2, reps );    }    }
          String outputFile = sub_dir + "/" + SCAFFOLDS_OUT + ".edits";
          BinaryWriter::writeFile( outputFile.c_str( ), edits );    }
     if (WRITE_ASSEMBLY)
     {    cout << Date( ) << ": installing patches" << endl;
          for ( int s = 0; s < scaffolds.isize( ); s++ )
          {    superb& S = scaffolds[s];
               for ( int p = S.Ngaps( ) - 1; p >= 0; p-- )
               {    if ( !patched[s][p] ) continue;
                    int m1 = S.Tig(p), m2 = S.Tig(p+1);
                    efasta M1 = tigse[m1], M2 = tigse[m2];
                    int trim1 = results[s][p].second, trim2 = results[s][p].third;
                    if ( trim1 > M1.isize( ) || trim2 > M2.isize( ) ) continue;
                    M1.resize( M1.isize( ) - trim1 );
                    M2 = M2.substr( trim2, M2.isize( ) - trim2 );
                    tigse[m1] = M1 + results[s][p].first + M2;
                    S.SetLen( p, tigse[m1].Length1( ) );
                    tigse[m2].resize(0);
                    int gap = 0, dev = 0;
                    bool nextgap = ( p+1 < S.Ngaps( ) );
                    if (nextgap) 
                    {    gap = S.Gap(p+1), dev = S.Dev(p+1);    }
                    S.RemoveTigByPos(p+1);
                    if (nextgap) 
                    {    S.SetGap(p, gap), S.SetDev(p, dev);    }    }    }
          Assembly A( scaffolds, tigse );
          A.remove_unused_contigs( );
          A.WriteAll( sub_dir + "/" + SCAFFOLDS_OUT );    }
     cout << Date( ) << ": done" << endl;    }
