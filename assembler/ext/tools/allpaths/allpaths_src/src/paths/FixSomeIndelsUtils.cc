///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "Basevector.h"
#include "PrintAlignment.h"
#include "STLExtensions.h"
#include "String.h"
#include "feudal/QualNibbleVec.h"
#include "graph/Digraph.h"
#include "paths/FixSomeIndelsUtils.h"
#include "paths/ReadLoc.h"
#include "polymorphism/DumbCall.h"

void SelectReads( const Bool VERBOSE,
		  const Bool USE_QUALS,
		  const int rep_start,
		  const int rep_stop,
		  const vec< triple<int64_t,int,int> > &RALIGNS_this,
		  const vec< triple<int64_t,int,int> > &JRALIGNS_this,
		  const String &run_dir,
		  const String &TIGS,			  
		  const BaseVecVec &bases,
		  const BaseVecVec &jbases,
		  const VecQualNibbleVec &quals,
		  const VecQualNibbleVec &jquals,
		  vecbasevector &reads,
		  VecQualNibbleVec &readsq,
		  ostream *p_rout )
{
     const int readlen = 101;   // NEED TO GET RID OF THIS!
     const int min_read_length = 50;
			    
     if (VERBOSE) 
          *p_rout << Date() << ": selecting fragment reads" << endl;
     int count = 0;
     for ( size_t i = 0; i < RALIGNS_this.size(); i++ )
     {    int start = RALIGNS_this[i].third;
          Bool rc = (start < 0);
          if (rc) start = -start-1;
          int stop = start + readlen;
          if (start > rep_start || stop < rep_stop) continue;
          int64_t rid = RALIGNS_this[i].first;
          BaseVecVec aread;
          if (TIGS == "") aread.push_back(bases[rid]);
          else aread.ReadOne(run_dir + "/frag_reads_filt_cpd.fastb", rid);
          if ( aread[0].isize( ) < min_read_length ) continue;
          if (rc) aread[0].ReverseComplement();
          reads.push_back_reserve(aread[0]);
          if ( quals.size( ) > 0 )
          {    if ( TIGS == "" ) readsq.push_back_reserve( quals[rid] );
               else
               {    vecqualvector areadq;
                    areadq.ReadOne( run_dir 
                         + "/frag_reads_filt_cpd.qualb", rid );
                    readsq.push_back_reserve( 
                         QualNibbleVec( areadq[0] ) );    }
               if (rc) readsq.back( ).ReverseMe( );    }
          count++;
          if (VERBOSE) 
          {    aread[0].Print(*p_rout, ToString(count) 
                    + "/" + ToString(rid));    }    }
     if (VERBOSE) PRINT_TO(*p_rout, count);
     if (VERBOSE) *p_rout << Date() << ": selecting jump reads" << endl;
     for (size_t i = 0; i < JRALIGNS_this.size(); i++) 
     {    int start = JRALIGNS_this[i].third;
          Bool rc = (start < 0);
          if (rc) start = -start-1;
          int stop = start + readlen;
          if ( start > rep_start || stop < rep_stop ) continue;
          int64_t rid = JRALIGNS_this[i].first;
          BaseVecVec aread;
          if (TIGS == "") aread.push_back(jbases[rid]);
          else aread.ReadOne(run_dir + "/jump_reads_filt_cpd.fastb", rid);
          if ( aread[0].isize( ) < min_read_length ) continue;
          if (rc) aread[0].ReverseComplement();
          aread[0].ReverseComplement(); // UGH!  BACKWARDS!!!
          reads.push_back_reserve(aread[0]);
          if ( USE_QUALS )
          {    if ( TIGS == "" ) readsq.push_back_reserve( jquals[rid] );
               else
               {    vecqualvector areadq;
                    areadq.ReadOne( run_dir 
                         + "/jump_reads_filt_cpd.qualb", rid );
                    readsq.push_back_reserve( 
                         QualNibbleVec( areadq[0] ) );    }
               if (rc) readsq.back( ).ReverseMe( );
               readsq.back( ).ReverseMe( );    } // UGH!  BACKWARDS!!!
          if (VERBOSE) 
          {    aread[0].Print(*p_rout, ToString(count++) + "/"
                    + ToString(rid));    }    }
     if (VERBOSE) PRINT_TO(*p_rout, count);     }

void FindBestHomeForReads( const Bool VERBOSE,
			   const String& title,
			   const vecbasevector& reads, 
			   const VecQualNibbleVec& quals,
			   const FastaVec& c,
			   const int nc1,
			   const int nc2, 
			   const double min_core_identity,
			   vec< vec<int> >& ERRS,
			   ostream *p_rout )
{    
     int max_rerrs = nc2 - int(ceil( double(nc2) * min_core_identity ));
     for ( size_t i = 0; i < reads.size( ); i++ )
     {    vec<int> starts;
          const basevector& R = reads[i];
          const QualNibbleVec& Q = quals[i];
          for ( int j = 0; j <= R.isize( ) - nc2; j++ )
          {    int e = 0;
               for ( int l = 0; l < nc2; l++ )
               {    if ( nc1+l < (int) c.size( ) && as_base(R[j+l]) != c[nc1+l] )
                    {    if ( ++e > max_rerrs ) break;    }    }
               if ( e <= max_rerrs ) starts.push_back( nc1 - j );    }
          int best_e = 1000000000, best_s = -1, best_q = 1000000000;
          for ( int si = 0; si < starts.isize( ); si++ ) 
          {    int s = starts[si], e = 0, q = 0;
               if ( s < 0 || s + R.isize( ) > (signed) c.size( ) ) continue;
               for ( int j = 0; j < R.isize( ); j++ )
               {    if ( as_base(R[j]) != c[s+j] )
                    {    q += Q[j];
                         if ( ++e >= best_e ) break;    }    } 
               if ( e < best_e )
               {    best_e = e; best_s = s; best_q = q;    }    }
          if ( best_s < 0 ) ERRS[i].push_back(-1000);
          else
          {    int q = best_q;
               int score = -1000;
               if ( q <= 100 ) score = -q;
               ERRS[i].push_back(score);
               if (VERBOSE) 
               {    avector<int> gaps(1), lengths(1);
                    gaps(0) = 0;
                    lengths(0) = R.size( );
                    alignment a;
                    a.Set( 0, best_s, 0, gaps, lengths );
                    *p_rout << "i = " << ", " << title << ", ";
                    PRINT2_TO( *p_rout, score, q );
                    FastaVec r(R);
                    PrintVisualAlignment( 
                         False, *p_rout, r, c, a );    }    }    }    }

void ShowFlakyRegions( const int n,
		       const int tig,
		       const double min_ref_frac,
		       const bvec &TIG,
		       const vecbvec &frag_reads,
		       const vec<dumbcall> &calls,
		       const vec<read_loc> &locs,
		       ostream &rout )
{    vec<int> flakes;
     for ( int i = 0; i < n; i++ )
     {    double ref_frac = double( calls[i].Count( TIG[i] ) ) 
               / double( calls[i].CountAll( ) );
          if ( ref_frac < min_ref_frac ) flakes.push_back(i);    }
     const int group_sep = 50;
     const int fflank = 10;
     for ( int j = 0; j < flakes.isize( ); j++ )
     {    int j2;
          for ( j2 = j + 1; j2 < flakes.isize( ); j2++ )
               if ( flakes[j2] - flakes[j2-1] > group_sep ) break;
          int left = flakes[j] - 2*fflank, right = flakes[j2-1] + 2*fflank;
          if ( left < 0 || right >= n )
          {    j = j2 - 1;
               continue;    }
          rout << "\nFLAKY REGION " << tig << "." << left << "-"
               << right << "\n";
          basevector fr;
          fr.SetToSubOf( TIG, left, right - left );
          rout << fr.ToString( ) << "\n";

           // Generate pkmer graph.

           vec< triple<basevector,int,int> > pkmers;
          vec< pair<int,int> > edges;
          const int pk = 12;
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 ) 
               {    Sort(pkmers);
                    vec< triple<basevector,int,int> > pkmers2;
                    for ( int l = 0; l < pkmers.isize( ); l++ )
                    {    int m = pkmers.NextDiff(l);
                         if ( m - l >= 2 ) pkmers2.push_back( pkmers[l] );
                         l = m - 1;    }
                    for ( int l = left; l <= right - pk; l++ )
                    {    basevector b;
                         b.SetToSubOf( TIG, l, pk );
                         pkmers2.push( b, l, l + pk );    }
                    UniqueSort(pkmers2);
                    pkmers = pkmers2;    }
               for ( int l = 0; l < locs.isize( ); l++ )
               {    const read_loc& rl = locs[l];
                    basevector b = frag_reads[ rl.ReadId( ) ];
                    if ( rl.Rc( ) ) b.ReverseComplement( );
                    if ( rl.Stop( ) + rl.Bandwidth( ) <= left ) continue;
                    if ( rl.Start( ) - rl.Bandwidth( ) >= right ) continue;
                    if ( pass == 1 )
                    {    rout << "\n";
                         rl.PrintVisualLoc( False, rout, "frag read " 
                              + ToString( rl.ReadId( ) ), b, TIG );    }
                    align a;
                    rl.GetAlign( a, b, TIG );
     
                          // Generate pkmers.
     
                          vec<int> to_p2( b.size( ) + 1 );
                    int p1 = a.pos1( ), p2 = a.pos2( );
                    to_p2[p1] = p2;
                    for ( int l = 0; l < a.Nblocks( ); l++ ) 
                    {    if ( a.Gaps(l) > 0 ) p2 += a.Gaps(l);
                         if ( a.Gaps(l) < 0 ) 
                         {    for ( int m = 0; m < -a.Gaps(l); m++ )
                                   to_p2[++p1] = p2;    }
                         for ( int x = 0; x < a.Lengths(l); x++ ) 
                              to_p2[++p1] = ++p2;    }
                    if ( pass == 1 )
                    {    for (int l = a.pos1( ); l <= a.Pos1( ) - pk; l++)
                         {    basevector x;
                              x.SetToSubOf( b, l, pk );
                              pkmers.push( 
                                   x, to_p2[l], to_p2[l+pk] );    }    }
                    else
                    {    vec<int> r( b.isize( ) + 1, -1 );
                         for (int l = a.pos1( ); l <= a.Pos1( ) - pk; l++)
                         {    basevector x;
                              x.SetToSubOf( b, l, pk );
                              triple<basevector,int,int>
                                   X( x, to_p2[l], to_p2[l+pk] );
                              r[l] = BinPosition( pkmers, X );    }
                         for (int l = a.pos1( ); l < a.Pos1( ) - pk; l++)
                         {    if ( r[l] >= 0 && r[l+1] >= 0 )
                                   edges.push( r[l], r[l+1] );    }
                              }    }    }
          int start = -1, stop = -1;
          for ( int l = left; l < right - pk; l++ )
          {    basevector b1, b2;
               b1.SetToSubOf( TIG, l, pk ), b2.SetToSubOf( TIG, l+1, pk );
               triple<basevector,int,int> X1( b1, l, l + pk );
               triple<basevector,int,int> X2( b2, l + 1, l + 1 + pk );
               int r1 = BinPosition( pkmers, X1 );
               if ( l == left ) start = r1;
               int r2 = BinPosition( pkmers, X2 );
               if ( l == right - pk - 1 ) stop = r2;
               edges.push( r1, r2 );    }
          UniqueSort(edges);
          int N = pkmers.size( );
          vec< vec<int> > from(N), to(N);
          for ( int l = 0; l < edges.isize( ); l++ )
          {    int v = edges[l].first, w = edges[l].second;
               from[v].push_back(w), to[w].push_back(v);    }
          for ( int l = 0; l < N; l++ )
          {    Sort( from[l] ), Sort( to[l] );    }
          digraph G( from, to );
          vec< vec<int> > paths;
          int maxpaths = 100;
          Bool succeed = G.AllPaths( start, stop, paths, maxpaths );
          rout << "\nPATHS:\n";
          if (succeed)
          {    vec<basevector> bpaths;
               for ( int l = 0; l < paths.isize( ); l++ )
               {    const vec<int>& p = paths[l];
                    basevector b = pkmers[ p[0] ].first;
                    for ( int j = 1; j < p.isize( ); j++ )
                         b.push_back( pkmers[ p[j] ].first.back( ) );
                    bpaths.push_back(b);
                    String s = b.ToString( );
                    PRINT2_TO( rout, l, s );    }    }

           // Report calls.

           rout << "\n";
          for ( int u = left; u < right; u++ )
          {    int ncalls = 0, ncalls_agree = 0;
               for ( int k = 0; k < 6; k++ )
               {    ncalls += calls[u].base[k];
                    if ( k == TIG[u] ) 
                         ncalls_agree += calls[u].base[k];    }
               double agree = double(ncalls_agree)/double(ncalls);
               rout << u << " ";
               for ( int k = 0; k < 6; k++ )
               {    for ( int l = 0; l < calls[u].base[k]; l++ )
                    {    if ( k < 4 ) rout << as_base(k);
                         if ( k == 4 ) rout << "D";
                         if ( k == 5 ) rout << "I";    }    }
               if ( ncalls > 0 )
               {    rout << " (" << as_base(TIG[u]) << " " 
                         << PERCENT_RATIO(3, ncalls_agree, ncalls)
                         << ")";    }
               rout << "\n";    }
          j = j2 - 1;    }    }





void SelectReads2( const Bool VERBOSE,
		  const int rep_start,
		  const int rep_stop,
		  const vec< triple<int64_t,int,int> > &RALIGNS_this,
		  const BaseVecVec &bases,
		  const VecQualNibbleVec &quals,
		  vecbasevector &reads,
		  VecQualNibbleVec &readsq,
		  ostream *p_rout )
{
     const int readlen = 101;   // NEED TO GET RID OF THIS!
     const int min_read_length = 50;
			    
     if (VERBOSE) *p_rout << Date() << ": selecting reads" << endl;
     int count = 0;
     for ( size_t i = 0; i < RALIGNS_this.size(); i++ )
     {    int start = RALIGNS_this[i].third;
          Bool rc = (start < 0);
          if (rc) start = -start-1;
          int stop = start + readlen;
          if (start > rep_start || stop < rep_stop) continue;
          int64_t rid = RALIGNS_this[i].first;
          BaseVec read = bases[rid];
          if ( read.isize( ) < min_read_length ) continue;
          if (rc) read.ReverseComplement();
          reads.push_back_reserve(read);
          if ( quals.size( ) > 0 )
          {    readsq.push_back_reserve( quals[rid] );
	       if (rc) readsq.back( ).ReverseMe( );   
	  }
	  count++;
          if (VERBOSE) 
          {    read.Print(*p_rout, ToString(count) 
                    + "/" + ToString(rid));    }    }
     if (VERBOSE) PRINT_TO(*p_rout, count);     }


void FindBestHomeForReads2( const Bool VERBOSE,
			   const String& title,
			   const vecbasevector& reads, 
			   const VecQualNibbleVec& quals,
			   const FastaVec& c,
			   const int nc1,
			   const int nc2, 
			   const double min_core_identity,
			   vec< vec<int> >& ERRS,
			   ostream *p_rout )
{ 
  const int MAX_SCORE = 100000;  
  int max_rerrs = nc2 - int(ceil( double(nc2) * min_core_identity ));
  for ( size_t i = 0; i < reads.size( ); i++ )
  {    
    vec<int> starts;
    const basevector& R = reads[i];
    const QualNibbleVec& Q = quals[i];
    for ( int j = 0; j <= R.isize( ) - nc2; j++ ) {    
      int e = 0;
      for ( int l = 0; l < nc2; l++ ) {    
	if ( nc1+l < (int) c.size( ) && as_base(R[j+l]) != c[nc1+l] )
	  if ( ++e > max_rerrs ) break;    
      }
      if ( e <= max_rerrs ) starts.push_back( nc1 - j );    
    }
    //*p_rout << "starts.size()=  " << starts.size() << endl;
    int best_e = 1000000000, best_s = -1000000000, best_q = 1000000000;
    for ( int si = 0; si < starts.isize( ); si++ ) {    
      int s = starts[si], e = 0, q = 0;
      //if ( s < 0 || s + R.isize( ) > (signed) c.size( ) ) continue;
      for ( int j = 0; j < R.isize( ); j++ ) {    
	if ( s + j < 0 || s + j >= (signed) c.size() ) continue;
	if ( as_base(R[j]) != c[s+j] ) {    
	  q += Q[j];
	  if ( ++e >= best_e ) break;    
	}    
      } 
      if ( q < best_q ) { best_e = e; best_s = s; best_q = q; }    
    }
    if ( best_s == -1000000000 ) ERRS[i].push_back( -MAX_SCORE );
    else {    
      int q = best_q;
      int score = - MAX_SCORE;
      if ( q <= MAX_SCORE) score = -q;
      ERRS[i].push_back(score);
      if (VERBOSE) {    
	*p_rout << "i = " << i << ", " << title << ", ";
	PRINT2_TO( *p_rout, score, q );
	*p_rout << "  best_e,s,q "<< best_e << " " << best_s << " " << best_q << endl;
	fastavector cDisplay;
	if ( best_s > 0 ) {
	  cDisplay.SetToSubOf( c, best_s, Min( c.size() - best_s, R.size() ) );
	  *p_rout << R.ToString() << endl;
	  *p_rout << cDisplay.ToString() << endl;
	} else {
	  *p_rout << R.ToString() << endl;
	  *p_rout << String(-best_s, ' ');
	  cDisplay.SetToSubOf( c, 0, Min( c.size(), R.size() - best_s ) );
	  *p_rout << cDisplay.ToString() << endl;
	}
	avector<int> gaps(1), lengths(1);
	gaps(0) = 0;
	lengths(0) = R.size( );
	//alignment a;
	//a.Set( 0, best_s, 0, gaps, lengths );
	//FastaVec r(R);
	//PrintVisualAlignment(False, *p_rout, r, c, a );    
      }    
    }    
  }    
}
