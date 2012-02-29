// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "pairwise_aligners/AlignTwoBasevectors.h"
#include "AnAssemblyClass.h"
#include "AnAssemblyClassUtil.h"
#include "AnnotatedContig.h"
#include "Basevector.h"
#include "ContigEvent.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "math/Functions.h"
#include "FindGaps.h"
#include "Quality.h"
#include "ReadLocation.h"
#include "ScoreAlignment.h"
#include "Set.h"
#include "ShortVector.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "STLExtensions.h"
#include "SuperMap.h"
#include "VecInserts.h"

void LoadAllAligns( const String &aligns_file,
		    int n_reads,
		    vec<brief_align> &all_aligns,
		    vec<int> &all_aligns_index )
{
  READX( aligns_file, all_aligns );

  all_aligns_index.resize( n_reads, -1 );
  for (int ii=(int)all_aligns.size( )-1; ii>=0; ii--)
    all_aligns_index[ all_aligns[ii].Id1() ] = ii;
}

void LoadVecSupers( const String &ann_supers_file, vec<super> &supers )
{
  supers.clear();

  READ( ann_supers_file, vec<annotated_supercontig>, asupers );
  supers.resize( asupers.size( ) );
  for ( unsigned int i = 0; i < supers.size( ); i++ ) {
    annotated_supercontig& as = asupers[i];
    supers[i].mtig.resize( as.NumContigs() );
    supers[i].gap.resize( as.NumContigs() - 1 );
    for ( int j = 0; j < as.NumContigs( ); j++ ) {
      supers[i].mtig[j] = as.Contig(j).ID( );
      if ( j < as.NumContigs( ) - 1 )
	supers[i].gap[j] = as.Gap(j);
    }
  }

  Destroy( asupers );
}

void ReorderSuper( assembly& A, int s )
{
  super &sup1 = A.supers[s];

  map<int,int> position_to_contig_map;

  int start = 0;
  for ( unsigned int i = 0; i < sup1.mtig.size(); ++i )
  {
    if ( i > 0 )
    {
      start += A.Len( sup1.mtig[i-1] );
      start += sup1.gap[i-1];
    }
    position_to_contig_map[start] = sup1.mtig[i];
  }

  sup1.mtig.clear();
  sup1.gap.clear();

  map<int,int>::iterator new_order_iter = position_to_contig_map.begin();
  map<int,int>::iterator last_iter = new_order_iter;

  while ( new_order_iter != position_to_contig_map.end() )
  {
    if ( last_iter != new_order_iter )
    {
      sup1.gap.push_back( new_order_iter->first -
			  last_iter->first -
			  A.Len( last_iter->second ) );
    }

    A.mtigs_to_super_pos[ new_order_iter->second ] = sup1.mtig.size();
    sup1.mtig.push_back( new_order_iter->second );

    last_iter = new_order_iter++;
  }
}

void ReorderSuperMaxGap( assembly& A, int s )
{
  const vecbasevector &contigs = A.mtig;
  super &old_super = A.supers[s];
  super_map s_map( old_super, contigs );

  vec<int> sorted_contigs( s_map.GetNumberContigs( ) );
  for (int ii=0; ii<(int)sorted_contigs.size( ); ii++)
    sorted_contigs[ii] = ii;

  order_contig_pos_MaxGap sorter( &s_map );
  if ( !is_sorted( sorted_contigs.begin( ), sorted_contigs.end( ), sorter ) ) {
    sort( sorted_contigs.begin( ), sorted_contigs.end( ), sorter );

    super new_super;
    new_super.mtig.resize( old_super.mtig.size( ) );
    new_super.gap.resize( old_super.gap.size( ) );

    for (int jj=0; jj<(int)new_super.mtig.size( ); jj++)
    {
      new_super.mtig[jj] = s_map.GetContigId( sorted_contigs[jj] );
      A.mtigs_to_super_pos[ new_super.mtig[jj] ] = jj;
    }

    for (int jj=0; jj<(int)new_super.gap.size( ); jj++)
      new_super.gap[jj]
	= s_map.GetStartOnSupercontigPos( sorted_contigs[jj+1] )
	- ( s_map.GetStopOnSupercontigPos( sorted_contigs[jj] ) + 1 );

    old_super = new_super;
  }
}

void CheckSuperOverlaps( ostream &log, assembly &Ass, int scg )
{
  // Start and stop of contigs on supercontig.
  super the_super = Ass.supers[ scg ];

  int n_contigs = Ass.SuperSize( scg );
  vec<int> start( n_contigs );
  vec<int> stop( n_contigs );

  start[0] = 0;
  stop[0] = -1 + Ass.Len( the_super.mtig[0] );
  for (int cg=1; cg<n_contigs; cg++) {
    start[cg] = 1 + stop[cg-1] + the_super.gap[cg-1];
    stop[cg] = -1 + start[cg] + Ass.Len( the_super.mtig[cg] );
  }

  // Main loop over all contigs.
  align al;

  for (int cg=0; cg<n_contigs-1; cg++) {
    for (int ii=cg+1; ii<n_contigs; ii++) {
      // If contig ii is declared not to overlap contig cg.
      if ( start[ii] >= stop[cg] )
	break;

      int cg_id = the_super.mtig[ cg ];
      int ii_id = the_super.mtig[ ii ];

      int cg_length = Ass.Len( cg_id );
      int ii_length = Ass.Len( ii_id );

      int declared_overlap
	= ( stop[ii] <= stop[cg] ) ? ii_length : 1 + stop[cg] - start[ii];

      int observed_overlap = Overlap( Ass, cg_id, ii_id, al );

      // Log.
      String sz_cg
	= ToString( cg_id ) + " (length=" + ToString( cg_length ) + ")";
      String sz_ii
	= ToString( ii_id ) + " (length=" + ToString( ii_length ) + ")";
      String sz_align
	= ToString( declared_overlap ) + " ("
	+ ToString( observed_overlap ) + ")";

      log << "super_" << scg << "\t"
	  << sz_cg << "   "
	  << sz_ii << "   "
	  << sz_align << "\n";
    }
  }

}


Bool Connected( assembly& a, int m1, int m2, Bool strong )
{    vec<int>& m1reads = a.reads_orig_index[m1];
     for ( unsigned int i = 0; i < m1reads.size( ); i++ )
     {    read_location& r1 = a.reads_orig[ m1reads[i] ];
          int id1 = r1.ReadId( );
          if ( a.pairs_index[id1] >= 0 )
          {    read_pairing& p = a.pairs[ a.pairs_index[id1] ];
               int id2 = p.Partner(id1);
               int si = a.simple_reads_orig_index[id2];
               if ( si >= 0 )
               {    read_location& r2 = a.reads_orig[si];
                    if ( strong && r1.OrientationOnContig( ) != ForwardOr ) continue;
                    if ( strong && r2.OrientationOnContig( ) != ReverseOr ) continue;
                    if ( r2.Contig( ) == m2 &&
                         r1.OrientationOnContig( ) != r2.OrientationOnContig( ) )
                         return True;    }    }    }
     return False;    }


int ConnectionCount( assembly& a, int m1, int m2 )
{    int answer = 0;
     vec<int>& m1reads = a.reads_orig_index[m1];
     for ( unsigned int i = 0; i < m1reads.size( ); i++ )
     {    read_location& r1 = a.reads_orig[ m1reads[i] ];
          int id1 = r1.ReadId( );
          if ( a.pairs_index[id1] >= 0 )
          {    read_pairing& p = a.pairs[ a.pairs_index[id1] ];
               int id2 = p.Partner(id1);
               int si = a.simple_reads_orig_index[id2];
               if ( si >= 0 )
               {    read_location& r2 = a.reads_orig[si];
                    if ( r2.Contig( ) == m2
                         && r1.OrientationOnContig( ) == ForwardOr
                         && r2.OrientationOnContig( ) == ReverseOr )
                         ++answer;
	       }
	  }
     }
     return answer;
}

int IllogicalLinksCount( assembly &a, int m1, int m2 )
{
  int answer(0);
  vec<int>& m1reads = a.reads_orig_index[m1];

  for ( unsigned int i = 0; i < m1reads.size( ); i++ )
  {
    read_location& r1 = a.reads_orig[ m1reads[i] ];
    int id1 = r1.ReadId( );
    if ( a.pairs_index[id1] >= 0 )
    {
      read_pairing& p = a.pairs[ a.pairs_index[id1] ];
      int id2 = p.Partner(id1);
      int si = a.simple_reads_orig_index[id2];
      if ( si >= 0 )
      {
	read_location& r2 = a.reads_orig[si];
	if ( r2.Contig( ) == m2 &&
             ( ( r1.OrientationOnContig( ) == ForwardOr &&
                 r2.OrientationOnContig( ) == ForwardOr ) ||
               ( r1.OrientationOnContig( ) == ReverseOr &&
                 r2.OrientationOnContig( ) == ReverseOr ) ) )
	  ++answer;
      }
    }
  }
  return answer;
}


Bool ReadsOverlap( const assembly& a, int m1, int m2, Bool RC )
{    const vec<int>& ri = a.reads_orig_index[m1];
     for ( unsigned int i = 0; i < ri.size( ); i++ )
     {    const read_location& r1 = a.reads_orig[ ri[i] ];
          int id1 = r1.ReadId( );
          int ind = a.all_aligns_index[id1];
          for ( int j = ind; ind >= 0 && j < (int) a.all_aligns.size( ); j++ )
          {    const brief_align& ap = a.all_aligns[j];
               if ( ap.Id1( ) > id1 ) break;
               int id2 = ap.Id2( );
               int i2 = a.simple_reads_orig_index[id2];
               if ( i2 < 0 ) continue;
               const read_location& r2 = a.reads_orig[i2];
               if ( r2.Contig( ) != m2 ) continue;
               if ( ( RC + r1.OrientationOnContig( )
                    + r2.OrientationOnContig( ) + ap.Rc2( ) ) % 2 == 0 )
                    return True;    }    }
     return False;    }

int OverlapUpdate( vec< set<int> > &no_overlap,
     assembly& a, int m1, int m2, align& al, int kmer_size,
     int min_mutmer, Bool RC, int mode, int min_align, int min_align_see )
{
  // String m1_lost = a.mtig_aligns_lost[m1] ? " lost" : "";
  // String m2_lost = a.mtig_aligns_lost[m2] ? " lost" : "";
  // String is_RC = ( RC ) ? " RC" : "";

  if ( Member( no_overlap[m1], m2 ) && Member( no_overlap[m2], m1 ) ) return 0;

  int overlap = Overlap(  a, m1, m2, al, kmer_size, min_mutmer,
			  RC, mode, min_align, min_align_see );

  if ( overlap > 0 ) {
    // TODO - update mtig_aligns and mtig_aligns_index.
  }
  else
  {    no_overlap[m1].insert(m2);
       no_overlap[m2].insert(m1);    }

  return overlap;
}


int Overlap( const assembly& a, int m1, int m2, align& al, int kmer_size,
     int min_mutmer, Bool RC, int mode, int min_align, int min_align_see,
     Bool require_read_evidence )
{
     min_align_see = Min( min_align, min_align_see );

     // Overlap often gets called twice in a row with the same arguments.
     // What follows is a stupid way to handle this (and also a slightly
     // risky way).

     static int last_m1 = -1, last_m2 = -1;
     static int last_pos1 = -1, last_pos2 = -1;
     static Bool last_RC;
     static int return_value;
     if ( m1 == last_m1 && m2 == last_m2 && RC == last_RC
          && al.pos1( ) == last_pos1 && al.pos2( ) == last_pos2 )
          return return_value;
     last_m1 = m1;
     last_m2 = m2;
     last_RC = RC;

     const basevector &b1 = a.mtig[m1], &b2 = a.mtig[m2];
     const qualvector &q1 = a.mtig_qual[m1], &q2 = a.mtig_qual[m2];

     // If overlapping very large sequences, make some performance enhancements.

     longlong length_product = (longlong) b1.size( ) * (longlong) b2.size( );
     longlong max_product = (longlong) 500 * (longlong) 1000 * (longlong) 1000;
     Bool alt_method = ( length_product <= max_product );
     if ( length_product > max_product && kmer_size <= 12 ) kmer_size = 24;
     int stretch = ( length_product <= max_product ? 15 : 3 );

     Bool end_alignments_lost = a.mtig_aligns_lost[m1] || a.mtig_aligns_lost[m2];

     // Attempt to show that there is no alignment between m1 and m2, by looking
     // for good alignments between reads in m1 and reads in m2.
     //
     // Note that this could be expedited by sorting all_aligns first by id,
     // and then by score.
     //
     // If there is a good alignment between reads, estimate the offset between
     // m1 and m2 (if an alignment between them exists).

     Bool alignment_found = False;
     static vec<int> offsets;
     offsets.clear( );
     const vec<int>& ri = a.reads_orig_index[m1];
     for ( unsigned int i = 0; i < ri.size( ); i++ )
     {    const read_location& r1 = a.reads_orig[ ri[i] ];
          int id1 = r1.ReadId( );
          int ind = a.all_aligns_index[id1];
          for ( int j = ind; ind >= 0 && j < (int) a.all_aligns.size( ); j++ )
          {    const brief_align& ap = a.all_aligns[j];
               if ( ap.Id1( ) > id1 ) break;
               int id2 = ap.Id2( );
               int i2 = a.simple_reads_orig_index[id2];
               if ( i2 < 0 ) continue;
               const read_location& r2 = a.reads_orig[i2];
               if ( r2.Contig( ) != m2 ) continue;

               if ( ( RC + r1.OrientationOnContig( )
                    + r2.OrientationOnContig( ) + ap.Rc2( ) ) % 2 != 0 )
                    continue;
               int start2;
               if ( !RC ) start2 = r2.StartOnContig( );
               else start2 = r2.LengthOfContig( ) - r2.StopOnContig( ) - 1;
               int offset = ap.Offset( );
               if ( r1.OrientationOnContig( ) == ReverseOr )
                    offset = r1.LengthOfRead( ) - r2.LengthOfRead( ) - offset;
               int o = offset + r1.StartOnContig( ) - start2;

               if ( (o > r1.LengthOfContig( ) && o > r2.LengthOfContig( )) //WWW
                    || r1.LengthOfContig( ) != (int) b1.size( ) //WWW
                    || r2.LengthOfContig( ) != (int) b2.size( ) ) //WWW
               {    cout << "Warning: strange data encountered:" << endl; //WWW
                    PRINT2_TO(cout, m1, m2 ); //WWW
                    PRINT2_TO(cout, b1.size( ), b2.size( ) ); //WWW
                    PRINT2_TO(cout, r1.LengthOfContig( ), // WWW
                         r2.LengthOfContig( ) ); //WWW
                    PRINT2_TO(cout, r1.StartOnContig( ), r1.StopOnContig( ) ); //WWW
                    PRINT2_TO(cout, r2.StartOnContig( ), r2.StopOnContig( ) ); //WWW
                    PRINT2_TO(cout, r1.LengthOfRead( ), r2.LengthOfRead( ) ); //WWW
                    PRINT_TO(cout,  int( r1.OrientationOnContig( ) ) ); //WWW
                    PRINT_TO(cout,  int( r2.OrientationOnContig( ) ) ); //WWW
                    PRINT_TO(cout,  int(RC) ); //WWW
                    PRINT_TO(cout,  int( ap.Rc2( ) ) ); //WWW
                    // PRINT2_TO( cout,  ap.a.pos1( ), ap.a.Pos1( ) ); //WWW
                    // PRINT2_TO( cout,  ap.a.pos2( ), ap.a.Pos2( ) ); //WWW
                    PRINT_TO( cout, o); //WWW
                    continue;    } //WWW

               alignment_found = True;
               offsets.push_back(o);    }    }

     if ( !alignment_found && require_read_evidence ) return -1;

     // Eliminate spurious offsets.  More specifically, if an offset implies an
     // overlap of n kb, then we expect there to be at least n offsets, in total.
     // Otherwise, the offset is regarded as spurious.  Obviously, this condition
     // could be refined, and it would probably be worthwhile to do so.

     offset_check:
     for ( unsigned int i = 0; i < offsets.size( ); i++ )
     {    int off = offsets[i];
          int over;
          if ( off >= 0 ) over = Min( b1.size( ) - off, b2.size( ) );
          else over = Min( b2.size( ) + off, b1.size( ) );
          int required_offsets = over / 1000;
          if ( (int) offsets.size( ) < required_offsets )
          {    offsets.erase( offsets.begin( ) + i );
               goto offset_check;    }    }
     if ( offsets.size( ) == 0 ) alignment_found = False;

     int min_offset = 0, max_offset = 0;
     if (alignment_found)
     {    Sort(offsets);
          min_offset = offsets.front( );
          max_offset = offsets.back( );    }

     Bool end_overlap_known = False, end_overlap_ambiguous = False;
     int end_offset_value = 0;
     Bool end_rc2_value = False;
     if ( !end_alignments_lost )
     {    int ci = a.mtig_aligns_index[m1];
          if ( ci >= 0 )
          {    for ( unsigned int i = ci; i < a.mtig_aligns.size( ); i++ )
               {    if ( a.mtig_aligns[i].id1 != m1 ) break;
                    if ( a.mtig_aligns[i].id2 > m2 ) break;
                    if ( a.mtig_aligns[i].id2 == m2 && a.mtig_aligns[i].rc2 == RC )
                    {    if ( i < a.mtig_aligns.size( ) - 1
                              && a.mtig_aligns[i+1].id1 == a.mtig_aligns[i].id1
                              && a.mtig_aligns[i+1].id2 == a.mtig_aligns[i].id2
                              && a.mtig_aligns[i+1].rc2 == a.mtig_aligns[i].rc2 )
                         {    end_overlap_ambiguous = True;
                              break;    }
                         end_overlap_known = True;
                         end_offset_value = a.mtig_aligns[i].offset;
                         end_rc2_value = a.mtig_aligns[i].rc2;

                         if ( !alignment_found ||
                              ( Abs( end_offset_value - min_offset ) <= 15
                                   && Abs( end_offset_value - max_offset ) <= 15 ) )
                              goto smithwat;

                         break;    }    }    }    }

     if ( !end_alignments_lost && !end_overlap_known && !end_overlap_ambiguous
          && alignment_found && max_offset - min_offset <= 15 )
     {    end_offset_value = (min_offset + max_offset)/2;
          goto smithwat;    }

     if ( !end_overlap_known && offsets.size( ) >= 3
          && max_offset - min_offset <= 20 )
     {    end_offset_value = (min_offset + max_offset)/2;
          goto smithwat;    }

     if ( !end_alignments_lost && !alignment_found && !end_overlap_ambiguous
          && !end_overlap_known && b1.size( ) >= 500 && b2.size( ) >= 500 )
     {    return_value = -1;
          return -1;    }

     if ( alignment_found && !end_alignments_lost && !end_overlap_known
          && !end_overlap_ambiguous && kmer_size < 24 ) kmer_size = 24;

     if ( 0 == 1 )
     {    smithwat:
          if ( end_offset_value > (int) b1.size( ) // WWW
               || end_offset_value < - (int) b2.size( ) ) // WWW
          {    cout << "Warning: end_offset_value is illegal" << endl; // WWW
               PRINT2_TO( cout,  b1.size( ), b2.size( ) ); // WWW
               PRINT_TO( cout, end_offset_value); // WWW
               return_value = -1; // WWW
               return -1;    } // WWW
          if ( !RC )
          {    int errs;
               SmithWatBandedA( b1, b2, end_offset_value, 20, al, errs );
               float score = ScoreAlignment( al, b1, q1, b2, q2 );
               int errors = ActualErrors( b1, b2, al );
               int overlap = al.Pos1( ) - al.pos1( );
               // cout << "using offset = " << al.Offset( ) << "\n"; // QQQ
               // PRINT3_TO( cout,  score, errors, overlap ); // QQQ
               // PrintVisualAlignment( True, cout, b1, b2, al, q1, q2 ); // QQQ
               if ( (score <= 20000 || errors <= 15) && overlap >= 50 )
               {    last_pos1 = al.pos1( );
                    last_pos2 = al.pos2( );
                    return_value = overlap;
                    return overlap;    }
               else
               {    return_value = -1;
                    return -1;    }    }
          else
          {    static basevector b2r;
               b2r.Setsize( b2.size( ) );
               b2r.ReverseComplement(b2);
               static qualvector q2r;
               q2r = Reverse(q2);
               int errs;
               SmithWatBandedA( b1, b2r, end_offset_value, 20, al, errs );
               float score = ScoreAlignment( al, b1, q1, b2r, q2r );
               int errors = ActualErrors( b1, b2r, al );
               int overlap = al.Pos1( ) - al.pos1( );
               // cout << "using offset = " << al.Offset( ) << "\n"; // QQQ
               // PRINT3_TO( cout,  score, errors, overlap ); // QQQ
               // PrintVisualAlignment( True, cout, b1, b2r, al, q1, q2r ); // QQQ
               if ( (score <= 20000 || errors <= 15) && overlap >= 50 )
               {    last_pos1 = al.pos1( );
                    last_pos2 = al.pos2( );
                    return_value = overlap;
                    return overlap;    }
               else
               {    return_value = -1;
                    return -1;    }    }    }

     int rc = RC;

     // If we haven't found an alignment, check for an alignment between the
     // last 750 bases.  Exit if we can't find it.

     if ( !alignment_found && b1.size( ) >= 750 && b2.size( ) >= 750 )
     {    static basevector end1(750), end2(750);
          int thismode;
          if ( !RC )
          {    end1.SetToSubOf( b1, b1.size( ) - 750, 750 );
               end2.SetToSubOf( b2, 0, 750 );
               thismode = 8;    }
          else
          {    end1.SetToSubOf( b1, b1.size( ) - 750, 750 );
               end2.SetToSubOf( b2, b2.size( ) - 750, 750 );
               thismode = 9;    }
          int answer =
               AlignTwoBasevectors( end2, end1, al, min_align_see, 100000, 0.25,
                    0, rc, thismode, kmer_size, stretch, 1, q1, q2, 20000, 15, cout,
                    False, 100, 10000, 1000, 100, alt_method, 40, min_mutmer );
          if ( answer < min_align ) answer = 0;
          if ( answer > 0 && al.pos1( ) == 0 )
          {    return_value = answer;
               al.Flip( );
               al.AddToPos1( b1.size( ) - 750 );
               last_pos1 = al.pos1( );
               last_pos2 = al.pos2( );
               return answer;    }
          if ( answer <= 0 )
          {    if ( !RC )
               {    end1.SetToSubOf( b1, 0, 750 );
                    end2.SetToSubOf( b2, b2.size( ) - 750, 750 );    }
               else
               {    end1.SetToSubOf( b1, 0, 750 );
                    end2.SetToSubOf( b2, 0, 750 );    }
               answer = AlignTwoBasevectors( end2, end1, al, min_align_see, 100000,
                         0.25, 0, rc, thismode, kmer_size, stretch, 1, q1, q2,
                         20000, 15, cout, False, 100, 10000, 1000, 100,
                         alt_method, 40, min_mutmer );
               if ( answer < min_align ) answer = 0;
               if ( answer <= 0 )
               {    last_pos1 = al.pos1( );
                    last_pos2 = al.pos2( );
                    return_value = answer;
                    return answer;    }
               if ( answer > 0 && al.pos2( ) == 0 )
               {    return_value = answer;
                    al.Flip( );
                    al.AddToPos2( b2.size( ) - 750 );
                    last_pos1 = al.pos1( );
                    last_pos2 = al.pos2( );
                    return answer;    }    }    }

     if ( !RC )
     {    int answer =
               AlignTwoBasevectors( b1, b2, al, min_align_see, 100000, 0.25, 0,
               rc, mode, kmer_size, stretch, 1, q1, q2, 20000, 15, cout, False,
               100, 10000, 1000, 100, alt_method, 40, min_mutmer );
          if ( answer < min_align ) answer = 0;
          last_pos1 = al.pos1( );
          last_pos2 = al.pos2( );
          return_value = answer;
          return answer;    }
     static basevector b2r;
     b2r.Setsize( b2.size( ) );
     b2r.ReverseComplement(b2);
     static qualvector q2r;
     q2r = Reverse(q2);
     int answer =
          AlignTwoBasevectors( b1, b2r, al, min_align_see, 100000, 0.25, 0, rc,
          mode, kmer_size, stretch, 1, q1, q2r, 20000, 15, cout, False, 100, 10000,
          1000, 100, alt_method, 40, min_mutmer );
     if ( answer < min_align ) answer = 0;
     last_pos1 = al.pos1( );
     last_pos2 = al.pos2( );
     return_value = answer;
     return answer;    }

void ReverseMtig( assembly& a, int m )
{    a.mtig[m].ReverseComplement( );
     ReverseThis( a.mtig_qual[m] );
     for ( unsigned int j = 0; j < a.reads_index[m].size( ); j++ )
          a.reads[ a.reads_index[m][j] ].Reverse( );
     for ( unsigned int j = 0; j < a.reads_orig_index[m].size( ); j++ )
          a.reads_orig[ a.reads_orig_index[m][j] ].Reverse( );

     // Revise mtig_aligns entries.

     int ind = a.mtig_aligns_index[m];
     if ( ind >= 0 )
     {    for ( int j = ind; j < (int) a.mtig_aligns.size( ); j++ )
          {    brief_align& b1 = a.mtig_aligns[j];
               if ( b1.id1 != m ) break;
               if ( a.mtig_aligns_lost[ b1.id1 ] || a.mtig_aligns_lost[ b1.id2 ] )
                    continue;
               b1.offset = a.mtig[b1.id1].size( ) - a.mtig[b1.id2].size( )
                    - b1.offset;
               b1.rc2 = !b1.rc2;
               int ind2 = a.mtig_aligns_index[ b1.id2 ];
               for ( int k = ind2; k < (int) a.mtig_aligns.size( ); k++ )
               {    brief_align& b2 = a.mtig_aligns[k];
                    if ( b2.id1 != b1.id2 || b2.id2 > b1.id1 ) break;
                    if ( b2.id2 == b1.id1 ) b2.rc2 = !b2.rc2;    }    }    }    }

void MergeMtigsHead( const basevector& b1, const basevector& b2,
     const qualvector& q1, const qualvector& q2, const align& al,
     basevector& c, qualvector& q, vec<int>& to1, vec<int>& to2 )
{
     // Do sanity check.

     ForceAssertEq( b1.size( ), q1.size( ) );
     ForceAssertEq( b2.size( ), q2.size( ) );
     ForceAssertLt( al.pos2( ), (int) b2.size( ) );
     ForceAssertLt( al.pos1( ), (int) b1.size( ) );

     // Merge the contigs themselves.

     MergeTwoBaseVectors( b1, b2, al, c, q1, q2, q );

     // Set up map for translating positions to the merged contig.
     // This code mirrors MergeTwoBasevectors.

     to1.resize_and_set( b1.size( ), -1 ), to2.resize( b2.size( ), -1 );
     int pos1 = al.pos1( ), pos2 = al.pos2( );
     int Pos1 = al.Pos1( ), Pos2 = al.Pos2( );
     int nblocks = al.Nblocks( );
     const avector<int> &gaps = al.Gaps( ), &lengths = al.Lengths( );
     int i = 0, p1 = pos1, p2 = pos2;
     if ( pos2 == 0 )
     {    for ( ; i < pos1; i++ )
          {    to1[i] = i;    }    }
     else
     {    for ( ; i < pos2; i++ )
          {    to2[i] = i;    }    }
     for ( int j = 0; j < nblocks; j++ )
     {    if ( gaps(j) > 0 )
          {    for ( int x = 0; x < gaps(j); x++ )
               {    to2[p2++] = i++;    }    }
          if ( gaps(j) < 0 )
          {    for ( int x = 0; x < -gaps(j); x++ )
               {    to1[p1++] = i++;    }    }
          for ( int x = 0; x < lengths(j); x++ )
          {    to1[p1++] = i, to2[p2++] = i;
               ++i;    }    }
     for ( int j = 0; j < (int) b1.size( ) - Pos1; j++ )
     {    if ( j + Pos1 >= (int) b1.size( ) ) break;
          to1[j + Pos1] = i++;
          if ( to1[j + Pos1] >= (int) c.size( ) )
               to1[j + Pos1] = (int) c.size( ) - 1;    }
     for ( int j = 0; j < (int) b2.size( ) - Pos2; j++ )
     {    if ( j + Pos2 >= (int) b2.size( ) ) break;
          to2[j + Pos2] = i++;
          if ( to2[j + Pos2] >= (int) c.size( ) )
               to2[j + Pos2] = (int) c.size( ) - 1;    }
     for ( int j = 0; j < (int) b1.size( ); j++ )
          ForceAssertGe( to1[j], 0 );
     for ( int j = 0; j < (int) b2.size( ); j++ )
          ForceAssertGe( to2[j], 0 );    }

void TryMergeMtigs( const assembly& a, int m1, int m2, const align& al,
     vec<read_location>& new_tig, basevector& new_tig_bases )
{    qualvector q;
     static vec<int> to1, to2, m1reads;
     MergeMtigsHead( a.mtig[m1], a.mtig[m2], a.mtig_qual[m1], a.mtig_qual[m2],
          al, new_tig_bases, q, to1, to2 );
     const basevector &b1 = a.mtig[m1], &b2 = a.mtig[m2];
     new_tig.clear( ), m1reads.clear( );
     for ( int i = 0; i < a.reads_orig_index[m1].isize( ); i++ )
     {    read_location r = a.reads_orig[a.reads_orig_index[m1][i]];
          m1reads.push_back( r.ReadId( ) );
          int start = Max( 0, r.Start( ) );
          start = Min( start, (int) new_tig_bases.size( ) - 1 );
          r.SetStartOnContig( to1[start] );
          r.SetLengthOfContig( new_tig_bases.size( ) );
          new_tig.push_back(r);    }
     UniqueSort(m1reads);
     for ( int i = 0; i < a.reads_orig_index[m2].isize( ); i++ )
     {    read_location r = a.reads_orig[ a.reads_orig_index[m2][i] ];
          if ( BinMember( m1reads, r.ReadId( ) ) ) continue;
          int start = Max( 0, r.Start( ) );
          start = Min( start, (int) to2.size( ) - 1 );
          r.SetStartOnContig( to2[start] );
          r.SetLengthOfContig( new_tig_bases.size( ) );
          r.SetContig(m1);
          new_tig.push_back(r);    }    }

void MergeMtigs( assembly& a, int m1, int m2, align& al )
{
     basevector c;
     qualvector q;
     static vec<int> to1, to2;
     MergeMtigsHead( a.mtig[m1], a.mtig[m2], a.mtig_qual[m1], a.mtig_qual[m2],
          al, c, q, to1, to2 );
     basevector &b1 = a.mtig[m1], &b2 = a.mtig[m2];
     qualvector &q1 = a.mtig_qual[m1], &q2 = a.mtig_qual[m2];

     // Update mtig_aligns.  We do the simplest thing, which is to kill
     // all relevant entries.

     a.mtig_aligns_lost[m1] = True;
     a.mtig_aligns_lost[m2] = True;

     // Housekeep.

     c.Swap( b1 );
     q.Swap( q1 );
     a.ClearMtig(m2);

     // Build read locations.

     static vec<int> m1reads;
     m1reads.clear( );
     for ( unsigned int i = 0; i < a.reads_index[m1].size( ); i++ )
     {    read_location& r = a.reads[ a.reads_index[m1][i] ];
          m1reads.push_back( r.ReadId( ) );
          int start = Max( 0, r.Start( ) );
          start = Min( start, (int) b1.size( ) - 1 );
          r.SetStartOnContig( to1[start] );
          r.SetLengthOfContig( b1.size( ) );    }
     UniqueSort(m1reads);
     for ( unsigned int i = 0; i < a.reads_index[m2].size( ); i++ )
     {    read_location& r = a.reads[ a.reads_index[m2][i] ];
          if ( BinMember( m1reads, r.ReadId( ) ) ) r.Kill( );
          else
          {    int start = Max( 0, r.Start( ) );
               start = Min( start, (int) to2.size( ) - 1 );
               r.SetStartOnContig( to2[start] );
               r.SetLengthOfContig( b1.size( ) );
               r.SetContig(m1);
               r.ForceInBounds( );
               a.reads_index[m1].push_back( a.reads_index[m2][i] );    }    }
     a.reads_index[m2].clear( );
     m1reads.clear( );
     for ( unsigned int i = 0; i < a.reads_orig_index[m1].size( ); i++ )
     {    read_location& r = a.reads_orig[a.reads_orig_index[m1][i]];
          m1reads.push_back( r.ReadId( ) );
          int start = Max( 0, r.Start( ) );
          start = Min( start, (int) b1.size( ) - 1 );
          r.SetStartOnContig( to1[start] );
          r.SetLengthOfContig( b1.size( ) );    }
     UniqueSort(m1reads);
     for ( unsigned int i = 0; i < a.reads_orig_index[m2].size( ); i++ )
     {    read_location& r = a.reads_orig[ a.reads_orig_index[m2][i] ];
          if ( BinMember( m1reads, r.ReadId( ) ) ) r.Kill( );
          else
          {    int start = Max( 0, r.Start( ) );
               start = Min( start, (int) to2.size( ) - 1 );
               r.SetStartOnContig( to2[start] );
               r.SetLengthOfContig( b1.size( ) );
               r.SetContig(m1);
               a.reads_orig_index[m1].push_back(
                    a.reads_orig_index[m2][i] );    }    }
     a.reads_orig_index[m2].clear( );    }

void RecomputeGaps( assembly& a, int supercontig_id, Bool reorder, Bool verbose )
{    super& s = a.supers[supercontig_id];
     int n = s.mtig.size( );
     vec<int> pos(n);
     vec<Bool> defined( n, False );
     pos[0] = 0;
     defined[0] = True;
     Bool progress = True;
     while(progress)
     {    progress = False;
          for ( int j = 1; j <= n-1; j++ )
          {    for ( int i = 0; i < n - j; i++ )
               {    int m1 = s.mtig[i], m2 = s.mtig[i+j];
                    int len1 = a.mtig[m1].size( );
                    if ( defined[i] != defined[i+j] && Connected( a, m1, m2 ) )
                    {    progress = True;
                         int g = ExpectedGap2( a, m1, m2 );
                         if ( defined[i] )
                         {    pos[i+j] = pos[i] + len1 + g;
                              defined[i+j] = True;    }
                         else
                         {    pos[i] = pos[i+j] - len1 - g;
                              defined[i] = True;    }
                         if ( j >= 2 ) goto reset_j;    }    }    }
          reset_j: continue;    }
     for ( int j = 0; j < n; j++ )
     {    if ( !defined[j] )
          {    if (verbose)
               {    cout << "Warning: In supercontig " << supercontig_id
                         << ", no evidence found to link " << "mtig " << s.mtig[0]
                         << " to mtig " << s.mtig[j] << "\n";
                    cout << "Value will be estimated based on layout "
                         << "information.\n";
                    cout << "It is likely that this problem arose because "
                         << "consensus "
                         << "sequence was not found for a contig.\n";    }
               int estimated_pos_j = 0;
               for ( int i = 0; i < j; i++ )
                    estimated_pos_j += a.mtig[ s.mtig[i] ].size( ) + s.gap[i];
               pos[j] = estimated_pos_j;
               defined[j] = True;    }    }

     for ( int j = 1; j < n; j++ )
     {    if ( pos[j] < pos[j-1] )
          {    int m1 = s.mtig[j], m2 = s.mtig[j-1];
               if (verbose)
               {    cout << "In supercontig " << supercontig_id
                         << ", it looks like mtig " << m1
                         << " should be placed before mtig " << m2 << "\n";
                    cout << "(by " << pos[j-1] - pos[j]
                         << " bases).  The lengths " << "of these mtigs are "
                         << a.mtig[m1].size( ) << " and " << a.mtig[m2].size( )
                         << ", respectively.\n";    }    }    }

     // Reorder the mtigs in the supercontig.

     if (reorder)
     {    vec< vec<int> > pos_id;
          for ( int j = 0; j < n; j++ )
          {    vec<int> v(2);
               v[0] = pos[j];
               v[1] = s.mtig[j];
               pos_id.push_back(v);    }
          Sort(pos_id);
          for ( int j = 0; j < n; j++ )
          {    pos[j] = pos_id[j][0];
               s.mtig[j] = pos_id[j][1];    }    }

     // Compute gaps, then do further reordering, with the goal of minimizing
     // the sum of the absolute values of the negative gaps between contigs in
     // the supercontig, but only do pairwise swaps to try to achieve this goal.

     for ( int j = 0; j < n-1; j++ )
          s.gap[j] = pos[j+1] - pos[j] - a.mtig[ s.mtig[j] ].size( );

     int neg_gap_total = 0;
     for ( int j = 0; j < n-1; j++ )
          if ( s.gap[j] < 0 ) neg_gap_total -= s.gap[j];

     if (reorder)
     {    int neg_gap_total_old = neg_gap_total + 1;
          while( neg_gap_total_old > neg_gap_total )
          {    neg_gap_total_old = neg_gap_total;
               for ( int l = 0; l < n-2; l++ )
               {    swap( pos[l], pos[l+1] );
                    swap( s.mtig[l], s.mtig[l+1] );
                    for ( int j = 0; j < n-1; j++ )
                         s.gap[j]
                              = pos[j+1] - pos[j] - a.mtig[ s.mtig[j] ].size( );
                    neg_gap_total = 0;
                    for ( int j = 0; j < n-1; j++ )
                         if ( s.gap[j] < 0 ) neg_gap_total -= s.gap[j];
                    if ( neg_gap_total < neg_gap_total_old )
                    {    if (verbose)
                         {    cout << "swapping positions of mtigs " << s.mtig[l]
                                   << " and " << s.mtig[l+1] << "\n";    }
                         break;    }
                    neg_gap_total = neg_gap_total_old;
                    swap( pos[l], pos[l+1] );
                    swap( s.mtig[l], s.mtig[l+1] );    }    }    }

     for ( int j = 0; j < n; j++ )
          a.mtigs_to_super_pos[ s.mtig[j] ] = j;

     for ( int j = 0; j < n-1; j++ )
          s.gap[j] = pos[j+1] - pos[j] - a.mtig[ s.mtig[j] ].size( );    }

void RecomputeGaps( assembly& a, Bool reorder, Bool verbose )
{    for ( int k = 0; k < (int) a.supers.size( ); k++ )
          RecomputeGaps( a, k, reorder, verbose );    }

void PrintLinkAccuracy( assembly& a )
{    cout << "\nlink accuracy between contigs within supercontigs"
          << " having more than one contig\n";
     for ( int k = 0; k < (int) a.supers.size( ); k++ )
     {    super& s = a.supers[k];
          int n = s.mtig.size( );
          if ( n == 1 ) continue;
          cout << "\nsupercontig " << k << "\n";

          // Compute contig positions.

          vec<int> len(n), pos(n);
          for ( int i = 0; i < n; i++ )
               len[i] = a.mtig[ s.mtig[i] ].size( );
          pos[0] = 0;
          for ( int i = 1; i < n; i++ )
               pos[i] = pos[i-1] + len[i-1] + s.gap[i-1];

          // Find all pairs of contigs which are directly linked.
          // Compute expected gap and s.d. for it.  How close to positions?

          float max_dev = 0;
          for ( int i = 0; i < n; i++ )
          {    int m1 = s.mtig[i];
               for ( int j = i + 1; j < n; j++ )
               {    int m2 = s.mtig[j];
                    if ( Connected( a, m1, m2, True ) )
                    {    int gap_ave, gapdev_ave;
                         Expect( a, m1, m2, gap_ave, gapdev_ave, True );
                         int given_gap = pos[j] - pos[i] - len[i];
                         float dev = float(given_gap - gap_ave) / float(gapdev_ave);
                         max_dev = Max( max_dev, dev );
                         cout << i+1
                              << " (mtig " << m1 << ") to " << j+1
                              << " (mtig " << m2 << ") of " << n
                              << ", given gap = " << given_gap
                              << ", " << setprecision(3)
                              << float(given_gap - gap_ave) / float(gapdev_ave)
                              << " dev's from predicted\n";    }    }    }
          cout << "maximum deviation: " << max_dev << "\n";    }
     cout << "\n";    }

void Reverse( assembly& a, int s )
{    super& sup = a.supers[s];
     for ( unsigned int i = 0; i < sup.mtig.size( ); i++ )
     {    int m = sup.mtig[i];
          ReverseMtig( a, m );
          a.mtigs_to_super_pos[m] = sup.mtig.size( ) - i - 1;    }
     ReverseThis( sup.mtig );
     ReverseThis( sup.gap );    }

void SpreadOverHoles( assembly& A, int s, vec<int>& spread )
{    super& sup = A.supers[s];
     int n = sup.mtig.size( );
     static vec<int> min_left_contig, max_left_contig, min_right_contig,
          max_right_contig, min_left_pos, max_left_pos, min_right_pos, max_right_pos;
     const int ALMOST_INFINITY = 1000 * 1000 * 1000;
     min_left_contig.resize_and_set( n-1, ALMOST_INFINITY );
     max_left_contig.resize_and_set( n-1, -1 );
     min_right_contig.resize_and_set( n-1, ALMOST_INFINITY );
     max_right_contig.resize_and_set( n-1, -1 );
     min_left_pos.resize_and_set( n-1, ALMOST_INFINITY );
     max_left_pos.resize_and_set( n-1, -1 );
     min_right_pos.resize_and_set( n-1, ALMOST_INFINITY );
     max_right_pos.resize_and_set( n-1, -1 );
     spread.resize_and_set( n-1, 0 );
     for ( int p1 = 0; p1 < (int) sup.mtig.size( ); p1++ )
     {    int m1 = sup.mtig[p1];
          int s1 = s;
          const vec<int>& m1reads = A.reads_orig_index[m1];
          for ( unsigned int i = 0; i < m1reads.size( ); i++ )
          {    const read_location& r1 = A.reads_orig[ m1reads[i] ];
               if ( r1.OrientationOnContig( ) != ForwardOr ) continue;
               int id1 = r1.ReadId( );
               if ( A.pairs_index[id1] >= 0 )
               {    const read_pairing& p = A.pairs[ A.pairs_index[id1] ];
                    int id2 = p.Partner(id1);
                    int si = A.simple_reads_orig_index[id2];
                    if ( si >= 0 )
                    {    const read_location& r2 = A.reads_orig[si];
                         if ( r2.OrientationOnContig( ) != ReverseOr ) continue;
                         int m2 = r2.Contig( );
                         int s2 = A.mtigs_to_supers[m2];
                         if ( s1 != s2 ) continue;
                         int p2 = A.mtigs_to_super_pos[m2];
                         for ( int j = p1; j < p2; j++ )
                         {    if ( p1 < min_left_contig[j] )
                              {    min_left_contig[j] = p1;
                                   min_left_pos[j] = r1.StartOnContig( );   }
                              else if ( p1 == min_left_contig[j] )
                                   min_left_pos[j]
                                        = Min( min_left_pos[j],
                                             r1.StartOnContig( ) );
                              if ( p1 > max_left_contig[j] )
                              {    max_left_contig[j] = p1;
                                   max_left_pos[j] = r1.StopOnContig( );   }
                              else if ( p1 == max_left_contig[j] )
                                   max_left_pos[j]
                                        = Max( max_left_pos[j], r1.StopOnContig( ) );
                              if ( p2 < min_right_contig[j] )
                              {    min_right_contig[j] = p2;
                                   min_right_pos[j] = r2.StartOnContig( );   }
                              else if ( p2 == min_right_contig[j] )
                                   min_right_pos[j]
                                        = Min( min_right_pos[j],
                                             r2.StartOnContig( ) );
                              if ( p2 > max_right_contig[j] )
                              {    max_right_contig[j] = p2;
                                   max_right_pos[j] = r2.StopOnContig( );   }
                              else if ( p2 == max_right_contig[j] )
                                   max_right_pos[j] = Max( max_right_pos[j],
                                        r2.StopOnContig( ) );
                                        }    }    }    }    }
     for ( int i = 0; i < n-1; i++ )
     {    if ( min_left_contig[i] == ALMOST_INFINITY ) continue;
          if ( min_right_contig[i] == ALMOST_INFINITY ) continue;
          int left_spread, right_spread;
          if ( min_left_contig[i] == max_left_contig[i] )
               left_spread = max_left_pos[i] - min_left_pos[i];
          else
          {    left_spread =
                    A.mtig[ sup.mtig[ min_left_contig[i] ] ].size( )
                         - min_left_pos[i] + max_left_pos[i];
               for ( int j = min_left_contig[i] + 1; j < max_left_contig[i]; j++ )
                    left_spread += A.mtig[ sup.mtig[j] ].size( );    }
          if ( min_right_contig[i] == max_right_contig[i] )
               right_spread = max_right_pos[i] - min_right_pos[i];
          else
          {    right_spread =
                    A.mtig[ sup.mtig[ min_right_contig[i] ] ].size( )
                         - min_right_pos[i] + max_right_pos[i];
               for ( int j = min_right_contig[i] + 1; j < max_right_contig[i]; j++ )
                    right_spread += A.mtig[ sup.mtig[j] ].size( );    }
          spread[i] = Min( left_spread, right_spread );    }    }

bool IsSuperConnected( assembly &A, int super_id )
{
  const vec<read_location> &locs = A.reads_orig;

  vec_inserts links( &A );
  links.AddContigs( A.supers[super_id].mtig );

  vec<int> orbit;
  equiv_rel the_equivalence( A.supers[super_id].mtig.size( ) );
  map<int, int>::iterator iter;

  for (int ii=0; ii<(int)links.size( ); ii++) {
    const insert_ends &ins = links[ii];
    int cg1 = locs[ ins.Loc1( ) ].Contig( );
    int cg2 = locs[ ins.Loc2( ) ].Contig( );
    if ( cg1 == cg2 )
      continue;

    int cg_pos1 = -1;
    int cg_pos2 = -1;

    iter = A.mtigs_to_super_pos.find( cg1 );
    if ( iter != A.mtigs_to_super_pos.end( ) )
      cg_pos1 = iter->second;
    iter = A.mtigs_to_super_pos.find( cg2 );
    if ( iter != A.mtigs_to_super_pos.end( ) )
      cg_pos2 = iter->second;

    ForceAssert( cg_pos1 > -1 && cg_pos2 > -1 );

    the_equivalence.Join( cg_pos1, cg_pos2 );
  }

  the_equivalence.Orbit( 0, orbit );

  return ( orbit.size( ) == A.supers[super_id].mtig.size( ) );
}

bool IsContigFlimsy( assembly &A, int contig_id )
{
  vec<orientation> orient;

  // Check paired reads on contig.
  const vec<int> &locs_index = A.reads_orig_index[ contig_id ];

  for (int ii=0; ii<(int)locs_index.size( ); ii++) {
    int read_id = A.reads_orig[ locs_index[ii] ].ReadId( );
    if ( A.pairs_index[ read_id ] < 0 )
      continue;
    orientation ori = A.reads_orig[ locs_index[ii] ].OrientationOnContig( );

    if ( orient.size( ) < 1 ) {
      orient.push_back( ori );
      continue;
    }

    if ( orient[0] != ori )
      return false;
  }

  // It is flimsy.
  return true;
}

bool IsSuperDegenerate( assembly &A, int super_id, int sc_info_min_length )
{
  // If gapped length of super is big enough.
  int super_length = A.SuperLen( super_id );
  if ( super_length > sc_info_min_length )
    return false;

  // If super contains at least one non-flimsy contig.
  for (int ii=0; ii<(int)A.supers[ super_id ].mtig.size( ); ii++) {
    int contig_id = A.supers[ super_id ].mtig[ ii ];
    if ( !IsContigFlimsy( A, contig_id ) )
      return false;
  }

  // It is degenerate.
  return true;
}

