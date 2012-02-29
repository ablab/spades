///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include "Alignment.h"
#include "Basevector.h"
#include "math/Functions.h"
#include "IndexedAlignmentPlusVector.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "Set.h"
#include "ShortVector.h"
#include "system/System.h"
#include "system/Types.h"
#include "Vec.h"

void alignment::Write( ostream& out ) const
{    int pos1, pos2;
     avector<int> gaps, lengths;
     packalign::Unpack( pos1, pos2, gaps, lengths );
     BinWrite( out, pos1 );
     BinWrite( out, pos2 );
     BinWrite( out, errors_ );
     BinWrite( out, lengths.length );
     for ( unsigned int i = 0; i < lengths.length; i++ )
     {    BinWrite( out, gaps(i) );
          BinWrite( out, lengths(i) );    }    }

void alignment::Read( istream& in )
{    int pos1, pos2, errors;
     BinRead( in, pos1 );
     BinRead( in, pos2 );
     BinRead( in, errors );
     int nblocks;
     BinRead( in, nblocks );
     avector<int> gaps(nblocks), lengths(nblocks);
     for ( int i = 0; i < nblocks; i++ )
     {    BinRead( in, gaps(i) );
          BinRead( in, lengths(i) );    }
     Set( pos1, pos2, errors, gaps, lengths );    }

void alignment::Write( ostream& out, int id1, int id2, Bool rc )
{    Write(out);
     BinWrite( out, id1 );
     BinWrite( out, id2 );
     BinWrite( out, rc );    }

void alignment::Read( istream& in, int& id1, int& id2, Bool& rc )
{    Read(in);
     BinRead( in, id1 );
     BinRead( in, id2 );
     BinRead( in, rc );    }

void alignment_plus::Write( ostream& out ) const
{    int id1 = Id1( ), id2 = Id2( );
     BinWrite( out, id1 );
     BinWrite( out, id2 );
     BinWrite( out, score );
     Bool rc2 = Rc2( );
     BinWrite( out, rc2 );
     int pos1, pos2;
     avector<int> gaps, lengths;
     a.packalign::Unpack( pos1, pos2, gaps, lengths );
     BinWrite( out, pos1 );
     BinWrite( out, pos2 );
     int e = a.Errors( );
     BinWrite( out, e );
     BinWrite( out, lengths.length );
     for ( unsigned int i = 0; i < lengths.length; i++ )
     {    BinWrite( out, gaps(i) );
          BinWrite( out, lengths(i) );    }    }

void alignment_plus::Read( istream& in )
{    int id1, id2;
     Bool rc2;
     BinRead( in, id1 );
     BinRead( in, id2 );
     BinRead( in, score );
     BinRead( in, rc2 );
     SetId1(id1);
     SetId2(id2);
     SetRc2(rc2);
     int pos1, pos2, errors;
     BinRead( in, pos1 );
     BinRead( in, pos2 );
     BinRead( in, errors );
     int nblocks;
     BinRead( in, nblocks );
     avector<int> gaps(nblocks), lengths(nblocks);
     for ( int i = 0; i < nblocks; i++ )
     {    BinRead( in, gaps(i) );
          BinRead( in, lengths(i) );    }
     a.Set( pos1, pos2, errors, gaps, lengths );    }

void alignment::Compactify( int len1, int len2 )
{    
     int pos1, pos2, errors, ni;
     avector<int> gaps, lengths;
     Unpack( pos1, pos2, errors, gaps, lengths );

     // Clean up alignments that start with gaps.  This should be very rare.
     // The goal here is just to get a valid alignment, not necessarily a good one.

     int first_gap;
     check_first_gap:
     ni = 0;
     if ( gaps.length == 0 ) goto finish;
     first_gap = gaps(0);
     if ( first_gap < 0 )
     {    first_gap = -first_gap; // but on second read
          if ( first_gap <= (int) pos2 )
          {    pos2 -= first_gap;
               lengths(0) += first_gap;
               gaps(0) = 0;    }
          else
          {    first_gap -= pos2;
               gaps(0) = -first_gap;
               if ( pos2 > 0 )
               {    lengths.Prepend(pos2);
                    gaps.Prepend(0);    }
               pos2 = 0;    }    }
     else if ( first_gap > 0 )
     {    if ( first_gap <= (int) pos1 )
          {    pos1 -= first_gap;
               lengths(0) += first_gap;
               gaps(0) = 0;    }
          else
          {    first_gap -= pos1;
               gaps(0) = first_gap;
               if ( pos1 > 0 )
               {    lengths.Prepend(pos1);
                    gaps.Prepend(0);    }
               pos1 = 0;    }    }

     // Remove an initial length of zero, if it exists.

     if ( lengths(0) == 0 )
     {    for ( unsigned int i = 1; i < gaps.length; i++ )
          {    gaps(i-1) = gaps(i);
               lengths(i-1) = lengths(i);    }
          gaps.resize( gaps.length - 1 );
          lengths.resize( lengths.length - 1 );
          goto check_first_gap;    }

     for ( unsigned int i = 1; i < gaps.length; i++ )
     {    if ( lengths(ni) != 0 && gaps(i) != 0 )
          {    ++ni;
               gaps(ni) = gaps(i);
               lengths(ni) = lengths(i);    }
          else
          {    lengths(ni) += lengths(i);
               if ( gaps(ni) * gaps(i) < 0 )
                    lengths(ni) += Min( Abs(gaps(ni)), Abs(gaps(i)) );
               gaps(ni) += gaps(i);    }    }

     if ( lengths(ni) == 0 ) 
     {    if ( ni == 0 )
          {    gaps.resize(0);
               lengths.resize(0);
               Set( pos1, pos2, errors, gaps, lengths );
               return;    }
          if ( gaps(ni) < 0 )
          {    int Pos2 = pos2;
               for ( int j = 0; j < ni; j++ )
               {    if ( gaps(j) > 0 ) Pos2 += gaps(j);
                    Pos2 += lengths(j);    }
               int npos2 = len2 - Pos2;
               Assert( npos2 >= 0 ); // XXX
               if ( npos2 == 0 ) --ni;
               else if ( npos2 >= -gaps(ni) )
               {    lengths(ni-1) -= gaps(ni);
                    --ni;    }
               else
               {    gaps(ni) += npos2;
                    lengths(ni) = npos2;    }    }
          else 
          {    int Pos1 = pos1;
               for ( int j = 0; j < ni; j++ )
               {    if ( gaps(j) < 0 ) Pos1 -= gaps(j);
                    Pos1 += lengths(j);    }
               int npos1 = len1 - Pos1;
               Assert( npos1 >= 0 ); // XXX
               if ( npos1 == 0 ) --ni;
               else if ( npos1 >= gaps(ni) )
               {    lengths(ni-1) += gaps(ni);
                    --ni;    }
               else
               {    gaps(ni) -= npos1;
                    lengths(ni) = npos1;    }    }    }

     gaps.resize(ni+1);
     lengths.resize(ni+1);
     if ( gaps(0) > 0 ) // shouldn't happen
     {    pos2 += gaps(0);
          gaps(0) = 0;    }
     else if ( gaps(0) < 0 ) // shouldn't happen
     {    pos1 -= gaps(0);
          gaps(0) = 0;    }

     finish:

     // Make alignments proper if only off by one base at the beginning.
     // (This is really covering up a bug in some other piece of code.)

     if ( (pos1 == 1 && pos2 > 0) || (pos1 > 0 && pos2 == 1) )
     {    --pos1;
          --pos2;
          ++lengths(0);    }
     
     Bool zero_gap = False;
     for ( unsigned int i = 1; i < gaps.length; i++ )
          if ( gaps(i) == 0 ) zero_gap = True;
     if (zero_gap) 
     // {    cout << "zero gap found, repeating\n"; // XXX
          // PRINT2( pos1, pos2 ); // XXX
          // cout << "gaps/lengths:"; // XXX
          // for ( unsigned int i = 0; i < gaps.length; i++ ) // XXX
          //     cout << " " << gaps(i) << "/" << lengths(i); // XXX
          // cout << "\n"; // XXX
          goto check_first_gap;    
          //     } // XXX

     Set( pos1, pos2, errors, gaps, lengths );    
     // cout << "leaving compactify" << endl; // ZZZ
          }

void alignment_plus::Print( Bool abbreviate, ostream& out, 
     const basevector& rd1, const basevector& rd2, const basevector& rd2rc, 
     const qualvector& q1, const qualvector& q2, const qualvector& q2rc, 
     int begin, Bool one_frame, Bool highlight_score )
{    out << "\nAlignment between sequences " << Id1( ) << " and ";
     if (Rc2( )) out << "rc of ";
     out << Id2( ) << ".\n\n";
     out << a.pos1( ) << " --> " << a.Pos1( ) << " (of "
          << rd1.size( ) << ")      ";
     out << a.pos2( ) << " --> " << a.Pos2( ) << " (of " << rd2.size( ) << ")\n";
     out << "\n";
     out << "score for this alignment = " << (highlight_score ? START_GREEN : "")
          << score << (highlight_score ? END_ESCAPE : "" ) << "\n";
     if ( !Rc2( ) ) 
          PrintVisualAlignment( abbreviate, out, rd1, rd2, a, q1, q2, begin,
               one_frame );
     else PrintVisualAlignment( abbreviate, out, rd1, rd2rc, a, q1, q2rc, begin,
               one_frame );    }

void alignment_plus::Print( Bool abbreviate, ostream& out, const basevector& rd1, 
          const basevector& rd2, const qualvector& q1, const qualvector& q2, 
          int begin, Bool one_frame, Bool highlight_score )
{    out << "\nAlignment between reads " << Id1( ) << " and ";
     if (Rc2( )) out << "rc of ";
     out << Id2( ) << ".\n\n";
     out << a.pos1( ) << " --> " << a.Pos1( ) << " (of "
          << rd1.size( ) << ")      ";
     out << a.pos2( ) << " --> " << a.Pos2( ) << " (of " << rd2.size( ) << ")\n";
     out << "\n";
     out << "score for this alignment = " << (highlight_score ? START_GREEN : "")
          << score << (highlight_score ? END_ESCAPE : "" ) << "\n";
     if ( !Rc2( ) ) 
          PrintVisualAlignment( abbreviate, out, rd1, rd2, a, q1, q2, begin,
               one_frame );
     else 
     {    basevector rd2rc;
          rd2rc = rd2;
          rd2rc.ReverseComplement( );
          qualvector q2rc;
          q2rc = q2;
          ReverseThis(q2rc);
          PrintVisualAlignment( abbreviate, out, rd1, rd2rc, a, q1, q2rc, begin,
               one_frame );    }    }

void alignment_plus::Swap( int rd1length, int rd2length )
{    a.Flip( );
     Bool rc2 = Rc2( );
     int temp = Id2( );
     SetId2( Id1( ) );
     SetId1(temp);
     SetRc2(rc2);
     if ( Rc2( ) ) a.ReverseThis( rd2length, rd1length );    }

void alignment_plus::SetToSwapOf( const alignment_plus& x, int rd1length, 
     int rd2length )
{    SetId2( x.Id1( ) );
     SetId1( x.Id2( ) );
     SetRc2( x.Rc2( ) );
     score = x.score;
     if ( Rc2( ) ) a.SetToReverseFlipOf( x.a, rd2length, rd1length );
     else a.SetToFlipOf( x.a );
     a.SetErrors( x.a.Errors( ) );    }

void alignment_plus::SetToSwapOf( const align& x, int id1, int id2, Bool rc2,
     float s, int rd1length, int rd2length )
{    SetId2( id1 );
     SetId1( id2 );
     SetRc2( rc2 );
     score = s;
     if ( Rc2( ) ) a.SetToReverseFlipOf( x, rd2length, rd1length );
     else a.SetToFlipOf( x );    }

vector<int> alignment::MutationsGap1Gap2( const basevector& rd1, 
     const basevector& rd2 )
{    vector<int> answer(3, 0);
     avector<int> gaps, lengths;
     int pos1, pos2, errors;
     Unpack( pos1, pos2, errors, gaps, lengths );
     int p1 = pos1, p2 = pos2;
     for ( unsigned int j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) > 0 ) 
          {    answer[1] += gaps(j);
               p2 += gaps(j);    }
          if ( gaps(j) < 0 )
          {    answer[2] -= gaps(j);
               p1 -= gaps(j);    }
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( rd1[p1] != rd2[p2] ) ++answer[0];
               ++p1; 
               ++p2;    }    }
     return answer;    }

int alignment::Mutations( const basevector& rd1, const basevector& rd2,
     const qualvector& q1, int min_score )
{    int answer = 0;
     avector<int> gaps, lengths;
     int pos1, pos2, errors;
     Unpack( pos1, pos2, errors, gaps, lengths );
     unsigned int j, p1 = pos1, p2 = pos2;
     for ( j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) > 0 ) p2 += gaps(j);
          if ( gaps(j) < 0 ) p1 -= gaps(j);
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( rd1[p1] != rd2[p2] && q1[p1] >= min_score ) ++answer;
               ++p1; 
               ++p2;    }    }
     return answer;    }

int alignment::Indels( const basevector& rd1, const basevector& rd2,
     const qualvector& q1, int min_score )
{    int answer = 0;
     avector<int> gaps, lengths;
     int pos1, pos2, errors;
     Unpack( pos1, pos2, errors, gaps, lengths );
     int p1 = pos1, p2 = pos2;
     for ( unsigned int j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) > 0 ) 
          {    if ( p1+1 < (int) q1.size( ) && q1[p1] >= min_score 
                    && q1[p1+1] >= min_score ) answer += gaps(j);
               p2 += gaps(j);    }
          if ( gaps(j) < 0 ) 
          {    if ( -gaps(j) + p1 < (int) q1.size( ) )
               {    int k;
                    for ( k = 0; k <= -gaps(j); k++ )
                         if ( q1[ p1 + k ] < min_score ) break;
                    if ( k == -gaps(j) + 1 ) answer -= gaps(j);    }
               p1 -= gaps(j);    }
          p1 += lengths(j);
          p2 += lengths(j);    }
     return answer;    }

vector<int> alignment::MutationsGap1Gap2( const basevector& rd1, 
     int from1, int to1, const basevector& rd2, int from2, int to2 )
{    vector<int> answer(3, 0);
     avector<int> gaps, lengths;
     int pos1, pos2, errors;
     Unpack( pos1, pos2, errors, gaps, lengths );
     int p1 = pos1, p2 = pos2;
     for ( unsigned j = 0; j < gaps.length; ++j ) 
     {    if ( gaps(j) > 0 ) 
          {    if ( p1 >= from1 && p1 < to1 && p2 >= from2 && p2 < to2 )
	            answer[1] += gaps(j);
               p2 += gaps(j);    }
          if ( gaps(j) < 0 ) 
          {    if ( p1 >= from1 && p1 < to1 && p2 >= from2 && p2 < to2 )
	            answer[2] -= gaps(j);
               p1 -= gaps(j);    }
          for ( int x = 0; x < lengths(j); ++x ) 
          {     if ( p1 >= from1 && p1 < to1 && p2 >= from2 && p2 < to2 &&
	             rd1[p1] != rd2[p2] )
	             ++answer[0];
                ++p1;
                ++p2;    }    }
     return answer;    }

int CorrelatePositions( const alignment& a, int x1 )
{    int pos1, pos2, errors;
     avector<int> gaps, lengths;
     a.Unpack( pos1, pos2, errors, gaps, lengths );
     if ( x1 < pos1 ) return OffTheEnd;
     if ( x1 == pos1 ) return pos2;
     for ( unsigned int j = 0; j < lengths.length; j++ )
     {    if ( gaps(j) > 0 ) pos2 += gaps(j);
          if ( gaps(j) < 0 )
               for ( int x = 0; x < -gaps(j); x++ )
                    if ( x1 == pos1++ ) return AtGap;
          if ( x1 < pos1 + lengths(j) ) return pos2 + x1 - pos1;
          else
          {    pos1 += lengths(j);
               pos2 += lengths(j);    }    }
     return OffTheEnd;    }

int ErrorsAt( const alignment& a, int x1, const basevector& rd1, 
     const basevector& rd2 )
{    int pos1, pos2, errors;
     avector<int> gaps, lengths;
     a.Unpack( pos1, pos2, errors, gaps, lengths );
     if ( x1 <= pos1 ) return 0;
     for ( unsigned int j = 0; j < lengths.length; j++ )
     {    if ( gaps(j) > 0 ) 
          {    if ( x1 == pos1 ) return gaps(j);
               pos2 += gaps(j);    }
          if ( gaps(j) < 0 )
               for ( int x = 0; x < -gaps(j); x++ )
                    if ( x1 == pos1++ ) return 1;
          if ( x1 < pos1 + lengths(j) ) 
          {    if ( rd1[x1] == rd2[ pos2 + x1 - pos1 ] ) return 0;
               else return 1;     }
          pos1 += lengths(j);
          pos2 += lengths(j);    }
     return 0;    }

alignment_plus::alignment_plus( int read_id1,
				int read_id2,
				int rd1length,
				int rd2length,
				Bool if_read2_is_rc,
				const alignment& a_arg,
				float score_arg )
  : score(score_arg),
    a(a_arg),
    read_id1_(read_id1),
    read_id2_(read_id2)
{
     if (if_read2_is_rc) read_id2_ = -read_id2_ - 1;

     /*
     // Perform a rudamentary test for integrity of the alignment.

     int fail = 0;

     int p1, p2;
     int errors;
     avector<int> gaps, lengths;
     a.Unpack( p1, p2, errors, gaps, lengths );

     if ( (int) p1 > rd1length ) fail = 1;
     if ( (int) p2 > rd2length ) fail = 2;

     if ( gaps.length == 0 ) fail = 3;
     if ( gaps(0) != 0 ) fail = 4;
     if ( lengths( lengths.length - 1 ) == 0 ) fail = 5;

     if (fail)
     {    cout << "\nalignment_plus constructor: alignment found which fails test "
               << fail << "\n";
          PRINT( Id1( ) ); PRINT( Id2( ) ); PRINT(p1); PRINT(p2); PRINT(rd1length);
          PRINT(rd2length); PRINT(int(Rc2( ))); PRINT(gaps.length);
          cout << "gaps/lengths:";
          for ( unsigned int i = 0; i < gaps.length; i++ )
               cout << " " << gaps(i) << "/" << lengths(i);
          cout << endl;
          ForceAssert( 0 == 1 );    }    
     */

          }

istream& operator>>( istream& in, vec<alignment_plus>& vecAligns )
{    
  VecAlignmentPlusReader( &in ).ReadAll( vecAligns );
  return in;
}

ostream& operator<<( ostream& out, const vec<alignment_plus>& vecAligns )
{    
  VecAlignmentPlusWriter( &out ).Write( vecAligns );
  return out;
}

void WriteAppend( const String &filename, const vec<alignment_plus> &vecAligns )
{
  const bool createIndex = false;
  VecAlignmentPlusWriter( filename, createIndex ).Append( vecAligns );
}

void WriteAlignsWithIndex( const String &filename, const vec<alignment_plus> &vecAligns )
{
  VecAlignmentPlusWriter( filename ).Write( vecAligns );
}

void WriteAppendWithIndex( const String &filename, const vec<alignment_plus> &vecAligns )
{
  VecAlignmentPlusWriter( filename ).Append( vecAligns );
}

// DepthOfCoverage: given a collection of reads, assembling to a contig of known
// length, determine the depth of coverage of this stretch by all reads which match 
// it well.

Float DepthOfCoverage( const vec<int>& reads, int contig_length,
     const vec<int>& read_lengths, const vec<alignment_plus>& all_aligns, 
     const vec<int>& all_aligns_index )
{    
     // First find all reads which match at least one of the given reads with score 
     // <= 100, and which is not in the original set of reads.

     set<int> original_reads, new_overlaps;
     for ( unsigned int j = 0; j < reads.size( ); j++ )
          original_reads.insert( reads[j] );
     for ( unsigned int j = 0; j < reads.size( ); j++ )
     {    int id = reads[j];
          int ind = all_aligns_index[id];
          for ( unsigned int k = ind; ind >= 0 && k < all_aligns.size( ); k++ )
          {    const alignment_plus& ap = all_aligns[k];
               if ( ap.Id1( ) > id ) break;
               if ( ap.score <= 100 && !Member( original_reads, ap.Id2( ) ) ) 
                    new_overlaps.insert( ap.Id2( ) );    }    }

     // For each of the reads, determine its leftmost and rightmost bases
     // which match a base in one of the given reads.  Compute the total number
     // of bases spanned.

     int total_bases = 0;
     for ( unsigned int i = 0; i < reads.size( ); i++ )
          total_bases += read_lengths[ reads[i] ];
     for ( set<int>::iterator i = new_overlaps.begin( ); i != new_overlaps.end( );
          i++ )
     {    int id = *i;
          int ind = all_aligns_index[id];
          int pos1 = 10000000, Pos1 = -10000000;
          for ( unsigned int k = ind; ind >= 0 && k < all_aligns.size( ); k++ )
          {    const alignment_plus& ap = all_aligns[k];
               if ( ap.Id1( ) > id ) break;
               if ( !Member( original_reads, ap.Id2( ) ) ) continue;
               pos1 = Min( pos1, (int) ap.a.pos1( ) );
               Pos1 = Max( Pos1, (int) ap.a.Pos1( ) );
               if ( pos1 == 0 && Pos1 == read_lengths[id] ) break;    }
          total_bases += Pos1 - pos1;    }

     return Float(total_bases)/Float(contig_length);    }

int FindPos1( int p1, int p2, const avector<int>& gaps, const avector<int>& lengths )
{    for ( unsigned int j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) < 0 ) p1 -= gaps(j);
          p1 += lengths(j);    }
     return p1;    }

int FindPos2( int p1, int p2, const avector<int>& gaps, const avector<int>& lengths )
{    for ( unsigned int j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) > 0 ) p2 += gaps(j);
          p2 += lengths(j);    }
     return p2;    }

// RequireProper tests an alignment_plus to see if it is proper.  It does not
// object to "dead" alignments, those with no gaps or lengths.  Return False
// if the alignment is dead or improper.

Bool RequireProper( const alignment_plus& ap, const vecbasevector& EE, 
     int test_no, Bool fatal )
{    int pos1, pos2, errors;
     avector<int> gaps, lengths;
     ap.a.Unpack( pos1, pos2, errors, gaps, lengths );
     if ( gaps.length == 0 ) return False;
     int Pos1 = FindPos1( pos1, pos2, gaps, lengths );
     int Pos2 = FindPos2( pos1, pos2, gaps, lengths );
     int id1 = ap.Id1( ), id2 = ap.Id2( );
     Bool zero_gap = False, bad_length = False;
     for ( unsigned int i = 1; i < gaps.length; i++ )
          if ( gaps(i) == 0 ) zero_gap = True;
     for ( unsigned int i = 0; i < lengths.length; i++ )
          if ( lengths(i) < 0 
               || ( i < lengths.length - 1 && lengths(i) == 0 ) ) bad_length = True;
     if ( zero_gap || bad_length || gaps(0) != 0 ||
          Pos1 > (int) EE[id1].size( ) || Pos2 > (int) EE[id2].size( ) ||
          !( (pos1 == 0 && Pos2 == (int) EE[id2].size( ))
          || (pos2 == 0 && Pos1 == (int) EE[id1].size( ))
          || (pos1 == 0 && Pos1 == (int) EE[id1].size( ))
          || (pos2 == 0 && Pos2 == (int) EE[id2].size( )) ) )
     {    cout << "\nWarning: test " << test_no << " for proper alignment failed\n";
          cout << "(This should not be a problem.)" << endl;
          cout << ( ap.Rc2( ) ? "reverse" : "forward" ) << " alignment\n";
          PRINT2( id1, id2 );
          PRINT2( pos1, Pos1 );
          PRINT2( pos2, Pos2 );
          PRINT2( EE[id1].size( ), EE[id2].size( ) );
          cout << "actual errors = ";
          if ( !ap.Rc2( ) ) 
               cout << ActualErrors( EE[id1], EE[id2], ap.a ) << "\n";
          else 
          {    basevector EErcid2 = EE[id2];
               EErcid2.ReverseComplement( );
               cout << ActualErrors( EE[id1], EErcid2, ap.a ) << "\n";    }
          cout << "gaps/lengths:";
          for ( unsigned int i = 0; i < lengths.length; i++ )
               cout << " " << gaps(i) << "/" << lengths(i);
          cout << "\n";
          if (fatal) ForceAssert( 0 == 1 );    
          return False;    }
     return True;    }

Bool RequireProper( const alignment_plus& ap, const vec<int>& EE_length,
     int test_no, Bool fatal )
{    int pos1, pos2, errors;
     avector<int> gaps, lengths;
     ap.a.Unpack( pos1, pos2, errors, gaps, lengths );
     if ( gaps.length == 0 ) return False;
     int Pos1 = FindPos1( pos1, pos2, gaps, lengths );
     int Pos2 = FindPos2( pos1, pos2, gaps, lengths );
     int id1 = ap.Id1( ), id2 = ap.Id2( );
     Bool zero_gap = False, bad_length = False;
     for ( unsigned int i = 1; i < gaps.length; i++ )
          if ( gaps(i) == 0 ) zero_gap = True;
     for ( unsigned int i = 0; i < lengths.length; i++ )
          if ( lengths(i) < 0 
               || ( i < lengths.length - 1 && lengths(i) == 0 ) ) bad_length = True;
     if ( zero_gap || bad_length || gaps(0) != 0 ||
          Pos1 > EE_length[id1] || Pos2 > EE_length[id2] ||
          !( (pos1 == 0 && Pos2 == EE_length[id2])
          || (pos2 == 0 && Pos1 == EE_length[id1])
          || (pos1 == 0 && Pos1 == EE_length[id1])
          || (pos2 == 0 && Pos2 == EE_length[id2]) ) )
     {    cout << "\nWarning: test " << test_no << " for proper alignment failed\n";
          cout << "(This should not be a problem.)" << endl;
          cout << ( ap.Rc2( ) ? "reverse" : "forward" ) << " alignment\n";
          PRINT2( id1, id2 );
          PRINT2( pos1, Pos1 );
          PRINT2( pos2, Pos2 );
          PRINT2( EE_length[id1], EE_length[id2] );
          cout << "gaps/lengths:";
          for ( unsigned int i = 0; i < lengths.length; i++ )
               cout << " " << gaps(i) << "/" << lengths(i);
          cout << "\n";
          if (fatal) ForceAssert( 0 == 1 );    
          return False;    }
     return True;    }

int TransitiveOffset( const alignment_plus& ap1, const alignment_plus& ap2,
     int rd1length, int rd2length, int rd3length )
{    int start = ap2.a.pos1( ) - ap2.a.pos2( ),
          stop = ap2.a.Pos1( ) + rd3length - ap2.a.Pos2( );
     if ( ap1.Rc2( ) ) start = rd2length - stop - 1;
     start = start - ap1.a.pos2( ) + ap1.a.pos1( );
     return start;    }

int Bandwidth( alignment& a )
{    const int add_to_bandwidth = 8; // heuristic
     int pos1, pos2;
     int errors;
     avector<int> gaps, lengths;
     a.Unpack( pos1, pos2, errors, gaps, lengths );
     int low = 0, high = 0, gap_total = 0;
     for ( unsigned int l = 0; l < gaps.length; l++ )
     {    gap_total += gaps(l);
          low = Min( low, gap_total );
          high = Max( high, gap_total );    }
     return Max( Abs(low), Abs(high) ) + add_to_bandwidth;    }

// There are versions of TrimAlignmentFront for class align and for 
// class alignment.  We should try to eliminate the one for class alignment.

void TrimAlignmentFront( alignment& a, int n )
{    avector<int> gaps, lengths;
     int pos1, pos2, errors;
     a.Unpack( pos1, pos2, errors, gaps, lengths );

     int pos1_start = pos1;
     unsigned int j;
     for ( j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) < 0 )
          {    pos1 -= gaps(j);
               if ( pos1 - pos1_start >= n )
               {    gaps(j) = 0;
                    break;    }    }
          if ( gaps(j) > 0 ) pos2 += gaps(j);
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( pos1 - pos1_start >= n ) 
               {    lengths(j) -= x;
                    gaps(j) = 0;
                    goto done;    }    
               ++pos1;
               ++pos2;    }
          gaps(j) = 0;    }
     done:

     if ( pos1 - pos1_start < n )
     {    cout << "Warning: TrimAlignmentFront found pos1 - pos1_start < n\n";
          PRINT2(pos1, pos1_start);
          PRINT(n);

          int p1, p2, e;
          avector<int> gapsx, lengthsx;
          a.Unpack( p1, p2, e, gapsx, lengthsx );

          PRINT2( a.pos1( ), a.Pos1( ) );
          PRINT2( a.pos2( ), a.Pos2( ) ); 
          cout << "gaps/lengths:"; 
          for ( unsigned int i = 0; i < lengthsx.length; i++ ) 
               cout << " " << gapsx(i) << "/" << lengthsx(i); 
          cout << "\n";     }
          
     if ( j > 0 )
     {    for ( unsigned int l = j; l < gaps.length; l++ )
          {    gaps(l-j) = gaps(l);
               lengths(l-j) = lengths(l);    }
          gaps.resize( gaps.length - j );
          lengths.resize( lengths.length - j );    }
     a = alignment( pos1, pos2, errors, gaps, lengths );    }


void TrimAlignmentFront( align& a, int n )
{    int pos1 = a.pos1( ), pos2 = a.pos2( );
     const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );
     int nblocks = a.Nblocks( );

     int pos1_start = pos1;
     int j;
     for ( j = 0; j < nblocks; j++ )
     {    if ( gaps(j) < 0 )
          {    pos1 -= gaps(j);
               if ( pos1 - pos1_start >= n )
               {    a.SetGap( j, 0 );
                    break;    }    }
          if ( gaps(j) > 0 ) pos2 += gaps(j);
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( pos1 - pos1_start >= n ) 
               {    a.AddToLength( j, -x );
                    a.SetGap( j, 0 );
                    goto done;    }    
               ++pos1;
               ++pos2;    }
          a.SetGap( j, 0 );    }
     done:

     if ( pos1 - pos1_start < n )
     {    cout << "Warning: TrimAlignmentFront found pos1 - pos1_start < n\n";
          PRINT2(pos1, pos1_start);
          PRINT(n);    }

     if ( j > 0 )
     {    for ( int l = j; l < nblocks; l++ )
          {    a.SetGap( l-j, gaps(l) );
               a.SetLength( l-j, lengths(l) );    }
          a.SetNblocks( nblocks - j );    }
     a.Setpos1(pos1);
     a.Setpos2(pos2);    }


void TrimAlignmentBack( align& a, int n )
{    int Pos1 = a.Pos1( ) - 1, Pos2 = a.Pos2( ) - 1;
     const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );
     int nblocks = a.Nblocks( );

     int Pos1_start = Pos1;
     int j;
     for ( j = nblocks-1; j >= 0; j-- )
     {    for ( int x = 0; x < lengths(j); x++ )
          {    if ( Pos1_start - Pos1 >= n ) 
               {    a.AddToLength( j, -x );
                    goto done;    }    
               --Pos1;
               --Pos2;    }
          a.SetLength( j, 0 );
          if ( gaps(j) < 0 ) Pos1 += gaps(j);
          if ( gaps(j) > 0 ) Pos2 -= gaps(j);
          a.SetGap( j, 0 );    }
     done:

     if ( j < 0 )
     {    cout << "Warning: TrimAlignmentBack trimmed entire alignment\n";
          PRINT2( Pos1, Pos1_start );
	  PRINT( n );    }

     a.SetNblocks( j+1 );    }


int MaxPerfectMatch( Bool rd1_is_rc, const align& a, const basevector& rd1, 
     const basevector& rd2 )
{    int answer = 0;
     int pos1 = a.pos1( ), pos2 = a.pos2( );
     const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );
     int nblocks = a.Nblocks( );
     int j, p1 = pos1, p2 = pos2;
     for ( j = 0; j < nblocks; j++ )
     {    if ( gaps(j) > 0 ) p2 += gaps(j);
          if ( gaps(j) < 0 ) p1 -= gaps(j);
          int count = 0;
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( !rd1_is_rc && rd1[p1] != rd2[p2] )
               {    answer = Max( answer, count );
                    count = 0;    }
               else if ( rd1_is_rc && 3 - rd1[ rd1.size( ) - p1 - 1 ] != rd2[p2] )
               {    answer = Max( answer, count );
                    count = 0;    }
               else ++count;
               ++p1; 
               ++p2;    }    
          answer = Max( answer, count );    }
     return answer;    }

// void WriteAppend( const String& f, const vec<alignment_plus>& v )
// {    ForceAssert( !IsRegularFile( f + ".gz" ) );
//      static longlong max_size_bound = longlong(10000000) * longlong(100000000);
//      if ( !IsRegularFile(f) )
//      {    ForceAssertLt( (longlong) v.size( ), max_size_bound );
//           Ofstream( out, f );
//           out << setfill( '0' ) << setw(15) << v.size( ) << "\n";
//           WriteAligns( out, v );    }
//      else
//      {    longlong n;
//           {    Ifstream( in, f );
//                in >> n;    }
//           ForceAssertLt( n + (longlong) v.size( ), max_size_bound );
//           ostrstream osize;
//           osize << setfill( '0' ) << setw(15) << n + v.size( ) << "\n";
//           int fd = Open( f, O_WRONLY );
//           WriteBytes( fd, osize.str( ), 16 );
//           close(fd);
//           ofstream out( f.c_str( ), ios::app );
//           WriteAligns( out, v );    }    }
