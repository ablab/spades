///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "ShortVector.h"
#include "math/Functions.h"

#define NO_COLOR

// Note: If you want to pipe output of these routines to some more-like program,
// I suggest piping to "less -r".

// ================================================================================
//
// PrintBlanks
// PrintBases
//
// ================================================================================

void PrintBlanks( ostream& out, int n )
{    for ( int i = 0; i < n; i++ )
          out << " ";    }

/*
template<class BASEVEC>
void PrintBases( ostream& out, int blanks, const BASEVEC& rd, int from, int to )
{    PrintBlanks(out, blanks);
     for ( int i = from; i < to; i++ )
          out << as_base(rd[i]);
     out << "\n";    }
*/

void PrintBases( ostream& out, int blanks, const basevector& rd, int from, int to )
{    PrintBlanks(out, blanks);
     for ( int i = from; i < to; i++ )
          out << as_base(rd[i]);
     out << "\n";    }

void PrintBases( ostream& out, int blanks, const fastavector& rd, int from, int to )
{    PrintBlanks(out, blanks);
     for ( int i = from; i < to; i++ )
          out << rd[i];
     out << "\n";    }

void PrintBases( ostream& out, int blanks, const vec<char>& rd, int from, int to )
{    PrintBlanks(out, blanks);
     for ( int i = from; i < to; i++ )
          out << rd[i];
     out << "\n";    }

void PrintScores1( ostream& out, int blanks, const qualvector& s, int from, int to )
{    PrintBlanks(out, blanks);
     for ( int i = from; i < to; i++ )
     {    if ( s[i] >= 100 ) out << 'H';
          else if ( s[i] >= 10 ) out << (unsigned char)(s[i]/10 + '0');
          else out << (unsigned char)(s[i] + '0');    }
     out << "\n";    }

void PrintScores2( ostream& out, int blanks, const qualvector& s, int from, int to )
{    PrintBlanks(out, blanks);
     for ( int i = from; i < to; i++ )
     {    if ( s[i] >= 100 ) out << 'H';
          else if ( s[i] >= 10 ) out << (unsigned char)(s[i]%10 + '0');
          else out << ' ';    }
     out << "\n";    }

char AsBase( const basevector& b, const int j )
{    return as_base( b[j] );    }

char AsBase( const fastavector& b, const int j )
{    return b[j];    }

char AsBase( const vec<char>& b, const int j )
{    return b[j];    }

Bool Equal( const basevector& b1, const int j1, const basevector& b2, const int j2 )
{    return b1[j1] == b2[j2];    }

Bool Equal( const vec<char>& b1, const int j1, const basevector& b2, const int j2 )
{    return toupper(b1[j1]) == as_base(b2[j2]);    }

Bool Equal( const fastavector& b1, const int j1, const basevector& b2, const int j2 )
{    return GeneralizedBase::fromChar(b1[j1])
          .matches( GeneralizedBase::fromChar( as_base( b2[j2] ) ) );    }

Bool Equal( const fastavector& b1, const int j1, 
     const fastavector& b2, const int j2 )
{    return GeneralizedBase::fromChar(b1[j1])
          .matches( GeneralizedBase::fromChar(b2[j2]) );    }

// ================================================================================
//
// PrintVisualAlignment: given two reads (rd1, rd2), and an alignment a of them,
// pretty print the overlapping region.
// If scores are provided, they are assumed to correspond to the first sequence.
//
// If abbreviate = True, abbreviate perfectly aligning regions whose rd1 quality
// scores all have score at least min_score_to_abbrev.  Also, scores below
// min_score_to_abbrev are flagged with a '-'.
//
// If abbreviate_poor = True, abbreviate poorly aligning regions, defined as those
// having mismatch rate >= min_fract_poor.
//
// ================================================================================

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignment( Bool abbreviate, ostream& out, const BASEVEC1& rd1,
     const BASEVEC2& rd2, const align& a, const qualvector& scores1,
     const qualvector& scores2, int begin, Bool one_frame, int min_score_to_abbrev,
     Bool abbreviate_poor, float min_fract_poor,
     Bool abbreviate_good, float max_fract_good, Bool print_heads_and_tails,
     const Bool CtoT_special )

{    int pw = 80; // output width
     const avector<int>& gaps = a.Gaps( );
     const avector<int>& lengths = a.Lengths( );
     int nblocks = a.Nblocks( );
     int pos1 = a.pos1( ), pos2 = a.pos2( );
     ForceAssertGe( a.pos1( ), 0 );
     ForceAssertGe( a.pos2( ), 0 );
     ForceAssertLt( a.pos1( ), (int) rd1.size( ) );
     ForceAssertLt( a.pos2( ), (int) rd2.size( ) );
     ForceAssertLe( a.Pos1( ), (int) rd1.size( ) );
     ForceAssertLe( a.Pos2( ), (int) rd2.size( ) );

     // Print some nucleotides (if any) before the match.

     if ( print_heads_and_tails && pos1 > 0 && pos2 > 0 )
     {    out << "\nHEADS:\n";
       //          for ( int l = 0; l < pw; l++ )
       //               out << '-';
          out << "\n";
          if ( pos1 <= pw && pos2 <= pw )
          {    if ( scores1.size( ) > 0 )
               {    PrintScores1( out, Max( pos2 - pos1, 0 ), scores1, 0, pos1 );
                    PrintScores2( out, Max( pos2 - pos1, 0 ), scores1, 0, pos1 );
                    out << "\n";    }
               PrintBases( out, Max( pos2 - pos1, 0 ), rd1, 0, pos1 );
               PrintBases( out, Max( pos1 - pos2, 0 ), rd2, 0, pos2 );
               if ( scores2.size( ) > 0 )
               {    out << "\n";
                    PrintScores1( out, Max( pos1 - pos2, 0 ), scores2, 0, pos2 );
                    PrintScores2( out, Max( pos1 - pos2, 0 ), scores2, 0, pos2 ); } }
          else if ( pos1 >= pw - 3 )
          {    if ( scores1.size( ) > 0 )
               {    out << "...";
                    PrintScores1( out, 0, scores1, pos1-(pw-3), pos1 );
                    out << "...";
                    PrintScores2( out, 0, scores1, pos1-(pw-3), pos1 );
                    out << "\n";    }
               out << "...";
               PrintBases( out, 0, rd1, pos1-(pw-3), pos1 );
               if ( pos2 >= pw - 3 )
               {    out << "...";
                    PrintBases( out, 0, rd2, pos2-(pw-3), pos2 );
                    if ( scores2.size( ) > 0 )
                    {    out << "\n...";
                         PrintScores1( out, 0, scores2, pos2-(pw-3), pos2 );
                         out << "...";
                         PrintScores2( out, 0, scores2, pos2-(pw-3), pos2 );   }   }
               else
               {    PrintBases( out, pw - pos2, rd2, 0, pos2 );
                    if ( scores2.size( ) > 0 )
                    {    out << "\n";
                         PrintScores1( out, pw - pos2, scores2, 0, pos2 );
                         PrintScores2( out, pw - pos2, scores2, 0, pos2 );  }  }  }
          else
          {    if ( scores1.size( ) > 0 )
               {    PrintScores1( out, pw - pos1, scores1, 0, pos1 );
                    PrintScores2( out, pw - pos1, scores1, 0, pos1 );
                    out << "\n";    }
               PrintBases( out, pw - pos1, rd1, 0, pos1 );
               out << "...";
               PrintBases( out, 0, rd2, pos2-(pw-3), pos2 );
               if ( scores2.size( ) > 0 )
               {    out << "\n";
                    PrintScores1( out, 0, scores2, pos2-(pw-3), pos2 );
                    PrintScores2( out, 0, scores2, pos2-(pw-3), pos2 );    }    }
	  //          for ( int l = 0; l < pw; l++ )
	  //               out << '-';
          out << "\n";    }

     // Print the match itself.

     int enuf = 0;
     for ( int j = 0; j < nblocks; j++ )
     {    if ( gaps(j) > 0 ) enuf += gaps(j);
          if ( gaps(j) < 0 ) enuf -= gaps(j);
          enuf += lengths(j);    }
     vector<char> v1(enuf), v2(enuf), e(enuf), s1(enuf), s2(enuf), t1(enuf),
          t2(enuf);

     int i = 0, p1 = pos1, p2 = pos2;
     for ( int j = 0; j < nblocks; j++ )
     {    if ( gaps(j) > 0 )
               for ( int x = 0; x < gaps(j); x++ )
               {    if ( p1 >= begin )
                    {    v1[i] = ' ';
                         v2[i] = AsBase( rd2, p2 );
                         e[i] = '|';
                         if ( scores1.size( ) > 0 )
                         {    s1[i] = ' ';
                              s2[i] = ' ';    }
                         if ( scores2.size( ) > 0 )
                         {    if ( scores2[p2] >= 100 ) t1[i] = t2[i] = 'H';
                              else if ( scores2[p2] >= 10 )
                              {    t1[i] = scores2[p2]/10 + '0';
                                   t2[i] = scores2[p2]%10 + '0';    }
                              else
                              {    t1[i] = scores2[p2] + '0';
                                   t2[i] = ' ';    }    }
                         ++i;    }
                    p2++;    }
          if ( gaps(j) < 0 )
               for ( int x = 0; x < -gaps(j); x++ )
               {    if ( p1 >= begin )
                    {    v2[i] = ' ';
                         v1[i] = AsBase( rd1,p1 );
                         if ( scores1.size( ) > 0 )
                         {    if ( scores1[p1] >= 100 ) s1[i] = s2[i] = 'H';
                              else if ( scores1[p1] >= 10 )
                              {    s1[i] = scores1[p1]/10 + '0';
                                   s2[i] = scores1[p1]%10 + '0';    }
                              else
                              {    s1[i] = scores1[p1] + '0';
                                   s2[i] = ' ';    }    }
                         if ( scores2.size( ) > 0 )
                         {    t1[i] = ' ';
                              t2[i] = ' ';    }
                         e[i++] = '|';    }
                    ++p1;    }
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( p1 >= begin )
               {    v1[i] = AsBase( rd1, p1 );
                    v2[i] = AsBase( rd2, p2 );
                    if ( scores1.size( ) > 0 )
                    {    if ( scores1[p1] >= 100 ) s1[i] = s2[i] = 'H';
                         else if ( scores1[p1] >= 10 )
                         {    s1[i] = scores1[p1]/10 + '0';
                              s2[i] = scores1[p1]%10 + '0';    }
                         else
                         {    s1[i] = scores1[p1] + '0';
                              s2[i] = ' ';    }    }
                    if ( scores2.size( ) > 0 )
                    {    if ( scores2[p2] >= 100 ) t1[i] = t2[i] = 'H';
                         else if ( scores2[p2] >= 10 )
                         {    t1[i] = scores2[p2]/10 + '0';
                              t2[i] = scores2[p2]%10 + '0';    }
                         else
                         {    t1[i] = scores2[p2] + '0';
                              t2[i] = ' ';    }    }
                    if ( !Equal( rd1, p1, rd2, p2 ) )
                    {    if ( !CtoT_special || rd1[p1] != 3 || rd2[p2] != 1 )
                              e[i] = '*';
                         else e[i] = '.';    }
                    else if ( CtoT_special && rd1[p1] == 1 && rd2[p2] == 1 )
                         e[i] = 'm';
                    else if ( scores1.size( ) > 0
                         && scores1[p1] < min_score_to_abbrev ) e[i] = '-';
                    else e[i] = ' ';
                    ++i;    }
               ++p1; ++p2;    }    }
     int rspace = MIN( int(rd1.size( )) - p1, int(rd2.size( )) - p2 );
     if ( abbreviate && !( (pos1 > 0 && pos2 > 0) || (rspace > 0) ) )
     {    int min_score = 1000;
          for ( int u = 0; u < (int) scores1.size( ); u++ )
               min_score = Min( min_score, (int) scores1[u] );
          if ( min_score >= min_score_to_abbrev && ActualErrors( rd1, rd2, a ) == 0 )
          {    out << "(perfect match of length " << i << ")\n\n";
               return;    }    }
     if ( print_heads_and_tails
          && ( (pos1 > 0 && pos2 > 0) || (rspace > 0) ) )
     {    out << "MATCH:\n";    }
     else out << "\n";
     int j = 0;
     if ( scores2.size( ) > 0 )
     {
          #ifndef NO_COLOR
               out << START_MAGENTA;
          #endif
          for ( int l = 0; l < pw; l++ )
               out << '-';
          #ifndef NO_COLOR
               out << END_ESCAPE;
          #endif
	  out << "\n";
     }
     Bool first_frame = True;
     Bool newline_last = False;
     while( pw*j < i )
     {
          if ( !first_frame && one_frame ) break;
          first_frame = False;

          if ( abbreviate && j > 0 )
          {    int matching_bases = 0;
               int k = 0;
               while( pw*j < i )
               {    for ( k = pw*j; k < pw*(j+1) && k < i; k++ )
                         if ( e[k] != ' ' ) break;
                    if ( k == pw*(j+1) || k == i )
                    {    matching_bases += k - pw*j;
                         ++j;    }
                    else break;    }
               if ( matching_bases > 0 )
               {    out << "(" << matching_bases << " matching bases)\n\n";
                    newline_last = True;    }
               if ( k == i ) break;    }

          if ( abbreviate_poor && abbreviate_good )
          {    FatalErr( "Currently only one of abbreviate_poor and abbreviate_good "
                    "may be used at a time." );    }

          if (abbreviate_poor)
          {    int mismatching_bases = 0;
               int k = 0;
               while( pw*j < i )
               {    int errors = 0, total_bases = 0;
                    for ( k = pw*j; k < pw*(j+1) && k < i; k++ )
                    {    ++total_bases;
                         if ( e[k] == '*' ) ++errors;    }
                    if ( float(errors)/float(total_bases) >= min_fract_poor )
                    {    mismatching_bases += k - pw*j;
                         ++j;    }
                    else break;    }
               if ( mismatching_bases > 0 )
               {    out << "(" << mismatching_bases << " poorly matching bases)\n\n";
                    newline_last = True;    }
               if ( k == i ) break;    }

          if (abbreviate_good)
          {    int mismatching_bases = 0;
               int k = 0;
               while( pw*j < i )
               {    int errors = 0, total_bases = 0;
                    for ( k = pw*j; k < pw*(j+1) && k < i; k++ )
                    {    ++total_bases;
                         if ( e[k] == '*' ) ++errors;    }
                    if ( float(errors)/float(total_bases) <= max_fract_good )
                    {    mismatching_bases += k - pw*j;
                         ++j;    }
                    else break;    }
               if ( mismatching_bases > 0 )
               {    out << "(" << mismatching_bases << " well matching bases)\n\n";
                    newline_last = True;    }
               if ( k == i ) break;    }

          if ( scores1.size( ) > 0 )
          {    for ( int k = pw*j; k < pw*(j+1) && k < i; k++ )
                    out << s1[k];
               out << "\n";
               for ( int k = pw*j; k < pw*(j+1) && k < i; k++ )
                    out << s2[k];
               out << "\n\n";    }
          for ( int k = pw*j; k < pw*(j+1) && k < i; k++ )
               out << e[k];
          out << "\n";
          for ( int k = pw*j; k < pw*(j+1) && k < i; k++ )
               out << v1[k];
          out << "\n";
          for ( int k = pw*j; k < pw*(j+1) && k < i; k++ )
               out << v2[k];
          out << "\n\n";
          if ( scores2.size( ) > 0 )
          {    for ( int k = pw*j; k < pw*(j+1) && k < i; k++ )
                    out << t1[k];
               out << "\n";
               for ( int k = pw*j; k < pw*(j+1) && k < i; k++ )
                    out << t2[k];
               out << "\n";    }
          newline_last = True;
          if ( scores2.size( ) > 0 )
          {
               #ifndef NO_COLOR
                    out << START_MAGENTA;
               #endif
               for ( int l = 0; l < pw; l++ )
                    out << '-';
               newline_last = False;
               #ifndef NO_COLOR
                    out << END_ESCAPE;
               #endif
                    }
          if ( pw * (j+1) < i ) 
          {    out << "\n";
               newline_last = True;    }
          ++j;    }
     if ( !newline_last ) out << "\n";

     // Print some nucleotides (if any) after the match.

     char tails1[85], tails2[85], s1x[85], s2x[85];
     if ( print_heads_and_tails && rspace > 0 )
     {    int tails1_ptr = 0, tails2_ptr = 0;

          if ( scores1.size( ) > 0 )
          {    i = 0;
               for ( j = p1; j < Min( p1+80, int(rd1.size( )), p1+rspace+5 ); j++ )
               {    if ( scores1[j] >= 100 ) s1x[i] = s2x[i] = 'H';
                    else if ( scores1[j] >= 10 )
                    {    s1x[i] = scores1[j]/10 + '0';
                         s2x[i] = scores1[j]%10 + '0';    }
                    else
                    {    s1x[i] = scores1[j] + '0';
                         s2x[i] = ' ';    }
                    ++i;    }
               if ( j != int(rd1.size() ) )
		    for (int tt = 0; tt < 3; tt++ )
		    {    s1x[i] = '.';
			 s2x[i] = '.';
			 ++i;    }    }

          for ( j = p1; j < Min( p1+80, int(rd1.size( )), p1+rspace+5 ); j++ )
               tails1[tails1_ptr++] = AsBase(rd1,j);
          if ( j != int(rd1.size( )) )
               for ( int tt = 0; tt < 3; tt++ )
                    tails1[tails1_ptr++] = '.';
          for ( j = p2; j < Min( p2+80, int(rd2.size( )), p2+rspace+5 ); j++ )
               tails2[tails2_ptr++] = AsBase(rd2,j);
          if ( j != int(rd2.size( )) )
               for ( int tt = 0; tt < 3; tt++ )
                    tails2[tails2_ptr++] = '.';
          out << "TAILS:\n\n";
          if ( scores1.size( ) > 0 )
          {    for ( j = 0; j < Min( 80, tails1_ptr ); j++ )
                    out << s1x[j];
               out << "\n";
               for ( j = 0; j < Min( 80, tails1_ptr ); j++ )
                    out << s2x[j];
               out << "\n\n";    }
          for ( j = 0; j < Min( 80, tails1_ptr ); j++ )
               out << tails1[j];
          out << "\n";
          for ( j = 0; j < Min( 80, tails2_ptr ); j++ )
               out << tails2[j];
          out << "\n\n";    }    }

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignment( Bool rd2_is_rc, Bool abbreviate, ostream& out,
     const BASEVEC1& rd1, BASEVEC2 rd2, const align& a,
     const qualvector& scores1, qualvector scores2, int begin, Bool one_frame,
     int min_score_to_abbrev, Bool abbreviate_poor, float min_fract_poor,
     Bool abbreviate_good, float max_fract_good )
{    if (rd2_is_rc)
     {    rd2.ReverseComplement( );
          scores2 = Reverse(scores2);    }
     PrintVisualAlignment( abbreviate, out, rd1, rd2, a, scores1, scores2,
          begin, one_frame, min_score_to_abbrev, abbreviate_poor,
          min_fract_poor );    }

template
void PrintVisualAlignment( Bool abbreviate, ostream& out, const basevector& rd1,
     const basevector& rd2, const align& a, const qualvector& scores1,
     const qualvector& scores2, int begin, Bool one_frame, int min_score_to_abbrev,
     Bool abbreviate_poor, float min_fract_poor,
     Bool abbreviate_good, float max_fract_good, Bool print_heads_and_tails,
     const Bool CtoT_special );

template
void PrintVisualAlignment( Bool abbreviate, ostream& out, const vec<char>& rd1,
     const basevector& rd2, const align& a, const qualvector& scores1,
     const qualvector& scores2, int begin, Bool one_frame, int min_score_to_abbrev,
     Bool abbreviate_poor, float min_fract_poor,
     Bool abbreviate_good, float max_fract_good, Bool print_heads_and_tails,
     const Bool CtoT_special );

template
void PrintVisualAlignment( Bool abbreviate, ostream& out, const fastavector& rd1,
     const basevector& rd2, const align& a, const qualvector& scores1,
     const qualvector& scores2, int begin, Bool one_frame, int min_score_to_abbrev,
     Bool abbreviate_poor, float min_fract_poor,
     Bool abbreviate_good, float max_fract_good, Bool print_heads_and_tails,
     const Bool CtoT_special );

template
void PrintVisualAlignment( Bool abbreviate, ostream& out, const fastavector& rd1,
     const fastavector& rd2, const align& a, const qualvector& scores1,
     const qualvector& scores2, int begin, Bool one_frame, int min_score_to_abbrev,
     Bool abbreviate_poor, float min_fract_poor,
     Bool abbreviate_good, float max_fract_good, Bool print_heads_and_tails,
     const Bool CtoT_special );

template void PrintVisualAlignment( Bool rd2_is_rc, Bool abbreviate, ostream& out, 
     const basevector& rd1, basevector rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), qualvector scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05 );

template void PrintVisualAlignment( Bool rd2_is_rc, Bool abbreviate, ostream& out, 
     const fastavector& rd1, fastavector rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), qualvector scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05 );

// ================================================================================
//
// PrintAlignment: print the data of an alignment "a"
//
// ================================================================================

void PrintAlignment(ostream& out, const basevector& rd1, const basevector& rd2,
     const alignment& a)
{    int pos1, pos2, errors;
     avector<int> gaps, lengths;
     a.Unpack( pos1, pos2, errors, gaps, lengths );
     out << "   Score (bigger=badder): " << errors << ".\n";
     out << "   The overlap on read 1 starts at position " << pos1 << ".\n";
     out << "   The overlap on read 2 starts at position " << pos2 << ".\n";
     for ( unsigned int j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) > 0 )
               out << "   Now there's a gap of " << gaps(j) << " on read 1.\n";
          else if ( gaps(j) < 0 )
               out << "   Now there's a gap of " << -gaps(j) << " on read 2.\n";
          out << "   Skip " << lengths(j) << " to the right on both reads.\n";    }

     int p1 = pos1, p2 = pos2, total_len = 0;
     for ( unsigned int j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) > 0 ) p2 += gaps(j);
          else if ( gaps(j) < 0 ) p1 -= gaps(j);
          p1 += lengths(j); p2 += lengths(j);
          total_len += lengths(j);    }
     out << "   [Pos1 = " << p1 << ", Pos2 = " << p2 << ", "
          << "rd1.size( ) = " << rd1.size( )
          << ", rd2.size( ) = " << rd2.size( ) << "]\n";

     out << "\n";    }

void PrintAlignment( Bool rd2_is_rc, ostream& out, const basevector& rd1,
     basevector rd2, const alignment& a)
{    if (rd2_is_rc) rd2.ReverseComplement( );
     PrintAlignment( out, rd1, rd2, a );    }

void PrintReadWithScores( basevector& B, qualvector& Q, ostream& out, const int pw )
{    Assert( B.size( ) == Q.size( ) );
     int n = B.size( );
     vector<char> b(2*n), s1(2*n), s2(2*n);
     for ( int j = 0; j < n; j++ )
     {    b[j] = as_base( B[j] );
          if ( Q[j] >= 10 )
          {    s1[j] = Q[j]/10 + '0';
               s2[j] = Q[j]%10 + '0';    }
          else
          {    s1[j] = Q[j] + '0';
               s2[j] = ' ';    }    }
     int j = 0;
     while( pw*j < n )
     {    for ( int K = pw*j; K < pw*(j+1) && K < n; K++ )
               out << b[K];
          out << "\n\n";
          for ( int K = pw*j; K < pw*(j+1) && K < n; K++ )
               out << s1[K];
          out << "\n";
          for ( int K = pw*j; K < pw*(j+1) && K < n; K++ )
               out << s2[K];
          out << "\n\n";
          ++j;    }    }

void PrintErrorsInAlignment(ostream& out, const basevector& rd1,
     const basevector& rd2, const alignment& a, const qualvector& scores1,
     const qualvector& scores2 )
{    avector<int> gaps, lengths;
     int pos1, pos2, errors;
     a.Unpack( pos1, pos2, errors, gaps, lengths );
     int p1 = pos1, p2 = pos2;
     out << "\nerror summary:\n\n";
     for ( unsigned int j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) > 0 )
               for ( int x = 0; x < gaps(j); x++ )
               {    out << "   (" << p1 << "," << p2 << ") = (gap,"
			<< as_base(rd2[p2]) << ") @ [ ," << int(scores2[p2])
                         << "]\n";
                    ++p2;    }
          if ( gaps(j) < 0 )
               for ( int x = 0; x < -gaps(j); x++ )
               {    out << "   (" << p1 << "," << p2 << ") = (" << as_base(rd1[p1])
                         << ",gap) @ (" << int(scores1[p1])
                         << ", )\n";
                    ++p1;    }
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( rd1[p1] != rd2[p2] )
                    out << "   (" << p1 << "," << p2 << ") = (" << as_base(rd1[p1])
                         << "," << as_base(rd2[p2]) << ") @ ["
                         << int(scores1[p1]) << "," << int(scores2[p2]) << ")\n";
               ++p1;
               ++p2;     }    }    }
