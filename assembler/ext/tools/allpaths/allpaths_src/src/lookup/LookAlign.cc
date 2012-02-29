///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <set>
#include <sstream>
#include <strstream>

#include "Alignment.h"
#include "Basevector.h"
#include "BlockAlign.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "ReadPairing.h"
#include "ScoreAlignment.h"
#include "ShortVector.h"
#include "STLExtensions.h"
#include "VecTemplate.h"
#include "lookup/LookAlign.h"
#include "math/Arith.h"
#include "math/Functions.h"
#include "random/Random.h"
#include "system/ParsedArgs.h"

String QUERY("QUERY");

void look_align::ResetFromAlign(const align & al, const basevector & b1,
				const basevector & b2) {
  a = al;
  vec<int> errs = a.MutationsGap1Gap2(b1, b2);
  mutations = errs[0];
  indels = errs[1] + errs[2];
}


look_align look_align::TrimmedTo1(int startOn1, int len,
				    const basevector & b1,
				    const basevector & b2) const {
    look_align ret(*this);
    ret.a = a.TrimmedTo1(startOn1, len);
    ret.query_length = len;
    basevector btrimmed;
    btrimmed.SetToSubOf(b1, startOn1, len);
    vec<int> errs = ret.a.MutationsGap1Gap2(btrimmed, b2);
    ret.mutations = errs[0];
    ret.indels = errs[1] + errs[2];
    return ret;
}

look_align_plus look_align::TrimmedTo1Plus(int startOn1, int len,
				    const basevector & b1,
				    const basevector & b2) const {
  look_align temp= TrimmedTo1(startOn1, len, b1, b2);
  ostringstream os;
  basevector btrimmed;
  btrimmed.SetToSubOf(b1, startOn1, len);
  temp.PrintParseable(os, &btrimmed, &b2);
  look_align_plus ret;
  ret.ReadParseable(os.str().c_str());
  return ret;
}

void look_align::PrintParseable( ostream& out, const basevector& query,
                                 const qualvector& query_qual, const basevector& target,
                                 int start_on_target, const vec<String>& seq_names,
                                 const vec<String>& target_names ) const
{
  out << QUERY << "\t" << seq_names[query_id] << "\t" << a.pos1( ) << "\t"
      << a.Pos1( ) << "\t" << query_length << "\t" << int(rc1) << "\t"
      << target_names[target_id] << "\t" << a.pos2( ) << "\t" << a.Pos2( )
      << "\t" << target_length << "\t" << a.Nblocks( );
  int p1 = a.pos1( ), p2 = a.pos2( );
  for ( int j = 0; j < a.Nblocks( ); j++ ) {
    int mismatches = 0;
    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
    for ( int x = 0; x < a.Lengths(j); x++ ) {
      if ( !rc1 && query[p1] != target[ p2 - start_on_target ] )
        ++mismatches;
      if ( rc1 && 3 - query[query_length - p1 - 1]
           != target[ p2 - start_on_target ] ) {
        ++mismatches;
      }
      ++p1, ++p2;
    }
    out << "\t" << a.Gaps(j) << "\t" << a.Lengths(j) << "\t"
        << mismatches;
  }
  out << "\n";
}

void look_align::PrintParseable( ostream& out,
				 const basevector *query,
				 const basevector *target, bool endl  ) const
{
  // Warning: if query or target are null, then the sum of all mutations
  // is saved in the first block. This is not ideal, but is better than
  // losing the information altogether.

  out << "QUERY\t"
      << query_id << "\t"
      << a.pos1( ) << "\t"
      << a.Pos1( ) << "\t"
      << query_length << "\t"
      << int(rc1) << "\t"
      << target_id << "\t"
      << a.pos2( ) << "\t"
      << a.Pos2( ) << "\t"
      << target_length << "\t"
      << a.Nblocks( );

  if ( query && target ) {
    int p1 = a.pos1( );
    int p2 = a.pos2( );
    int start_on_target = 0;
    for (int j=0; j<a.Nblocks( ); j++) {
      int mismatches = 0;
      if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
      if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
      for ( int x = 0; x < a.Lengths(j); x++ ) {
	if ( !rc1 && (*query)[p1] != (*target)[ p2 - start_on_target ] )
	  ++mismatches;

	if ( rc1 && (query_length > static_cast<unsigned int>(p1)) && 3 - (*query)[query_length - p1 - 1]
	     != (*target)[ p2 - start_on_target ] ) {
	  ++mismatches;
	}
	++p1, ++p2;
      }
      out << "\t" << a.Gaps(j)
	  << "\t" << a.Lengths(j)
	  << "\t"  << mismatches;
    }
  }
  else {
    for (int ii=0; ii<a.Nblocks( ); ii++) {
      out << "\t" << a.Gaps(ii)
	  << "\t" << a.Lengths(ii)
	  << "\t" << (0 == ii ? mutations : 0);//save on first block.
    }
  }

  if ( endl ) out << "\n";
}

Bool look_align::ReadParseableOrFail( const String& in )
{
  istrstream ins( in.c_str( ) );
  String field, target_id_s;
  ins >> field;
  ForceAssert( field == QUERY );
  int pos1, pos2, irc1, nblocks;
  ins >> query_id >> pos1 >> field >> query_length >> irc1 >> target_id_s;
  if ( !target_id_s.IsInt( ) ) return False;
  target_id = target_id_s.Int( );
  ins >> pos2 >> field >> target_length >> nblocks;
  rc1 = irc1;
  avector<int> g, l;
  g.Setsize(nblocks), l.Setsize(nblocks);
  mutations = 0, indels = 0;
  for ( int i = 0; i < nblocks; i++ ) {
    int mut;
    ins >> g(i) >> l(i) >> mut;
    if ( l(i) <= 0 ) return False;
    mutations += mut;
    indels += Abs( g(i) );
  }
  a.Set( pos1, pos2, g, l );
  return True;
}

void look_align::ReadParseable( const String& in )
{
  istrstream ins( in.c_str( ) );
  ReadParseable(ins);
}

void look_align::ReadParseable( istream & ins )
{
  String field, target_id_s;
  ins >> field;
  ForceAssert( field == QUERY );
  int pos1, pos2, irc1, nblocks;
  ins >> query_id >> pos1 >> field >> query_length >> irc1 >> target_id_s;
  if ( !target_id_s.IsInt( ) )
  {
    FatalErr( "Alignment of query " << query_id << " to target "
              << target_id_s << ": can't read non-numeric target name." );
  }
  target_id = target_id_s.Int( );
  ins >> pos2 >> field >> target_length >> nblocks;
  rc1 = irc1;
  avector<int> g, l;
  g.Setsize(nblocks), l.Setsize(nblocks);
  mutations = 0, indels = 0;
  for ( int i = 0; i < nblocks; i++ ) {
    int mut;
    ins >> g(i) >> l(i) >> mut;
    if ( l(i) <= 0 ) {
      FatalErr( "Alignment of query " << query_id << " to target "
                << target_id << ": block " << i << " of " << nblocks
                << " has size " << l(i) << "." );
    }
    mutations += mut;
    indels += Abs( g(i) );
  }
  a.Set( pos1, pos2, g, l );
}

void look_align_plus::ReadParseable( const String& in )
{
  istrstream ins( in.c_str( ) );
  String field, target_id_s;
  ins >> field;
  ForceAssert( field == QUERY );
  int pos1, pos2, irc1, nblocks;
  ins >> query_id >> pos1 >> field >> query_length >> irc1 >> target_id_s;
  if ( !target_id_s.IsInt( ) ) {
    FatalErr( "Alignment of query " << query_id << " to target "
              << target_id_s << ": can't read non-numeric target name." );
  }
  target_id = target_id_s.Int( );
  ins >> pos2 >> field >> target_length >> nblocks;
  rc1 = irc1;
  avector<int> g, l;
  g.Setsize(nblocks), l.Setsize(nblocks);
  mutations = 0, indels = 0;
  mutations_by_block.clear( );
  for ( int i = 0; i < nblocks; i++ ) {
    int mut;
    ins >> g(i) >> l(i) >> mut;
    if ( l(i) < 0 || (l(i) == 0 && g(i) >= 0 && i != (nblocks - 1)) ) {
      FatalErr( "Alignment of query " << query_id << " to target "
                << target_id << ": block " << i << " of " << nblocks
                << " has size " << l(i) << "." );
    }
    mutations_by_block.push_back(mut);
    mutations += mut;
    indels += Abs( g(i) );
  }
  a.Set( pos1, pos2, g, l );
}

void look_align::PrintParseableBrief( ostream& out, const basevector& query,
                                      const qualvector& query_qual, const basevector& target,
                                      int start_on_target, const vec<String>& seq_names,
                                      const vec<String>& target_names ) const
{
  out << QUERY << "\t" << seq_names[query_id] << "\t" << a.pos1( ) << "\t"
      << a.Pos1( ) << "\t" << query_length << "\t" << int(rc1) << "\t"
      << target_names[target_id] << "\t" << a.pos2( ) << "\t"
      << a.Pos2( ) << "\t" << target_length << "\t" << a.Nblocks( )
      << "\n";
}

void look_align::PrintReadableBrief( ostream& out, const String & query_name,
				     const String & target_name ) const {
  out << query_name << ( rc1 ? "rc" : "fw" ) << " vs "
      << target_name << ", " << mutations << " mismatches/"
      << indels << " indels (of " << query_length << "), from "
      << a.pos1( ) << "-" << a.Pos1( ) << " to " << a.pos2( ) << "-"
      << a.Pos2( ) << " (of " << target_length << ")\n";
}

void look_align::PrintReadableBrief( ostream& out, const basevector& query,
                                     const qualvector& query_qual,
				     const basevector& target,
                                     int start_on_target,
				     const vec<String>& seq_names,
                                     const vec<String>& target_names ) const
{
  PrintReadableBrief(out, seq_names[query_id], target_names[target_id]);
}

void look_align::PrintReadableBrief( ostream& out) const
{
  PrintReadableBrief(out, ToString(query_id), ToString(target_id));
}

int look_align::CountNqs( const basevector& query, const qualvector& query_qual,
                          const basevector& target, int start_on_target ) const
{
  const int NQS_floor1 = 30;
  const int NQS_floor2 = 25;
  const int NQS_radius = 5;
  const int NQS_extra = 1;
  int see = 0;
  int p1 = a.pos1( ), p2 = a.pos2( );
  for ( int j = 0; j < a.Nblocks( ); j++ ) {
    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
    for ( int x = 0; x < a.Lengths(j); x++ ) {
      if ( x >= NQS_radius && x < a.Lengths(j) - NQS_radius ) {
        int xp1 = ( !rc1 ? p1 : query_length - p1 - 1 );
        int xp2 = p2 - start_on_target;
        Bool mismatch = False;
        if ( !rc1 && query[xp1] != target[xp2] ) mismatch = True;
        if ( rc1 && 3 - query[xp1] != target[xp2] ) mismatch = True;
        if (mismatch) {
          if ( query_qual[xp1] >= NQS_floor1 ) {
            Bool range_ok = True;
            for ( int y = 1; y <= NQS_radius; y++ ) {
              if ( query_qual[xp1-y] < NQS_floor2 || query_qual[xp1+y] < NQS_floor2 ) {
                range_ok = False;
              }
            }
            if (range_ok) {
              int ndiff = 0;
              for ( int y = -NQS_radius; y <= NQS_radius; y++ ) {
                int yp1 = ( !rc1 ? (p1+y) : query_length - (p1+y) - 1 );
                int yp2 = (p2+y) - start_on_target;
                if ( y != 0 && !rc1 && query[yp1] != target[yp2] ) {
                  ++ndiff;
                }
                if ( y != 0 && rc1 && 3 - query[yp1] != target[yp2] ) {
                  ++ndiff;
                }
              }
              if ( ndiff <= NQS_extra )
                ++see;
            }
          }
        }
      }
      ++p1, ++p2;
    }
  }
  return see;
}

int look_align::CountNqsWithReturns( const basevector& query, const qualvector& query_qual,
				     const basevector& target, int start_on_target,
				     const vec<String>& seq_names,
				     const vec<String>& target_names,
				     int &qualSnpBase,
				     int &targetSnpBase )
{
  const int NQS_floor1 = 30;
  const int NQS_floor2 = 25;
  const int NQS_radius = 5;
  const int NQS_extra = 1;
  int see = 0;
  int p1 = a.pos1( ), p2 = a.pos2( );
  for ( int j = 0; j < a.Nblocks( ); j++ ) {
    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
    for ( int x = 0; x < a.Lengths(j); x++ ) {
      if ( x >= NQS_radius && x < a.Lengths(j) - NQS_radius ) {
        int xp1 = ( !rc1 ? p1 : query_length - p1 - 1 );
        int xp2 = p2 - start_on_target;
        Bool mismatch = False;
        if ( !rc1 && query[xp1] != target[xp2] ) mismatch = True;
        if ( rc1 && 3 - query[xp1] != target[xp2] ) mismatch = True;
        if (mismatch) {
          if ( query_qual[xp1] >= NQS_floor1 ) {
            Bool range_ok = True;
            for ( int y = 1; y <= NQS_radius; y++ ) {
              if ( query_qual[xp1-y] < NQS_floor2
                   || query_qual[xp1+y] < NQS_floor2 ) {
                range_ok = False;
              }
            }
            if (range_ok) {
              int ndiff = 0;
              for ( int y = -NQS_radius; y <= NQS_radius; y++ ) {
                int yp1 = ( !rc1 ? (p1+y) : query_length - (p1+y) - 1 );
                int yp2 = (p2+y) - start_on_target;
                if ( y != 0 && !rc1 && query[yp1] != target[yp2] ) {
                  ++ndiff;
                }
                if ( y != 0 && rc1 && 3 - query[yp1] != target[yp2] ) {
                  ++ndiff;
                }
              }
              if ( ndiff <= NQS_extra ) {
                qualSnpBase = xp1;
                targetSnpBase = xp2;
                ++see;
              }
            }
          }
        }
      }
      ++p1, ++p2;
    }
  }
  return see;
}

void look_align::PrintNqs( ostream& out, const basevector& query,
                           const qualvector& query_qual, const basevector& target,
                           int start_on_target, const vec<String>& seq_names,
                           const vec<String>& target_names ) const
{
  int see = CountNqs( query, query_qual, target, start_on_target );
  out << "see " << see << " NQS(30,25) differences, " << "rate = "
      << setprecision(4) << 100.0 * float(see) / float( a.Pos1( ) - a.pos1( ) )
      << "%\n";
}

Float look_align::QualScore( const basevector& query,
                             const qualvector& query_qual, const basevector& target,
                             int start_on_target, const qualvector& qual_g ) const
{
  qualvector target_q;
  if ( qual_g.size( ) > 0 )
    target_q.SetToSubOf( qual_g, start_on_target, target.size( ) );
  else target_q.resize(0);
  align al;
  al = a;
  al.Setpos2( a.pos2( ) - start_on_target );
  basevector Query;
  qualvector Query_qual;
  Query = query;
  Query_qual = query_qual;
  if (rc1) {
    Query.ReverseComplement( );
    Query_qual.ReverseMe( );
  }
  Float score = ScoreAlignment( al, Query, Query_qual, target, target_q );
  return score;
}

void look_align::PrintQualScore( ostream& out, const basevector& query,
                                 const qualvector& query_qual, const basevector& target,
                                 int start_on_target, const qualvector& qual_g ) const
{
  out << "quality score = " << QualScore( query, query_qual, target,
                                          start_on_target, qual_g ) << "\n";
}

Float look_align::RmrPercent( const basevector& query,
                              const qualvector& query_qual, const basevector& target,
                              int start_on_target ) const
{
  int p1 = a.pos1( ), p2 = a.pos2( );
  vec<int> perf;
  perf.clear( );
  for ( int j = 0; j < a.Nblocks( ); j++ )
  {
    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
    int run = 0;
    for ( int x = 0; x < a.Lengths(j); x++ ) {
      Bool mismatch = False;
      if ( !rc1 && query[p1] != target[ p2 - start_on_target ] ) {
        mismatch = True;
      }
      if ( rc1 && 3 - query[query_length - p1 - 1] != target[ p2 - start_on_target ] ) {
        mismatch = True;
      }
      if (mismatch) {
        perf.push_back( run + 1 );
        run = 0;
      }
      else ++run;
      ++p1, ++p2;
    }
    perf.push_back( run + 1 );
  }
  return Float(100) / WeightedMean(perf);
}

Float look_align::MMPercent( const basevector& query,
                             const qualvector& query_qual, const basevector& target,
                             int start_on_target ) const
{
  if ( query_length == 0 ) return Float(100);
  int p1 = a.pos1( ), p2 = a.pos2( );
  int matches = 0, mismatches = 0;
  for ( int j = 0; j < a.Nblocks( ); j++ ) {
    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
    for ( int x = 0; x < a.Lengths(j); x++ ) {
      Bool mismatch = False;
      if ( !rc1 && query[p1] != target[ p2 - start_on_target ] ) {
        mismatch = True;
      }
      if ( rc1 && 3 - query[query_length - p1 - 1] != target[ p2 - start_on_target ] ) {
        mismatch = True;
      }
      if (mismatch) ++mismatches;
      else ++matches;
      ++p1, ++p2;
    }
  }
  return Float(100) * Float( Max( 0, matches - mismatches ) ) / Float(query_length);
}

void look_align::PrintRmr( ostream& out, const basevector& query,
                           const qualvector& query_qual, const basevector& target,
                           int start_on_target ) const
{
  out << "rmr = " << setprecision(2)
      << RmrPercent( query, query_qual, target, start_on_target ) << "%\n";
}

void look_align::PrintMM( ostream& out, const basevector& query,
                          const qualvector& query_qual, const basevector& target,
                          int start_on_target ) const
{
  Float mm = MMPercent( query, query_qual, target, start_on_target );
  if ( mm == Float(100) ) out << "mm = 100%\n";
  else out << "mm = " << setprecision(2) << mm << "%\n";
}

void look_align::PrintRmrByBlock( ostream& out, const basevector& query,
                                  const qualvector& query_qual, const basevector& target,
                                  int start_on_target, const vec<String>& seq_names,
                                  const vec<String>& target_names ) const
{
  out << "rmr by block:\n";
  int p1 = a.pos1( ), p2 = a.pos2( );
  for ( int j = 0; j < a.Nblocks( ); j++ ) {
    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
    vec<int> perf;
    perf.clear( );
    int run = 0;
    for ( int x = 0; x < a.Lengths(j); x++ ) {
      Bool mismatch = False;
      if ( !rc1 && query[p1] != target[ p2 - start_on_target ] ) {
        mismatch = True;
      }
      if ( rc1 && 3 - query[query_length - p1 - 1] != target[ p2 - start_on_target ] ) {
        mismatch = True;
      }
      if (mismatch) {
        perf.push_back( run + 1 );
        run = 0;
      }
      else ++run;
      ++p1, ++p2;
    }
    perf.push_back( run + 1 );
    out << a.Lengths(j) << " bases, rmr = "
        << setprecision(2) << Float(100) / WeightedMean(perf)
        << "%\n";
  }
}

void look_align::PrintVisual( ostream& out, const basevector & query,
                              const basevector& target, const Bool abbr ) const
{    if (rc1)
     {    basevector b = query;
          b.ReverseComplement( );
          PrintVisualAlignment( abbr, out, b, target, a );    }
     else PrintVisualAlignment( abbr, out, query, target, a );    }

void look_align::PrintVisual( ostream& out, const fastavector & query,
                              const basevector& target, const Bool abbr ) const
{    if (rc1)
     {    fastavector b = query;
          b.ReverseComplement( );
          PrintVisualAlignment( abbr, out, b, target, a );    }
     else PrintVisualAlignment( abbr, out, query, target, a );    }

void look_align::PrintVisual( ostream& out, const basevector& query,
                              const qualvector& query_qual, const basevector& target,
                              int start_on_target, const Bool abbr ) const
{    basevector rd1;
     rd1 = query;
     if (rc1) rd1.ReverseComplement( );
     qualvector q1;
     q1 = query_qual;
     if (rc1) q1.ReverseMe( );
     align ar;
     ar = a;
     ar.Setpos2( a.pos2( ) - start_on_target );
     PrintVisualAlignment( abbr, out, rd1, target, ar, q1 );     }

void look_align::PrintVisual( ostream& out, const fastavector& query,
                              const qualvector& query_qual, const basevector& target,
                              int start_on_target, const Bool abbr ) const
{    fastavector rd1;
     rd1 = query;
     if (rc1) rd1.ReverseComplement( );
     qualvector q1;
     q1 = query_qual;
     if (rc1) q1.ReverseMe( );
     align ar;
     ar = a;
     ar.Setpos2( a.pos2( ) - start_on_target );
     PrintVisualAlignment( abbr, out, rd1, fastavector(target), ar, q1 );     }

void look_align::PrintVisual( ostream& out, const basevector& query,
     const qualvector& query_qual, const fastavector& target,
     int start_on_target, const Bool abbr ) const
{    basevector rd1;
     rd1 = query;
     if (rc1) rd1.ReverseComplement( );
     qualvector q1;
     q1 = query_qual;
     if (rc1) q1.ReverseMe( );
     align ar;
     ar = a;
     ar.Setpos2( a.pos2( ) - start_on_target );
     PrintVisualAlignment( abbr, out, fastavector(rd1), target, ar, q1 );     }

void look_align::PrintVisual( ostream& out, const fastavector& query,
     const qualvector& query_qual, const fastavector& target,
     int start_on_target, const Bool abbr ) const
{    fastavector rd1;
     rd1 = query;
     if (rc1) rd1.ReverseComplement( );
     qualvector q1;
     q1 = query_qual;
     if (rc1) q1.ReverseMe( );
     align ar;
     ar = a;
     ar.Setpos2( a.pos2( ) - start_on_target );
     PrintVisualAlignment( abbr, out, fastavector(rd1), target, ar, q1 );     }

void look_align::PrintVisual( ostream& out, const basevector& query,
                              const qualvector& query_qual, 
                              const qualvector& target_qual, 
                              const basevector& target,
                              int start_on_target, const Bool abbr ) const
{    basevector rd1;
     rd1 = query;
     if (rc1) rd1.ReverseComplement( );
     qualvector q1;
     q1 = query_qual;
     if (rc1) q1.ReverseMe( );
     align ar;
     ar = a;
     ar.Setpos2( a.pos2( ) - start_on_target );
     PrintVisualAlignment( abbr, out, rd1, target, ar, q1, target_qual );     }

void look_align::PrintVisual( ostream& out, const fastavector& query,
                              const qualvector& query_qual, 
                              const qualvector& target_qual, 
                              const basevector& target,
                              int start_on_target, const Bool abbr ) const
{    fastavector rd1;
     rd1 = query;
     if (rc1) rd1.ReverseComplement( );
     qualvector q1;
     q1 = query_qual;
     if (rc1) q1.ReverseMe( );
     align ar;
     ar = a;
     ar.Setpos2( a.pos2( ) - start_on_target );
     PrintVisualAlignment( abbr, out, rd1, fastavector(target), ar, q1, target_qual );     }

void look_align::PrintVisual( ostream& out, const basevector& query,
                              const qualvector& query_qual, 
                              const qualvector& target_qual, 
                              const fastavector& target,
                              int start_on_target, const Bool abbr ) const
{    basevector rd1;
     rd1 = query;
     if (rc1) rd1.ReverseComplement( );
     qualvector q1;
     q1 = query_qual;
     if (rc1) q1.ReverseMe( );
     align ar;
     ar = a;
     ar.Setpos2( a.pos2( ) - start_on_target );
     PrintVisualAlignment( abbr, out, fastavector(rd1), target, ar, q1, target_qual );     }

void look_align::PrintVisual( ostream& out, const fastavector& query,
                              const qualvector& query_qual, 
                              const qualvector& target_qual, 
                              const fastavector& target,
                              int start_on_target, const Bool abbr ) const
{    fastavector rd1;
     rd1 = query;
     if (rc1) rd1.ReverseComplement( );
     qualvector q1;
     q1 = query_qual;
     if (rc1) q1.ReverseMe( );
     align ar;
     ar = a;
     ar.Setpos2( a.pos2( ) - start_on_target );
     PrintVisualAlignment( abbr, out, rd1, target, ar, q1, target_qual );     }

void look_align::PrintVisual( ostream& out, const basevector & query,
                              const basevector& target, const basevector& target_rc,
                              const Bool abbr, const Bool reverse_display ) const
{    if (rc1)
     {    basevector b = query;
          if ( !reverse_display )
          {    b.ReverseComplement( );
               PrintVisualAlignment( abbr, out, b, target, a );    }
          else
          {    align arc = a;
               arc.ReverseThis( query.size( ), target.size( ) );
               PrintVisualAlignment( abbr, out, query, target_rc, arc );    }    }
     else PrintVisualAlignment( abbr, out, query, target, a );    }

void look_align::PrintVisual( ostream& out, const basevector& query,
                              const qualvector& query_qual, const basevector& target,
                              const basevector& target_rc, int start_on_target,
                              const Bool abbr, const Bool reverse_display ) const
{    basevector rd1;
     qualvector q1;
     rd1 = query, q1 = query_qual;
     align ar;
     ar = a;
     if ( !reverse_display || !rc1 )
     {    if (rc1) rd1.ReverseComplement( );
          if (rc1) q1.ReverseMe( );
          ar.Setpos2( a.pos2( ) - start_on_target );
          PrintVisualAlignment( abbr, out, rd1, target, ar, q1 );    }
     else
     {    ar.Setpos2( a.pos2( ) - start_on_target );
          ar.ReverseThis( query.size( ), target.size( ) );
          PrintVisualAlignment( abbr, out, rd1, target_rc, ar, q1 );    }    }

void look_align_plus::WriteParseable( ostream& out ) const
{
  out << QUERY << "\t" << query_id << "\t" << a.pos1( ) << "\t"
      << a.Pos1( ) << "\t" << query_length << "\t" << int(rc1) << "\t"
      << target_id << "\t" << a.pos2( ) << "\t" << a.Pos2( )
      << "\t" << target_length << "\t" << a.Nblocks( );
  int p1 = a.pos1( ), p2 = a.pos2( );
  for ( int j = 0; j < a.Nblocks( ); j++ ) {
    out << "\t" << a.Gaps(j) << "\t" << a.Lengths(j) << "\t" << mutations_by_block[j];
  }
  out << "\n";
}

off_t look_align_plus::BinaryWrite(int fd) const {
  int bytes= look_align::BinaryWritePortable(fd);
  int size = mutations_by_block.size();
  bytes+= SafeWrite(fd, &size, sizeof(size));
  bytes += SafeWrite(fd, &mutations_by_block[0],
		     size * sizeof(mutations_by_block[0]));
  return bytes;
}

off_t look_align_plus::BinaryRead(int fd) {
  int bytes= look_align::BinaryReadPortable(fd);
  int size = 0;
  bytes+= read(fd, &size, sizeof(size));
  ForceAssert(0<=size && size<=100000);
  mutations_by_block.resize(size);
  bytes += read(fd, &mutations_by_block[0],
		size * sizeof(mutations_by_block[0]));
  return bytes;
}


/// Load look_aligns from a binary file.
///  q_ids: if not null, load only hits with given query_id (must be sorted)
///  t_ids: if not null, load only hits with given target_id (must be sorted)
void LoadLookAlignBinary( const String &filename,
			  vec<look_align> &hits,
			  const set<int> *q_ids,
			  const set<int> *t_ids ) {
  int fd = OpenForRead(filename);
  int size;
  ReadBytes(fd, &size, sizeof(int));
  ForceAssert( size>=0 && size < 100000000 );
  hits.clear();
  look_align la;
  for (int i=0; i !=size; ++i) {
    la.BinaryReadPortable(fd);
    if (q_ids && q_ids->find(la.query_id) == q_ids->end()) continue;
    if (t_ids && t_ids->find(la.target_id) == t_ids->end()) continue;
    hits.push_back(la);
  }
  Close(fd);
}

/// Load look_align_plusses from a binary file.
///  q_ids: if not null, load only hits with given query_id (must be sorted)
///  t_ids: if not null, load only hits with given target_id (must be sorted)
void LoadLookAlignPlusBinary( const String &filename,
			  vec<look_align_plus> &hits,
			  const set<int> *q_ids,
			  const set<int> *t_ids ) {
  int fd = OpenForRead(filename);
  LoadLookAlignPlusBinary(fd, hits, q_ids, t_ids);
  Close(fd);
}


/// Write look_aligns to a binary file.
///  q_ids: if not null, load only hits with given query_id
///  t_ids: if not null, load only hits with given target_id
void WriteLookAlignBinary( const String &filename,
			   const vec<look_align> &hits,
			   const set<int> *q_ids,
			   const set<int> *t_ids ) {
  int fd = OpenForWrite(filename);
  int size = hits.size();
  WriteBytes(fd, &size, sizeof(int));
  for (int i=0; i !=size; ++i) {
    if (q_ids && q_ids->find(hits[i].query_id) == q_ids->end()) continue;
    if (t_ids && t_ids->find(hits[i].target_id) == t_ids->end()) continue;
    hits[i].BinaryWritePortable(fd);
  }
  Close(fd);
}

/// Write look_align_plusses to a binary file.
///  q_ids: if not null, load only hits with given query_id
///  t_ids: if not null, load only hits with given target_id
void WriteLookAlignPlusBinary( const String &filename,
			   const vec<look_align_plus> &hits,
			   const set<int> *q_ids,
			   const set<int> *t_ids ) {
  int fd = OpenForWrite(filename);
  WriteLookAlignPlusBinary( fd, hits, q_ids, t_ids );
  Close(fd);
}


/// Write look_align_plusses to an open binary file.
///
/// Parameters:
///  fd - a file descriptor, open for writing
///  q_ids - if not null, load only hits with given query_id
///  t_ids - if not null, load only hits with given target_id
void WriteLookAlignPlusBinary( int fd,
			       const vec<look_align_plus> &hits,
			       const set<int> *q_ids,
			       const set<int> *t_ids ) {
  int size = hits.size();
  WriteBytes(fd, &size, sizeof(int));
  for (int i=0; i !=size; ++i) {
    if (q_ids && q_ids->find(hits[i].query_id) == q_ids->end()) continue;
    if (t_ids && t_ids->find(hits[i].target_id) == t_ids->end()) continue;
    hits[i].BinaryWrite(fd);
  }
}

/// Load look_align_plusses from an open binary file.
///  q_ids: if not null, load only hits with given query_id (must be sorted)
///  t_ids: if not null, load only hits with given target_id (must be sorted)
void LoadLookAlignPlusBinary( int fd,
			      vec<look_align_plus> &hits,
			      const set<int> *q_ids,
			      const set<int> *t_ids ) {
  int size;
  ReadBytes(fd, &size, sizeof(int));
  ForceAssert(size>=0 && size < 10000000);
  hits.clear();
  look_align_plus la;
  for (int i=0; i !=size; ++i) {
    la.BinaryRead(fd);
    if (q_ids && q_ids->find(la.query_id) == q_ids->end()) continue;
    if (t_ids && t_ids->find(la.target_id) == t_ids->end()) continue;
    hits.push_back(la);
  }
}


/**
   Write a vector of vectors of <look_align_pluses> to a binary file.
*/
void WriteVecLookAlignPlusVec( const String& file_name,
			       const vec_look_align_plus_vec& vecLookAlignPlusVec ) {
  int fd = OpenForWrite(file_name);
  for (int i=0; i < vecLookAlignPlusVec.isize(); i++) {
    WriteLookAlignPlusBinary( fd, vecLookAlignPlusVec[i] );
  }
  Close(fd);
}


/**
   Load a vector of vectors of <look_align_pluses> from a binary file.
*/
void LoadVecLookAlignPlusVec( const String& file_name,
			      vec_look_align_plus_vec& vecLookAlignPlusVec ) {
  int fd = OpenForWrite(file_name);
  for (int i=0; i < vecLookAlignPlusVec.isize(); i++) {
    look_align_plus_vec oneVec;
    LoadLookAlignPlusBinary( fd, oneVec );
    vecLookAlignPlusVec.push_back( oneVec );
  }
  Close(fd);
}


/*
 * LoadLookAlignPlus
 */
void LoadLookAlignPlus( const String &file_name,
			vec<look_align_plus> &hits,
			const vec<int> *q_ids,
			const vec<int> *t_ids )
{
  if ( !IsRegularFile( file_name ) )
    FatalErr ( "LoadLookAlignPlus: Can't find LookAlign file '" << file_name << "'.  Exiting.\n" );

  if ( q_ids )
    ForceAssert( is_sorted( q_ids->begin( ), q_ids->end( ) ) );
  if ( t_ids )
    ForceAssert( is_sorted( t_ids->begin( ), t_ids->end( ) ) );

  hits.clear( );
  hits.reserve( LineCount( file_name ) );

  fast_ifstream in( file_name.c_str( ) );
  look_align_plus hit;
  String a_line;
  Bool good=False;
  while ( 1 ) {
    good = getline_if_match( in, a_line, QUERY );
    if ( in.fail() ) break;
    if ( !good ) continue;

    hit.ReadParseable( a_line );
    if ( q_ids )
      if ( ! binary_search( q_ids->begin( ), q_ids->end( ), hit.query_id ) )
	continue;
    if ( t_ids )
      if ( ! binary_search( t_ids->begin( ), t_ids->end( ), hit.target_id ) )
	continue;
    hits.push_back( hit );
  }
}

int LookAlignOffset( const look_align_plus &hit )
{
  const block_align b_align( &hit );
  int best_stretch = 0;
  int best_idx = -1;
  for (int ii=0; ii<(int)b_align.size( ); ii++) {
    int new_stretch = b_align[ii].Len( ) - b_align[ii].Mut( );
    if ( new_stretch >= best_stretch ) {
      best_stretch = new_stretch;
      best_idx = ii;
    }
  }
  ForceAssert( best_idx > -1 );
  return b_align[best_idx].Begin2( ) - b_align[best_idx].Begin1( );
}

bool LookAlignOffset( const look_align_plus &hit1,
		      const look_align_plus &hit2,
		      int &offset )
{
  // int len1, int len2, const align &al1, const align &al2, int &amt )
  offset = 0;

  int len1 = hit1.query_length;
  int len2 = hit2.query_length;
  const align &al1 = hit1.a;
  const align &al2 = hit2.a;
  const block_align b_al1( &hit1 );
  const block_align b_al2( &hit2 );

  // Hanging end amounts.
  int max_hang = 12;

  int left_hang1 = al1.pos1( );
  int left_hang2 = al2.pos1( );
  int right_hang1 = len1 - al1.Pos1( );
  int right_hang2 = len2 - al2.Pos1( );

  // One read embedded in the other, then no hanging ends or leave.
  bool hang1 = ( left_hang1 > max_hang || right_hang1 > max_hang );
  bool hang2 = ( left_hang2 > max_hang || right_hang2 > max_hang );

  bool embedded = false;
  if ( al1.pos2( ) < al2.pos2( ) && al2.Pos2( ) < al1.Pos2( ) )
    embedded = true;
  if ( al2.pos2( ) < al1.pos2( ) && al1.Pos2( ) < al2.Pos2( ) )
    embedded = true;
  if ( embedded && ( hang1 || hang2 ) )
    return false;

  // Either ( read1 is on left and read2 on right ), or vice versa.
  bool order12 = ( al1.pos2( ) < al2.pos2( ) ) ? true : false;

  // No hanging ends on overlap region or leave.
  if ( order12 ) {
    if ( right_hang1 > max_hang || left_hang2 > max_hang )
      return false;
  }
  else {
    if ( right_hang2 > max_hang || left_hang1 > max_hang )
      return false;
  }

  // Search for best overlapping block.
  int best_overlap = 0;
  bool if_overlap = false;
  for (int ii=0; ii<(int)b_al1.size( ); ii++) {
    for (int jj=0; jj<(int)b_al2.size( ); jj++) {
      int loc_overlap = BlocksOverlap( b_al1[ii], b_al2[jj] );
      if ( loc_overlap > best_overlap ) {
	int offset_begin2 = b_al2[jj].Begin2( ) - b_al2[jj].Begin1( );
	int offset_begin1 = b_al1[ii].Begin2( ) - b_al1[ii].Begin1( );
	offset = offset_begin2 - offset_begin1;
	best_overlap = loc_overlap;
	if_overlap = true;
      }
    }
  }

  // Do not trust long overlaps if best block is short.
  int best_block_multiplier = 4;
  const int max_al_len = best_block_multiplier * best_overlap;

  int align_left = Max( al1.pos2( ), al2.pos2( ) );
  int align_right = Min( al1.Pos2( ), al2.Pos2( ) );
  int align_len = align_right - align_left;

  if ( align_len > max_al_len )
    if_overlap = false;

  return if_overlap;
}

int LookAlignOffsetOverlap( const look_align_plus &hit1,
			    const look_align_plus &hit2 )
{
  int offset;

  if ( LookAlignOffset( hit1, hit2, offset ) ) {
    int len1 = hit1.query_length;
    int len2 = hit2.query_length;
    int left = Max( 0, offset );
    int right = Min( offset + len2, len1 );

    return Max( 0, right - left );
  }

  return 0;
}

// return true if the beginning of a query (left side for fw orientation)
// aligns to the target
bool IsLeft( const look_align &align )
{
  if ( !align.rc1 )
    return ( align.a.pos1() == 0 && align.a.Pos1() < (int) align.query_length );
  else
    return ( align.a.pos1() > 0 && align.a.Pos1() == (int) align.query_length );
}

// return true if the end of a query (right side for fw orientation)
// aligns to the target
bool IsRight( const look_align &align )
{
  if ( !align.rc1 )
    return ( align.a.pos1() > 0 && align.a.Pos1() == (int) align.query_length );
  else
    return ( align.a.pos1() == 0 && align.a.Pos1() < (int) align.query_length );
}

// return true if neither end of the query
// aligns to the target
bool IsCenter( const look_align &align )
{
  return ( align.a.pos1() > 0 && align.a.Pos1() < (int) align.query_length );
}


// return the begin, end bases of the query which do not align to the target
pair<int,int> hangingEnd( const look_align &align )
{
  if ( !align.rc1 )
  {
    if ( IsLeft( align ) )
      return ( make_pair(align.a.Pos1(), align.query_length) );
    else if ( IsRight( align ) )
      return ( make_pair(0, align.a.pos1()) );
    else return ( make_pair( -1, -1 ) );
  }
  else
  {
    if ( IsLeft( align ) )
      return ( make_pair( align.query_length-align.a.pos1(), align.query_length) );
    else if ( IsRight( align ) )
      return ( make_pair(0, align.query_length-align.a.Pos1()) );
    else return ( make_pair( -1, -1 ) );
  }
}

// return the maximum hang amount (largest between head and tail)
int MaxHang( const look_align &hit )
{
  int len1 = hit.query_length;
  int len2 = hit.target_length;
  int head = Min( hit.a.pos1( ), hit.a.pos2( ) );
  int tail = Min( len1 - hit.a.Pos1( ), len2 - hit.a.Pos2( ) );

  return Max( head, tail );
}

// return size of largest gap
int MaxGap( const look_align &hit )
{
  int largest_gap = 0;
  for (int jj=0; jj<(int)hit.a.Nblocks( ); jj++)
    largest_gap = Max( largest_gap, Abs( hit.a.Gaps( jj ) ) );

  return largest_gap;
}

// a pair is logical iff hits belong to same target and have different orient.
bool IsLogicalPair( const look_align &hit1,
		    const look_align &hit2 )
{
  // Hits belong to different targets.
  if ( ! ( hit1.target_id == hit2.target_id ) )
    return false;

  // Same orientation.
  if ( ! ( hit1.rc1 != hit2.rc1 ) )
    return false;

  // Ok.
  return true;
}

// an insert is valid if it is logical and the two end reads are at the
// right distance from each other. The "right distance" is determined
// by one or both of the following tests (see also Stretch( ) below):
//   (t1)  -max_stretch <= stretch <= max_stretch (max_stretch>0), and:
//   (t2)   1/max_mult <= mult_stretch <= max_mult (max_mult>1)
// Remarks:
// * Test (t1) will be run if max_stretch > 0, test (t2) if max_mult > 1.
// * At least one of (t1) and (t2) must be run, if both are run, then a
//   pair is declared valid if it satisfies (t1) and/or (t2).
bool IsValidPair( const look_align &hit1,
		  const look_align &hit2,
		  const read_pairing &pair,
		  double max_stretch,
		  double max_mult )
{
  ForceAssert( max_mult > 1.0 || max_stretch > 0 );

  if ( ! IsLogicalPair( hit1, hit2 ) )
    return false;

  if ( max_stretch > 0 ) {
    float stretch = Stretch( hit1, hit2, pair );
    if ( - max_stretch <= stretch && stretch <= max_stretch )
      return true;
  }

  if ( max_mult > 1.0 ) {
    float mult_stretch = Stretch( hit1, hit2, pair, true );
    if ( 1.0 / max_mult <= mult_stretch && mult_stretch <= max_mult )
      return true;
  }

  return false;
}

// return observed separation (warning! If pair fails IsLogicalPair
// Stretch will assert)
int ObservedSeparation( const look_align &hit1,
			const look_align &hit2 )
{
  // Check input.
  ForceAssert( IsLogicalPair( hit1, hit2 ) );

  // Observed separation.
  int left = 0;
  int right = 0;
  if ( hit2.rc1 ) {
    right = hit2.a.pos2( ) - hit2.a.pos1( );
    left = hit1.a.Pos2( ) + ( hit1.query_length - hit1.a.Pos1( ) );
  }
  else {
    right = hit1.a.pos2( ) - hit1.a.pos1( );
    left = hit2.a.Pos2( ) + ( hit2.query_length - hit2.a.Pos1( ) );
  }
  int sep_observ = right - left;

  // Return.
  return sep_observ;
}

// return stretch (warning! If pair fails IsLogicalPair Stretch will assert)
float Stretch ( const look_align &hit1,
		const look_align &hit2,
		const read_pairing &pair,
		bool as_multiplier )
{
  // Observed separation (will also check input).
  int sep_observ = ObservedSeparation( hit1, hit2 );
  int sep_given = pair.sep;

  // Return stretch as multiplier (handle the special case sep_given = 0).
  if ( as_multiplier ) {
    if ( sep_given == 0 )
      return sep_observ == 0 ? 1.0 : SafeQuotient( sep_observ, 1 );
    return SafeQuotient( sep_observ, sep_given );
  }

  // Return conventional stretch.
  ForceAssert( pair.sd > 0 );
  return SafeQuotient( sep_observ - sep_given, pair.sd );

}

void Trim1Together(const basevector & b1, const basevector & b2,
		   const look_align & la, int startOn1, int len,
	      basevector & trimmedb1, look_align & trimmedla) {
  trimmedb1.SetToSubOf(b1,startOn1, len);
  trimmedla = la.TrimmedTo1(startOn1, len, b1, b2);
}

void Trim1Together(const basevector & b1, const basevector & b2,
		   const look_align & la, int startOn1, int len,
		   basevector & trimmedb1, look_align_plus & trimmedla){
  trimmedb1.SetToSubOf(b1,startOn1, len);
  trimmedla = la.TrimmedTo1Plus(startOn1, len, b1, b2);
}

uint CountLookAligns( const String & fname ) {
  fast_ifstream in(fname);
  uint count = 0;
  String line;
  while (true) {
    if (getline_if_match(in, line, QUERY)) ++count;
    if (in.fail()) break;
  }
  return count;
}

void LoadLookAligns( const String& file_name, vec<look_align>& aligns,
		      vec< vec<align_id_t> >& aligns_index, unsigned long nqueries ) {
  aligns.resize(CountLookAligns(file_name));
  aligns_index.clear( );
  aligns_index.resize(nqueries);
  int aligns_count = 0;
  String line;
  fast_ifstream in(file_name);
  bool good = false;
  while(1) {
    good = getline_if_match( in, line , QUERY);
    if ( in.fail( ) ) break;
    if ( good ) {
      aligns[aligns_count].ReadParseable(line);
      int id = aligns[aligns_count].query_id;
      aligns_index[id].push_back(aligns_count);
      ++aligns_count;
      // if ( !( aligns_count % 4000000 ) ) DPRINT( aligns_count );
    }
  }
}

void LoadLookAligns( const String& file_name, vec<look_align>& aligns,
     const Bool ignore_bads ) 
{
  aligns.resize(CountLookAligns(file_name));
  String line;
  int64_t aligns_count=0;
  fast_ifstream in(file_name);
  bool good = false;
  while(1) {
    good = getline_if_match( in, line, QUERY );
    if ( in.fail( ) ) break;
    if ( good ) 
    {    if ( !ignore_bads ) aligns[aligns_count++].ReadParseable(line);
         else
         {    Bool read_ok = aligns[aligns_count].ReadParseableOrFail(line);
              if (read_ok) aligns_count++;    }    }
  }
  aligns.resize(aligns_count);
}


void SetWritePrettyLookAligns(ofstream & os)
{
  if ( parsed_args::ARGS ) {
    parsed_args& command = (parsed_args&)(*parsed_args::ARGS);
    command.PrintTheCommandPretty(os);
  }
}


void WriteLookAligns( const String& file_name, const vec<look_align>& aligns ) 
{
  Ofstream(impltout, file_name );

  SetWritePrettyLookAligns(impltout);

  for ( align_id_t i=0; i != aligns.isize(); ++i ) {
      aligns[i].PrintParseable(impltout);
      aligns[i].PrintReadableBrief(impltout);
  }
}


ofstream & operator << (ofstream & os, const look_align & la)
{
  la.PrintParseable(os);
  la.PrintReadableBrief(os);
  return os;
}


void BuildLookAlignsIndex( const vec<look_align>& aligns,
			   vec< vec<align_id_t> >& aligns_index, int nqueries ) {
  aligns_index.clear( );
  aligns_index.resize(nqueries);
  for (int i = 0; i < aligns.isize(); i++ ) {
    int id = aligns[i].query_id;
    aligns_index[id].push_back(i);
  }
}

void RemoveOverlap(look_align & la1, look_align & la2, bool verbose) {

  if (verbose) {
    la1.PrintParseable(cout);
    la2.PrintParseable(cout);
  }
  int startOn1 = la1.rc1 ? la1.query_length - la1.a.Pos1()
    : la1.a.pos1();
  int startOn2 = la2.rc1 ? la2.query_length - la2.a.Pos1()
    : la2.a.pos1();

  if (startOn1 == startOn2 ) { //remove the shorter of the two alignments.
    if (la1.a.extent1() < la2.a.extent1()) swap(la1,la2);
    la2.a = align();
    la2.mutations = la2.indels = la2.nhits = 0;
    return;
  }
  else if (startOn1 > startOn2) {
    swap(la1, la2);
    swap(startOn1, startOn2);
  }
  if (verbose) PRINT2(startOn1, startOn2);
  if (startOn1 >= startOn2) {
    la1.PrintReadableBrief(cout);
    la2.PrintReadableBrief(cout);
    PRINT2(startOn1, startOn2);
    FatalErr("startOn1 should be less than startOn2");
  }

  //Calculate the endpoint for trimming
  int newend = startOn2 - 1;
  align a;
  //Adjust the end of la1. Note that the indels
  //and mutations members of la1 become invalid.
  if (!la1.rc1) {
    a = la1.a.TrimmedTo1(0,newend);
  } else {
    a = la1.a.TrimmedTo1(la1.query_length - newend, startOn2);
    a.Setpos1(la1.query_length - newend);
  }
  la1.a = a;
  if (verbose) {
    la1.PrintParseable(cout);
    la2.PrintParseable(cout);
  }
}

look_align FromAlignmentPlus(const alignment_plus & alp, int query_length,
			     int target_length) {
  int tid = alp.Id1();
  int qid = alp.Id2();
  alignment alg = alp.a;

  //look_aligns reverse the first sequence, whereas alignment_pluses reverse the
  //second. So we need to reverse the whole thing to do the transformation.
  if ( alp.Rc2() ) {
    alg.ReverseThis( query_length,  target_length);
  }
  align a;
  a.UnpackFrom(alg);
  return look_align( qid, tid, query_length, target_length, alp.Rc2(), a, 0, 0, 0 );
}

void GetBestAligns
(
 const vecbasevector& bases,             // the reads
 const vec<look_align>& aligns,          // alignments of reads to reference
 const vec< vec<int> >& aligns_index,    // index by read of alignments
 vec<int>& best                          // index of best alignment (returned)
 )
{
  best.resize_and_set( bases.size( ), -1 );
  for ( size_t id = 0; id < bases.size( ); id++ ) {
    int min_errors = 1000000000;
    vec<int> minerr;
    minerr.clear( );
    for ( int j = 0; j < aligns_index[id].isize( ); j++ ) {
      int errs = aligns[ aligns_index[id][j] ].Errors( );
      if ( errs < min_errors ) minerr.clear( );
      if ( errs <= min_errors ) {
        minerr.push_back(j);
        min_errors = errs;
      }
    }
    if ( minerr.nonempty( ) )
      best[id] = minerr[ ( minerr.solo( ) ? 0 : ( randomx( ) % minerr.size( ) ) ) ];
  }
}

template<int> void BinaryWrite3( const String& filename, const vec<int>& v );
template<int> void BinaryRead3( const String& filename, vec<int>& v, bool strict = false );


BINARY2_DEF(GaplessAlign);
BINARY3_DEF(GaplessAlign);


