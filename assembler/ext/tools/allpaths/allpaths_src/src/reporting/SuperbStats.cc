/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "ParseSet.h"
#include "Superb.h"
#include "Vec.h"
#include "math/Functions.h"
#include "math/NStatsTools.h"
#include "system/ParsedArgs.h"
#include "system/RunTime.h"
#include "system/System.h"

/**
 * StrRange
 *
 * It assumes lens is sorted (this is not checked).
 */
String StrRange( const vec<int> &lens )
{
  if ( lens.size( ) < 1 ) return "na";
  String strF = ToString( lens[0] );
  String strL = ToString( lens[ lens.size( )-1 ] );
  return "[" + strF + ", " + strL + "]";
}

/**
 * return string of Median value, return na if not apllicable
 * assuming sorted
 */
String StrMedian (const vec<int> &lens )
{
  return lens.size()>0 ? ToString(Median(lens)): "na";
}

/**
 * StrSafeN50
 *
 * It assumes lens is sorted (this is not checked). Remove 0's from
 * the vector, and call the standard N50 on what is left.
 */
String StrSafeN50( const vec<int> &lens )
{
  vec<int> lplus;
  lplus.reserve( lens.size( ) ) ;
  for (size_t ii=0; ii<lens.size( ); ii++)
    if ( lens[ii] > 0 )
      lplus.push_back( lens[ii] );

  if ( lplus.size( ) < 1 ) return "na";
  return ToString( N50( lplus ) );
}

/**
 * StrAiriness
 *
 * Average of airiness for scaffolds (optionally weighted by scaffold
 * length).
 */
String StrAiriness( const vec<double> &airiness, const vec<int> *glens = 0 )
{
  vec<double> wdata;
  if ( glens ) {
    double sum = BigSum( *glens );
    wdata.reserve( airiness.size( ) );
    for (int ii=0; ii<airiness.isize( ); ii++) {
      double weight = double( (*glens)[ii] ) / sum;
      wdata.push_back( weight * airiness[ii] );
    }
  }
  const vec<double> &data = glens ? wdata : airiness;
  return ToString( 100.0 * Mean( data ), 2 );
}

/**
 * SuperbStats
 *
 * Basic statistics for a superb file.
 *
 * SELECT: report stats for these super ids only (if given)
 * COMPLEMENT: report total for the complement of SELECT (if SELECT is given)
 * ARCHIVE: send output to file (<SUPERBS>.SuperbStats), rather than cout
 * NSTATS: if true, print N{98,75,50,25} for contigs and supers
 * AIRINESS: if true, print stats on airiness of supers (gaps in bases)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( SUPERBS );
  CommandArgument_String_OrDefault( SELECT, "" );
  CommandArgument_Bool_OrDefault( COMPLEMENT, False );
  CommandArgument_Bool_OrDefault( ARCHIVE, False );
  CommandArgument_Bool_OrDefault( AIRINESS, True);
  CommandArgument_Bool_OrDefault( NSTATS, False );
  EndCommandArguments;


  String archive_file =  SUPERBS + ".SuperbStats";
  ofstream archive_out;
  if ( ARCHIVE ) archive_out.open( archive_file.c_str( ) );
  ostream &out = ARCHIVE ? archive_out : * (ostream *) &cout;
  if ( ARCHIVE ) {
    cout << "Sending output to " << archive_file << "\n" << endl;
    PrintCommandPretty( out );
  }

  vec<superb> supers;
  ReadSuperbs( SUPERBS, supers );

  vec<Bool> select;
  {
    if ( SELECT == "" )
      select.resize( supers.size( ), True );
    else {
      vec<int> ids;
      ParseIntSet( SELECT, ids );
      select.resize( supers.size( ), COMPLEMENT ? True : False );
      for (int ii=0; ii<ids.isize( ); ii++)
	select[ ids[ii] ] = COMPLEMENT ? False : True;
    }
  }

  int n_selected = 0;
  int n_contigs = 0;
  for (int ii=0; ii<(int)supers.size( ); ii++)
    n_contigs += supers[ii].Ntigs( );

  vec<int> glens; // gapped scaffold lengths
  vec<int> ulens; // ungapped scafffold lengths
  vec<double> airiness; // gap_size / scaffold_length
  vec<int> clens; // contig lengths
  vec<int> gaps;  // gap sizes
  vec<int> gaps_pos; // >0 gaps 
  vec<int> gaps_neg; // <=0 gaps  
  glens.reserve( supers.size( ) );
  ulens.reserve( supers.size( ) );
  airiness.reserve( supers.size( ) );
  clens.reserve( n_contigs );
  gaps.reserve( n_contigs );
  gaps_pos.reserve( n_contigs );
  gaps_neg.reserve( n_contigs );
  longlong sum_overlap_gap = 0; // extremely overlapped gap lengths, larger than 5 * dev 

  for (int ii=0; ii<(int)supers.size( ); ii++) {
    if ( ! select[ii] ) continue;
    n_selected++;
    glens.push_back( supers[ii].TrueLength( ) );
    ulens.push_back( supers[ii].ReducedLength( ) );
    //contigs
    for ( int jj = 0; jj < supers[ii].Ntigs( ); jj++ )
      clens.push_back(  supers[ii].Len(jj) );
    //gaps
    double in_gaps = 0;
    for (int jj=0; jj<supers[ii].Ntigs( )-1; jj++) {
      int gap = supers[ii].Gap( jj );
      int dev = supers[ii].Dev( jj );
      gaps.push_back(gap);
      if ( gap > 0) {
	gaps_pos.push_back(gap);
	in_gaps += double( gap );
      }
      else {
	gaps_neg.push_back(-gap);
	if ( gap < -5 * dev) { // extremely negative gap lengths
	  sum_overlap_gap += - gap;
	}
      }
    }
    airiness.push_back( in_gaps / double(supers[ii].TrueLength( ) ) );
  }
  
  // extreme negative gap lengths per Mb
  longlong overlap_gap_Mb = sum_overlap_gap * 1000000 / BigSum(glens)  ;
  sort( glens.begin( ), glens.end( ) );
  sort( ulens.begin( ), ulens.end( ) );

  sort( clens.begin( ), clens.end( ) );
  sort( gaps.begin( ), gaps.end( ) );
  sort( gaps_pos.begin( ), gaps_pos.end( ) );
  sort( gaps_neg.begin( ), gaps_neg.end( ) );
  
  if ( supers.size() < 1 )
    return 0;

  cout << "negative gaps per Mb (after 5 dev): " << overlap_gap_Mb << "\n";
  if ( AIRINESS ) {
    String str_air = StrAiriness( airiness );
    String str_wair = StrAiriness( airiness, &glens );
    String str_base = "% of bases in gaps ";
    cout << str_base << "(average per scaffold): " << str_air << "\n"
	 << str_base << "(weighted average per scaffold): " << str_wair << "\n";
  }
  cout << endl;

  // Printable table.
  vec< vec<String> > table;
  vec<String> line;
  line.push_back( "" );
  line.push_back( "total_length" );
  line.push_back( "#" );
  line.push_back( "N50" );
  line.push_back( "range" );
  line.push_back( "median" );
  table.push_back( line );

  line.clear();
  line.push_back( "scaffold_gapped" );
  line.push_back( ToString( BigSum( glens ) ) );
  line.push_back( ToString( glens.size()) );
  line.push_back( StrSafeN50( glens ) );
  line.push_back( StrRange( glens ) );
  line.push_back( StrMedian( glens ) );
  table.push_back( line );

  line.clear();
  line.push_back( "scaffold_ungapped" );
  line.push_back( ToString( BigSum( ulens ) ) );
  line.push_back( ToString( ulens.size()) );
  line.push_back( StrSafeN50( ulens ) );
  line.push_back( StrRange( ulens ) );
  line.push_back( StrMedian( ulens ) );
  table.push_back( line );

  line.clear();
  line.push_back( "contigs" );
  line.push_back( ToString( BigSum( clens ) ) );
  line.push_back( ToString( clens.size()) );
  line.push_back( StrSafeN50( clens ) );
  line.push_back( StrRange( clens ) );
  line.push_back( StrMedian( clens ) );
  table.push_back( line );

  line.clear();
  line.push_back( "gap  >0");
  line.push_back( ToString( BigSum( gaps_pos ) ) );
  line.push_back( ToString( gaps_pos.size()) );
  line.push_back( StrSafeN50( gaps_pos ) );
  line.push_back( StrRange( gaps_pos ) );
  line.push_back( StrMedian( gaps_pos ) );
  table.push_back( line );

  line.clear();
  line.push_back( "gap <=0" );
  line.push_back( ToString( BigSum( gaps_neg ) ) );
  line.push_back( ToString( gaps_neg.size()) );
  line.push_back( StrSafeN50( gaps_neg ) );
  line.push_back( StrRange( gaps_neg ) );
  line.push_back( StrMedian( gaps_neg ) );
  table.push_back( line );

  line.clear();
  line.push_back( "gap all");
  line.push_back( ToString( BigSum( gaps ) ) );
  line.push_back( ToString( gaps.size()) );
  line.push_back( "na" );
  line.push_back( StrRange( gaps ) );
  line.push_back( StrMedian( gaps ) );
  table.push_back( line );
  
  // Print.
  PrintTabular( out, table, 2, "lrrrr" );
  out << endl;

  if ( NSTATS ) {
    PrintBasicNStats( "contigs", clens, cout );
    cout << endl;
    
    PrintBasicNStats( "ungapped_scaffolds", ulens, cout );
    cout << endl;
    
    PrintBasicNStats( "gapped_scaffolds", glens, cout );
    cout << endl;
  }

  cout << Date( ) << ": done" << endl;
}
