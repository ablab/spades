///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "paths/reporting/ReportWalkingRatesCore.h"

/**
 * ReportWalkingRatesCore
 */
void ReportWalkingRatesCore( ostream &out,
			     const bool VERBOSE,
			     const vec<int> &ids,
			     const vec< pair<int,int> > &info )
{
  // Counters.
  longlong tot_walked = 0;
  longlong tot_tried = 0;
  longlong n_empty = 0;
  
  // Parse info.
  for (int nhood_id=0; nhood_id<ids.isize( ); nhood_id++) {
    int n_walked = info[nhood_id].first;
    int n_tried = info[nhood_id].second;
    if ( n_tried < 0 ) {
      if ( VERBOSE )
	cout << nhood_id << "\tu_"
	     << ids[nhood_id] << "\tmissing info\n";
      continue;
    }
    tot_walked += n_walked;
    tot_tried += n_tried;
    if ( n_walked == 0 ) n_empty++;
    
    if ( ! VERBOSE ) continue;
   
    double ratio = -1;
    if ( n_tried > 0 ) ratio = 100.0 * SafeQuotient( n_walked, n_tried );
    String str_ratio = n_tried > 0 ? ToString( ratio, 1 ) + "%" : "na";
    out << nhood_id << "\tu_"
	 << ids[nhood_id] << "\tn_walked: "
	 << n_walked << " out of "
	 << n_tried << "\tsuccess_rate: "
	 << str_ratio << "\n";
  }
  if ( VERBOSE ) out << endl;

  // Print summary.
  double s_ratio = -1.0;
  if ( tot_tried > 0 ) s_ratio = 100.0 * SafeQuotient( tot_walked, tot_tried );

  out << "Insert walking statistics:\n"
      << "\n"
      << "  local nhoods :   " << ToStringAddCommas( ids.size( ) ) << "\n"
      << "  inserts total:   " << ToStringAddCommas( tot_tried ) << "\n"
      << "  inserts walked:  " << ToStringAddCommas( tot_walked ) << "\n"
      << "  percent walked:  " << ToString( s_ratio, 1 ) << "%\n"
      << "\n";
  
  if ( tot_tried < 1 ) {
    out << "No input found.\n" << endl;
    return;
  }
  
  const int n_bins = 11;
  vec<int> walked_histo( n_bins, 0 );
  for (int ii=0; ii<info.isize( ); ii++) {
    if ( info[ii].first < 0 ) continue;
    int nhood_walked = Min( info[ii].first, n_bins - 1 );
    walked_histo[nhood_walked] += 1;
  }
  
  vec< vec<String> > table;
  vec<String> line;
  
  line = MkVec( String( "n walked" ),
		String( "n nhoods" ),
		String( "%" ),
		String( "cumul %" ) );
  table.push_back( line );

  longlong cumul_total = 0;
  longlong n_nhoods = ids.size( );
  for (int ii=0; ii<n_bins; ii++) {
    cumul_total += walked_histo[ii];
    double bin_ratio = 100.0 * SafeQuotient( walked_histo[ii], n_nhoods );
    double cumul_ratio = 100.0 * SafeQuotient( cumul_total, n_nhoods );
    String str_plus = ( ii == n_bins - 1 ) ? "+" : "";

    line.clear( );
    line.push_back( ToStringAddCommas( ii ) + str_plus );
    line.push_back( ToStringAddCommas( walked_histo[ii] ) );
    line.push_back( ToString( bin_ratio, 1 ) );
    line.push_back( ToString( cumul_ratio, 1 ) );
     table.push_back( line );
   }
  
  out << "Histogram of inserts walked per nhood:\n\n";
  PrintTabular( out, table, 3, "rrrr" );
  out << endl;
  
}
