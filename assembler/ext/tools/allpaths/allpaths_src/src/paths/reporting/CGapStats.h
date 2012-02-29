///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef C_GAP_STATS
#define C_GAP_STATS

#include "SupersHandler.h"
#include "Vec.h"
#include "math/Functions.h"

/**
 * class CGapStats
 *
 * Report on gap accuracy.
 */
class CGapStats {

public:

  // supers are only used to reserve memory.
  CGapStats( const shandler *supers = 0 ) :
    sigma_ ( 0.0 ),
    n_supers_ ( 0 ),
    n_gaps_in_supers_ ( 0 ) {
    if ( supers ) {
      this->Reserve( supers->NContigs( ) );
      n_supers_ = supers->Size( );
      n_gaps_in_supers_ = supers->NContigs( ) - supers->Size( );
    }
  }
  
  void Reserve( int ngaps ) { gaps_.reserve( ngaps ); }
  
  void SetSigma( double sigma ) { sigma_ = sigma; }

  void SetNSupers( int n_supers ) { n_supers_ = n_supers; }

  void SetNGapsInSupers( int n_gaps ) { n_gaps_in_supers_ = n_gaps; }
  
  void Add( int observed, int gap_size, int gap_dev ) {
    gaps_.push_back( triple<int,int,int>( observed, gap_size, gap_dev ) );
  }

  // How many gaps in dataset.
  int Size( ) const { return gaps_.isize( ); }

  // Discrepancy with observed gap.
  int Delta( int ii ) const {
    int observed = gaps_[ii].first;
    int gsize = gaps_[ii].second;
    int gdev = gaps_[ii].third;
    
    int rad = sigma_ * (double)gdev;
    int winl = gsize - rad;
    int winr = gsize + rad;
    
    if ( winl > observed ) return ( winl - observed );   // gap is bigger
    if ( winr < observed ) return ( winr - observed );   // gap is smaller
    return 0;   // on target
  }

  // Sum of Abs( delta ) for all gaps.
  unsigned int TotalDelta( int &smaller, int &equal, int &bigger ) const {
    uint total = 0;
    smaller = 0;
    equal = 0;
    bigger = 0;
    
    for (size_t ii=0; ii<gaps_.size( ); ii++) {
      int delta = this->Delta( ii );
      total += Abs( delta );
      if ( delta < 0 ) smaller++;
      else if ( delta > 0 ) bigger++;
      else equal++;
    }

    return total;
  }
  
  // Set sigma and generate table-line report for the given sigma (or legend).
  vec<String> TableLine( const double sigma, bool legend = false ) const {
    sigma_ = sigma;
    
    vec<String> tline;
    if ( legend ) {
      tline.push_back( "sigma" );
      tline.push_back( "total_delta" );
      tline.push_back( "smaller" );
      tline.push_back( "on_target" );
      tline.push_back( "bigger" );

      return tline;
    }
    
    int smaller = 0;
    int equal = 0;
    int bigger = 0;
    uint tot = this->TotalDelta( smaller, equal, bigger );

    tline.push_back( ToString( sigma, 1 ) );

    tline.push_back( ToString( tot ) );

    int n_gaps = gaps_.isize( );
    double ratio = 100.0 * SafeQuotient( smaller, n_gaps );
    tline.push_back( ToString( smaller ) + " (" + ToString( ratio, 2 ) + "%)" );

    ratio = 100.0 * SafeQuotient( equal, n_gaps );
    tline.push_back( ToString( equal ) + " (" + ToString( ratio, 2 ) + "%)" );

    ratio = 100.0 * SafeQuotient( bigger, n_gaps );
    tline.push_back( ToString( bigger ) + " (" + ToString( ratio, 2 ) + "%)" );
    
    return tline;
  }

  void PrintReport( ostream &out ) const {
    int n_gaps = gaps_.isize( );
    
    out << "GAP STATISTICS (" << n_gaps
	<< " valid gaps, over a total of " << n_gaps_in_supers_
	<< " gaps in " << n_supers_
	<< " supers):\n" << endl;
   
    if ( n_gaps < 1 ) {
      out << " No valid gaps found.\n" << endl;
      return;
    }

    vec< vec<String> > table;
    table.push_back( this->TableLine( 0.0, true ) );
    for (int ii=0; ii<6; ii++)
      table.push_back( this->TableLine( ( ii * 1.0 ) ) );

    PrintTabular( out, table, 4, "rrrrr" );
    out << endl;
  }
  
  
private:

  mutable double sigma_;              // allowed stretch for reported gaps
  int n_supers_;                      // total number of supers
  int n_gaps_in_supers_;              // total numer of gaps in supers
  vec< triple<int,int,int> > gaps_;   // observed, gap_size, gap_dev
  
};

#endif
