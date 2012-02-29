///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__IDENTIFY_CIRCULAR_SCAFFOLDS_H
#define PATHS__IDENTIFY_CIRCULAR_SCAFFOLDS_H

#include "Fastavector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "math/HoInterval.h"
#include "paths/Alignlet.h"
#include "paths/ReadLoc.h"

/**
 * class CRing
 *
 * Data structure for circular scaffolds.
 */
class CRing {

public:

  CRing( ) :
    super_id_ ( -1 ), gap_size_( 0 ), gap_dev_ ( -1 ), n_links_ ( 0 )
  { }
	     
  CRing( int super_id, int gap_size, int gap_dev, int n_links ) :
    super_id_ ( super_id ),
    gap_size_ ( gap_size ),
    gap_dev_ ( gap_dev ),
    n_links_ ( n_links )
  { }

  int SuperId( ) const { return super_id_; }
  int GapSize( ) const { return gap_size_; }
  int GapDev( ) const { return gap_dev_; }
  int NLinks( ) const { return n_links_; }
  
  void GenerateLegend( vec<String> &legend ) const {
    legend.clear( );
    legend = MkVec( ToString( "scaffold_id" ),
		    ToString( "gap_size" ),
		    ToString( "gap_dev" ),
		    ToString( "n_links" ) );
  }
  
  void PrintableInfo( vec<String> &printable ) const {
    printable.clear( );
    printable = MkVec( ToString( super_id_ ),
		       ToString( gap_size_ ),
		       ToString( gap_dev_ ),
		       ToString( n_links_ ) );
  }
  
  friend bool operator< ( const CRing &left, const CRing &right ) {
    return ( left.SuperId( ) < right.SuperId( ) );
  }
  
  
private:

  int super_id_;
  int gap_size_;
  int gap_dev_;
  int n_links_;

};

/**
 * IdentifyCircularScaffolds
 *
 * Identify circular scaffolds, and disconnect circular scaffolds from
 * all other scaffolds (this is done by resetting all the indices of
 * reads pointing away from circular scaffolds to 0).
 *
 * is_circular: if super[ii] is tagged as circular
 * VERBOSE: 0 (no log) to 2 (max log)
 */
void IdentifyCircularScaffolds( const PairsManager &pairs,
				const vec<fastavector> &contigs,
				const vec<superb> &supers,
				const vec<alignlet> &aligns,
				vec<int> &index,
				vec<Bool> &is_circular,
				ostream *log = 0,
				int VERBOSE = 0 );

/**
 * IdentifyCircularScaffolds
 *
 * This implementation relies on readlocs, rather than alignlets and
 * index, and hence it does not reset indices to cut away circular
 * scaffolds from the other scaffolds. Output is stored as a vec of
 * triples (super_id, gap, dev).
 *
 * NUM_THREADS: if negative, use all available
 */
void IdentifyCircularScaffolds( const vec<PairsManager> &pairs,
				const vec<fastavector> &contigs,
				const vec<superb> &supers,
				read_locs_on_disk &locs_file,
				vec<CRing> &rings,
				int NUM_THREADS = -1,
				ostream *log = 0 );

#endif
