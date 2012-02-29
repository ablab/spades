// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#include "PairsHandler.h"
#include "ReadPairing.h"
#include "String.h"

/*
 * phandler
 * Constructor
 */
phandler::phandler( int n_reads, const String &pairs_file ) :
  n_reads_ ( n_reads )
{
  this->LoadFromFile( pairs_file );
}

/*
 * phandler
 * LoadFromFile
 */
void phandler::LoadFromFile( const String &pairs_file )
{
  READX( pairs_file, pairs_ );

  pairs_to_.clear( );
  pair_id_.clear( );
  pairs_to_.resize( n_reads_, -1 );
  pair_id_.resize( n_reads_, -1 );
  for (int ii=0; ii<(int)pairs_.size( ); ii++) {
    pairs_to_[pairs_[ii].id1] = pairs_[ii].id2;
    pairs_to_[pairs_[ii].id2] = pairs_[ii].id1;
    pair_id_[pairs_[ii].id1] = ii;
    pair_id_[pairs_[ii].id2] = ii;
  }
}

/*
 * phandler
 * GetPair
 */
const read_pairing *phandler::GetPair( int read_id ) const
{
  int pair_id = this->GetPairId( read_id );
  if ( pair_id < 0 )
    return 0;
  return &pairs_[pair_id];
}
