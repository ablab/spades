// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "TrimsFile.h"

/*
 * Load trims in the specified vectors.
 */
void LoadTrims( const String &trims_file,
		vec<int> &left_trims,
		vec<int> &right_trims )
{
  left_trims.clear( );
  right_trims.clear( );

  ForceAssert( IsRegularFile( trims_file ) );

  int n_reads = LineCount( trims_file );
  left_trims.reserve( n_reads );
  right_trims.reserve( n_reads );

  Ifstream( trims_stream, trims_file );
  
  for (int ii=0; ii<n_reads; ii++) {
    int left_amount;
    int right_amount;

    trims_stream >> left_amount >> right_amount;
    left_trims.push_back( left_amount );
    right_trims.push_back( right_amount );
  }

  trims_stream.close( );
}
