// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

// Generates indices relative to pairing.

#ifndef PAIRS_HANDLER_H
#define PAIRS_HANDLER_H

#include "ReadPairing.h"
#include "String.h"

class phandler {

public:

  phandler( int n_reads ) : n_reads_ ( n_reads ) { }

  phandler( int n_reads, const String &pairs_file );

  void LoadFromFile( const String &pairs_file );

  int GetPartnerId( int read_id ) const { return pairs_to_[read_id]; }

  int GetPairId( int read_id ) const { return pair_id_[read_id]; }

  // This returns null if read_id is not paired.
  const read_pairing *GetPair( int read_id ) const;
  
  // Actual pairs (mostly for backward compatibility).
  const vec<read_pairing> &Pairs( ) const { return pairs_; }

  // Map sending read_id to id of partner.
  const vec<int> &PairsTo( ) const { return pairs_to_; }

  // Similary, the map sending read_id to the id of pair containing it.
  const vec<int> &PairId( ) const { return pair_id_; }

  // Returns pair at ii.
  const read_pairing &operator[]( int ii ) const { return pairs_[ii]; }
  
  
private:

  int n_reads_;              // needed to resize pairs_to_ and pair_id
  vec<read_pairing> pairs_;  // the pairs
  vec<int> pairs_to_;        // id of partner (or -1)
  vec<int> pair_id_;         // id of pair (or -1)

};

#endif
