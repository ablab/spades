///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SUPERS_HANDLER_H
#define SUPERS_HANDLER_H

#include "String.h"
#include "Superb.h"

/**
 * class shandler
 *
 * Generates indices relative to supers.
 */
class shandler {

public:

  // n_contigs can be set to -1 (see .cc, LoadFromFile( ) for details).
  shandler( int n_contigs );

  // n_contigs can be set to -1 (see .cc, LoadFromFile( ) for details).
  shandler( int n_contigs, const String &supers_file );

  // If min_gap is given, reset gaps to >= *min_gap (stdevs are not changed).
  shandler( const vec<superb> &supers, int *min_gap = 0 );

  // If min_gap is given, reset gaps to >= *min_gap (stdevs are not changed).
  void LoadFromFile( const String &supers_file, int *min_gap = 0 );
  
  // Supers count.
  int Size( ) const { return supers_.size( ); }
  
  // Contigs count.
  int NContigs( ) const { return n_contigs_; }

  // True length of super id (TrueLength is defined in Superb.h).
  int TrueLength( int id ) const { return true_ends_[id] - true_begins_[id]; }

  // Returns id of super containing cg_id.
  int ToSuper( int cg_id ) const { return to_super_[cg_id]; }

  // Start on super of cg_id.
  int StartOnSuper( int cg_id ) const { return start_on_super_[cg_id]; }

  // Position in super of cg_id.
  int PosOnSuper( int cg_id ) const { return pos_on_super_[cg_id]; }
  
  // Returns whole to_super map (mostly for backward compatibility).
  const vec<int> &ToSuper( ) const { return to_super_; }

  // Full map for start on supers.
  const vec<int> &StartOnSuper( ) const { return start_on_super_; }

  // Full map for position on supers.
  const vec<int> &PosOnSuper( ) const { return pos_on_super_; }
  
  // Returns super at ii.
  const superb &operator[]( int ii ) const { return supers_[ii]; }

  // Return all supers.
  const vec<superb> &AllSupers( ) const { return supers_; }
  
  
private:

  void Setup( int *min_gap = 0 );
  
  
private:

  int n_contigs_;           // needed to resize to_super_ and other maps
  vec<superb> supers_;      // the actual supers
  vec<int> true_begins_;    // true begin (see Superb.h) of supers
  vec<int> true_ends_;      // true end (see Superb.h) of supers
  vec<int> to_super_;       // map contig id to id of super containing it
  vec<int> start_on_super_; // start on super of given contig
  vec<int> pos_on_super_;   // position in super of given contig
  
};

#endif
