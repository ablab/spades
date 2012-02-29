// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef VEC_INSERTS_H
#define VEC_INSERTS_H

#include "AnAssemblyClass.h"
#include "InsertEnds.h"



/*
 * class vec_inserts
 *
 * An organized vector-of-inserts class. The minimum amount of data
 * corresponds to the inserts in one contig.
 */
class vec_inserts {
  
public:
  
  vec_inserts( const assembly *the_assembly );
  
  void Clear( );

  void AddAllContigs( );
  
  void AddSupercontigs( const vec<int> &super_ids );
  
  void AddContigs( const vec<int> &contig_ids );
  
  unsigned int size( ) const { return inserts_.size( ); }
  
  const insert_ends& operator[] ( int ii ) const { return inserts_[ii]; }

  void FindGapLinks( int super_id, int pos, vec<int> &insert_ids ) const;
  
  void FindSpotLinks( int contig_id, int pivot, vec<int> &insert_ids ) const;

  
private:
  
  const assembly *the_assembly_;  // the assembly (raw data).
  
  vec<int> contig_ids_;           // contig ids (sorted).
  vec<insert_ends> inserts_;      // all inserts (sorted).
  
};



#endif
