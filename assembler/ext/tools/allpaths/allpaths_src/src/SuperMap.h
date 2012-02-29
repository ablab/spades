// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef SUPER_MAP_H
#define SUPER_MAP_H

#include "AnAssemblyClass.h"
#include "Basevector.h"
#include "Vec.h"



/*
 * class super_map
 *
 * One supercontig map. Beware: stop is the last base in a contig (and it
 * belongs to the contig). You must run Setup( ) after each change in the
 * supercontig (and at the beginning if you used the empty constructor).
 *
 * Some method appear in two forms, e.g. GetContigLength( ) and
 * GetContigLengthPos( ). The difference is that the first method is
 * passed the contig id; meanwhile the second the initial position of
 * the contig in the supercontig.
 */
class super_map {
  
public:
  
  super_map( ) { }

  super_map( const super &new_super, const vec<int> &contig_lengths );
  
  super_map( const super &new_super, const vecbasevector &contigs );
  
  // You must run Setup( ) if you change the supercontig.
  void Setup( const super &new_super, const vec<int> &contig_lengths );
  
  // You must run Setup( ) if you change the supercontig.
  void Setup( const super &new_super, const vecbasevector &contigs );

  // If contig_id belongs to the supercontig.
  bool InSupercontig( int contig_id ) const;

  int GetNumberContigs() const;

  int GetSupercontigLength( ) const;

  int GetUngappedSupercontigLength( ) const;

  int GetContigLength( int contig_id ) const;
 
  int GetContigLengthPos( int contig_pos ) const;

  vec<int> GetContigLengths( ) const;

  int GetPositionInSupercontig( int contig_id ) const;
  
  int GetContigId( int pos_in_supercontig ) const;

  vec<int> GetContigIds( ) const;
  
  int GetGap( int after_pos_in_supercontig ) const;

  vec<int> GetGaps( ) const;
  
  int GetStartOnSupercontig( int contig_id ) const;

  int GetStartOnSupercontigPos( int contig_pos ) const;

  int GetUngappedStartOnSupercontig( int contig_id ) const;

  int GetUngappedStartOnSupercontigPos( int contig_pos ) const;

  int GetStopOnSupercontig( int contig_id ) const;

  int GetStopOnSupercontigPos( int contig_pos ) const;

  int GetUngappedStopOnSupercontig( int contig_id ) const;

  int GetUngappedStopOnSupercontigPos( int contig_pos ) const;

  // Sort contigs by start on supercontig.
  void SortContigsByStart( vec<int> &contig_ids ) const;
  
  // Sort contigs by minimizing negative gaps.
  void SortContigsByMaxGap( vec<int> &contig_ids ) const;

  
private:

  void Setup( const super &new_super,
	      const vec<int> *contig_lengths,
	      const vecbasevector *contigs );
  
  // Return the position in ids_ of contig_id.
  int FindPos( int contig_id ) const;

  // clear all (really everything).
  void Clear( );
  
  // reserve memory for all vectors.
  void Reserve( unsigned int size );

  void AddContig( int contig_id, int start, int stop );
  
  
  /*
   * ids_: ids in the contig;
   * start_: start bases of contigs;
   * stop_: stop bases of contigs;
   * contig_to_pos_: map (contig id -> position in ids_).
   */
  vec<int> ids_;
  vec<int> start_;
  vec<int> stop_;
  map<int, int> contig_to_pos_;
  
};



/*
 * order_contig_pos_Start
 * ordering functor
 *
 * Order contigs by start on supercontig. Notice that left and right refer
 * to the original positions in the supercontig, and not to their id.
 */
struct order_contig_pos_Start
  : public binary_function<const int&, const int&, bool>
{
private:
  const super_map *super_;
  
public:
  order_contig_pos_Start( const super_map *super ) :
    super_( super ) {}
  
  bool operator() ( const int &left, const int &right ) {
    int left_start = super_->GetStartOnSupercontigPos( left );
    int right_start = super_->GetStartOnSupercontigPos( right );
    
    return ( left_start < right_start );
  }

};



/*
 * order_contig_pos_MaxGap
 * ordering functor
 *
 * Order contigs by minimizing negative gaps. Notice that left and
 * right refer to the original positions in the supercontig, and not
 * to their id.
 */
struct order_contig_pos_MaxGap
  : public binary_function<const int&, const int&, bool>
{
private:
  const super_map *super_;

public:
  order_contig_pos_MaxGap( const super_map *super ) :
    super_( super ) {}
  
  bool operator() ( const int &left, const int &right ) {
    int left_begin = super_->GetStartOnSupercontigPos( left );
    int left_end = 1 + super_->GetStopOnSupercontigPos( left );
    int right_begin = super_->GetStartOnSupercontigPos( right );
    int right_end = 1 + super_->GetStopOnSupercontigPos( right );
    
    return ( right_begin - left_end > left_begin - right_end );
  }

};



#endif
