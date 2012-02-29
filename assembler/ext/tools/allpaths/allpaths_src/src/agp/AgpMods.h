// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef AGPMODS_H
#define AGPMODS_H

#include "String.h"
#include "Vec.h"

// The class contig_range provides begin and end positions on a
// contig.  If the end position is before the begin position, the
// range is considered reversed.  If a length is provided, the range
// can be switched back and forth.

class contig_range
{
 private:
  int begin_pos_;
  int end_pos_;

 public:

  int BeginPos() const { return ( this->IsReversed() ? end_pos_ : begin_pos_ ); }
  int EndPos() const { return ( this->IsReversed() ? begin_pos_ : end_pos_ ); }

  bool IsReversed() const { return ( begin_pos_ > end_pos_ ); }
  void Reverse()
    {
      swap( begin_pos_, end_pos_ );
    }

  int RcBeginPos( int contig_length ) const
    {
      Assert( this->IsReversed() );
      return contig_length - begin_pos_ - 1;
    }

  int RcEndPos( int contig_length ) const
    {
      Assert( this->IsReversed() );
      return contig_length - end_pos_ - 1;
    }

  void SetBeginPos( const int begin_pos ) { begin_pos_ = begin_pos; }
  void SetEndPos( const int end_pos ) { end_pos_ = end_pos; }

  bool Overlap( const contig_range& other ) const
    {
      return ( this->EndPos() >= other.BeginPos() &&
	       this->BeginPos() <= other.EndPos() );
    }

  friend ostream& operator<< ( ostream& out, const contig_range& anchor );
};



// The class insertion_anchor contains the id of the anchoring contig
// and the ranges on the anchoring contig and the replacement contig
// that match.  The ranges have the same orientation if either both
// are reversed or both are not.  The anchor can be reversed
// (reversing the contained ranges) if the length of the anchor and
// the length of the replacement are provided.

class insertion_anchor
{
 private:
  int anchor_id_;
  contig_range anchor_range_;
  contig_range replacement_range_;

 public:
  int Id() const { return anchor_id_; }
  contig_range AnchorRange() const { return anchor_range_; }
  contig_range ReplacementRange() const { return replacement_range_; }
  
  void SetId( const int id )
    { anchor_id_ = id; }
  void SetAnchorRange( const contig_range& anchor_range ) 
    { anchor_range_ = anchor_range; }
  void SetReplacementRange( const contig_range& replacement_range ) 
    { replacement_range_ = replacement_range; }

  bool SameOrientation() const
    {
      return ( ( anchor_range_.IsReversed() && replacement_range_.IsReversed() ) ||
	       ( ! anchor_range_.IsReversed() && ! replacement_range_.IsReversed() ) );
    }

  bool AnchorIsReversed() const
    { return anchor_range_.IsReversed(); }

  bool ShouldRemoveAnchor( int anchor_length ) const
    {
      return ( anchor_range_.BeginPos() == 0 && 
	       anchor_range_.EndPos() == anchor_length - 1 );
    }

  void Reverse()
    {
      anchor_range_.Reverse();
      replacement_range_.Reverse();
    }

  friend istream& operator>> ( istream& in, insertion_anchor& anchor );
  friend ostream& operator<< ( ostream& out, const insertion_anchor& anchor );
};


// The class insertion specifies the name of the replacement contig
// and the two anchors in the assembly for that contig.  The entire
// insertion can be reversed if the lengths of the two anchors and the
// length of the replacement contig are provided.

class insertion
{
 private:
  String id_;
  insertion_anchor first_anchor_;
  insertion_anchor last_anchor_;

 public:
  String Id() const { return id_; }
  insertion_anchor FirstAnchor() const { return first_anchor_; }
  insertion_anchor LastAnchor() const { return last_anchor_; }
  
  void SetId( const String& id ) { id_ = id; }
  void SetFirstAnchor( const insertion_anchor& anchor ) { first_anchor_ = anchor; }
  void SetLastAnchor( const insertion_anchor& anchor ) { last_anchor_ = anchor; }

  void Reverse()
    {
      first_anchor_.Reverse();
      last_anchor_.Reverse();
      swap( first_anchor_, last_anchor_ );
    }

  bool IsFirstAnchor( const int contig_id ) const
    { return ( contig_id == first_anchor_.Id() ); }
  
  bool IsLastAnchor( const int contig_id ) const
    { return ( contig_id == last_anchor_.Id() ); }

  friend istream& operator>> ( istream& in, insertion& the_insertion );
  friend ostream& operator<< ( ostream& out, const insertion& the_insertion );
};


// The class contig_trim stores the id, the amount and the end of a
// contig to be trimmed.

class contig_trim
{
 private:
  int id_;
  int amount_;
  bool front_;

 public:
  contig_trim()
    : id_(-1) { }
  contig_trim( const int id, const int amount, const bool trim_front )
    : id_(id), amount_(amount), front_(trim_front) { }

  int Id() const { return id_; }
  int Amount() const { return amount_; }
  bool TrimFront() const { return front_; }
  bool TrimBack() const { return ! front_; }

  void SetId( const int id ) { id_ = id; }
  void SetAmount( const int amount ) { amount_ = amount; }
  void SetTrimFront( const bool trim_front ) { front_ = trim_front; }
  void SetTrimBack( const bool trim_back ) { front_ = ! trim_back; }

  friend istream& operator>> ( istream& in, contig_trim& trim );
  friend ostream& operator<< ( ostream& out, const contig_trim& trim );
};


// The class agp_mods allows input and output of formatted AGP modifications.

class agp_mods
{
 public:
  bool SuperShouldBeRemoved( int super_id ) const;
  const vec<int>& GetSupersToRemove() const { return supers_to_remove_; }

  bool ContigShouldBeUnplaced( int contig_id ) const;
  const vec<int>& GetContigsToUnplace() const { return contigs_to_unplace_; }

  const vec<insertion>& GetInsertions() const { return insertions_; }
  
  const vec<contig_trim>& GetTrims() const { return trims_; }

  void Read( const String& filename );
  void Write( const String& filename ) const;

 private:
  vec<int> supers_to_remove_;
  vec<int> contigs_to_unplace_;
  vec<insertion> insertions_;
  vec<contig_trim> trims_;
};

#endif  






