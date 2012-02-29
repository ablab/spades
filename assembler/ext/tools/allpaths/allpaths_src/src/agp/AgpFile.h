// Copyright (c) 2000, 2001 Whitehead Institute for Biomedical Research
// 

#ifndef AGPFILE_H
#define AGPFILE_H

#include "String.h"
#include "Vec.h"

#include "agp/AgpMods.h"

class agp_entry
{
 public:
  virtual void Print( ostream& out,
		      const String& chromosome,
		      longlong& base_position,
		      longlong& unit_position,
                      const String *base_name = 0 ) = 0;

  virtual ~agp_entry() { };
};

class agp_contig : public agp_entry
{
 public:
  enum contig_type { wgs_contig, 
                     finished_contig, 
                     draft_contig, 
                     blunt_joined_contig,
                     other_contig };

 private:
  String name_;
  int id_;
  int length_;
  int start_;
  int stop_;
  contig_type type_;
  bool is_rc_;

  char TypeChar_() const;

 public:
  agp_contig()
    : start_(0), stop_(-1) { }

  agp_contig( const String& name, 
	      const int length,
	      const int start, const int stop,
	      const bool is_rc )
    : name_(name), id_(-1), length_(length), 
    start_(start), stop_(stop), type_(finished_contig), is_rc_(is_rc) { }

  agp_contig( const int id,
	      const int length,
	      const int start, const int stop,
	      const bool is_rc )
    : id_(id), length_(length), 
    start_(start), stop_(stop), type_(wgs_contig), is_rc_(is_rc) { }

  String Name() const { return ( id_ < 0 ? name_ : String("contig_") + ToString(id_) ); }
  int Id() const { return id_; }
  int Length() const { return length_; }
  int Start() const { return start_; }
  int Stop() const { return stop_; }
  contig_type Type() const { return type_; }
  bool IsRC() const { return is_rc_; }

  int Extent() const { return stop_ - start_ + 1; }

  void SetName( const String& name ) { name_ = name; id_ = -1; }
  void SetId( const int id ) { id_ = id; }
  void SetLength( const int length ) { length_ = length; }
  void SetStart( const int start ) { start_ = start; }
  void SetStop( const int stop ) { stop_ = stop; }
  void SetType( const contig_type type ) { type_ = type; }
  void SetRC( const bool rc ) { is_rc_ = rc; }

  bool Check( int cg_length ) const;

  void Print( ostream& out,
	      const String& chromosome,
	      longlong& base_position,
	      longlong& unit_position,
              const String *base_name = 0 );
};

class agp_gap : public agp_entry
{
 private:
  String type_;
  int len_;
  bool is_bridged_;

 public:
  agp_gap()
    : len_(-1) { }
  agp_gap( const String& type, const int length, const bool is_bridged )
    : type_(type), len_(length), is_bridged_(is_bridged) { }

  String Type() const { return type_; }
  int Length() const { return len_; }
  bool IsBridged() const { return is_bridged_; }

  void SetType( const String& type ) { type_ = type; }
  void SetLength( const int length ) { len_ = length; }
  void SetBridged( const bool bridged ) { is_bridged_ = bridged; }

  void ExtendToAtLeast( const int minimum_length )
    {
      this->SetLength( max<int>( this->Length(), minimum_length ) );
    }

  void Print( ostream& out,
	      const String& chromosome,
	      longlong& base_position,
	      longlong& unit_position,
	      const String *base_name = 0 );
};

class agp_chromosome
{
 private:

  enum entry_type { contig_type, gap_type };
  
  typedef struct entry_pair {
    entry_type type;
    int index;
    
    entry_pair( entry_type the_type, int the_index )
      : type(the_type), index(the_index) { }
  } entry_pair;
  
  String name_;
  vec<agp_contig> contigs_;
  vec<agp_gap> gaps_;
  vec<entry_pair> entries_;
  
  
public:
  agp_chromosome() { }
  
  agp_chromosome( const String& name )
    : name_(name) { }

  String Name() const { return name_; }
  
  unsigned int GetEntriesCount( ) const { return entries_.size( ); }

  void GetEntry( int ii, agp_contig **contig, agp_gap **gap );

  void SetName( const String& name ) { name_ = name; }

  void AddSuper( const vec<int>& contig_ids, 
		 const vec<int>& contig_lens, 
		 const vec<int>& gaps,
		 const int minimum_gap,
		 const agp_mods& modifications,
		 const bool is_reversed );

  // Appends the given contig to the end of the chromosome.
  void AddContig( const agp_contig& contig );

  // If the last entry is a contig, appends the given gap to the end
  // of the chromosome.  If the last entry is a gap, two things can
  // happen.  If either or both the last entry and the given gap are
  // bridged and strict is true, the given gap is appended to the end
  // of the chromosome.  If neither is bridged or strict is false, the
  // existing gap is lengthened TO AT LEAST the length of the given
  // gap (and the gap is set to be unbridged).  This is to accomodate
  // the use of gaps of a certain size to indicate gaps of unknown
  // length.
  void AddGap( const agp_gap& gap, bool strict = true );

  bool AttemptInsertion( agp_contig& new_contig,
			 insertion ins,
			 const int minimum_gap_size );

  bool IsEmpty() const { return entries_.empty(); }

  // Load an agp_chromosome from file.
  void Load( const String &file_name );

  void Print( ostream& out,
	      const String *base_name = 0 );
};

#endif
