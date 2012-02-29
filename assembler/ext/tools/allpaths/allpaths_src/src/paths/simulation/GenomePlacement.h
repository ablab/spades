/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_SIMULATION_GENOMEPLACEMENT_H
#define PATHS_SIMULATION_GENOMEPLACEMENT_H

#include "system/Types.h"
#include "SemanticTypes.h"

#include <iostream>
#include <iomanip>

/*
   Class: genome_placement
   
   The placement (gap-free alignment) of a <base vector>, such as a <read> or a
   <unipath>, onto the <reference genome> . In the documentation of this class we'll usually call the thing
   placed on the genome a read, but it can be any basevector.

   One record of <.placements> files.  Produced by <PathsToLocs>.
*/
class genome_placement {

     public:

  genome_placement( ):
    read_id_( -1 ), kmer_count_( -1 ), genome_id_( -1 ), pos_( -1 ),
    Pos_( -1 ), rc_( False ), nplaces_( -1 ) { }

     /*
        Constructor: genome_placement constructor

	Create a genome placement.

	Parameters:

	   read_id - the id of the basevector placed on the genome -- that is,
	        the index of this basevector in some <vecbasevector>.
	   kmer_count - the number of kmers in this read or path.  easilly derivable
	      from the length of the read or path represented by read_id, but copied here
	      so we have a self-contained structure.
	   genome_id - the <genome id> for the <genome part> onto which the read is placed.
	   pos - index of the first base of the genome part onto which the read is placed.
	   Pos - index of the last base of the genome part onto which the read is placed.
	     Note that this index is inclusive: it is the index of the last base of the placement,
	     not of the base after the last base.
	   rc - whether the read itself aligns to the genome part, or the read's reverse complement.
	   nplaces - to how many places in the genome does this read or unipath align?  makes sense
	      only in the context of a _vector_ of genome_placements -- all genome_placements with
	      a given read_id have their nplaces set to the total number of genome_placements with
	      that read_id in that vector.
     */
     genome_placement( int read_id, int kmer_count, int genome_id, int pos, int Pos,
          Bool rc, int nplaces ) : read_id_(read_id), kmer_count_(kmer_count),
          genome_id_(genome_id), pos_(pos), Pos_(Pos), rc_(rc), nplaces_(nplaces) { }

  Bool IsValid() const {
    return read_id_ >= 0  &&  kmer_count_ >= 0  &&  genome_id_ >= 0 && pos_ >= 0
      && Pos_ >= 0 && nplaces_ >= 0;
  }

     friend bool operator<( const genome_placement& p1, const genome_placement& p2 )
     {    if ( p1.genome_id_ < p2.genome_id_ ) return True;
          if ( p1.genome_id_ > p2.genome_id_ ) return False;
          if ( p1.pos_ < p2.pos_ ) return True;
          return False;    }

     struct OrderByRead 
       : public binary_function<genome_placement,genome_placement,bool>
     {
       bool operator() ( const genome_placement& lhs, const genome_placement& rhs ) const {
         return ( lhs.GetReadId() < rhs.GetReadId() );
       }
     };

     friend bool operator==( const genome_placement& p1, const genome_placement& p2 )
     {    return ( p1.read_id_ == p2.read_id_ &&
                   p1.kmer_count_ == p2.kmer_count_ &&
                   p1.genome_id_ == p2.genome_id_ &&
                   p1.pos_ == p2.pos_ &&
                   p1.Pos_ == p2.Pos_ &&
                   p1.rc_ == p2.rc_ &&
                   p1.nplaces_ == p2.nplaces_ );
     }

     friend ostream& operator<<( ostream& out, const genome_placement& p )
     {    return out << "path " << setiosflags(ios::fixed) << setw(6) << p.read_id_
                     << " (" << setw(4) << p.kmer_count_ << " kmers)"
                     << " --> " << p.genome_id_ << "." << p.pos_ << "-" << p.Pos_
                     << ( p.rc_ ? " (rc)" : " (fw)" ) << " [" << p.nplaces_ << " places]"
                     << resetiosflags(ios::fixed) << "\n";    }

     // Method: GetReadId
     // Returns the id of the basevector which we're placing on the genome.  The id is the
     // index of this basevector in its <vecbasevector>, be that a list of reads or of unipaths
     // or of something else.
     int GetReadId() const { return read_id_; }

     // Method: GetKmerCount
     // Returns the number of kmers in the unipath that we're placing on the genome --
     // but isn't this the same information you can get from calling <GetReadId()>?
     // Probably, but this lets you get the length of the placed unipath without
     // pulling up the unipath file and looking up the unipath itself by its index in the file.
     int GetKmerCount() const { return kmer_count_; }

     // Method: GetGenomeId
     // Get the id of the <genome part> onto which we're placing our basevector.
     int GetGenomeId() const { return genome_id_; }

     // Method: GetStartOnGenome
     // Get the index at which the alignment of our basevector onto its genome part starts.
     int GetStartOnGenome() const { return pos_; }

     // Method: GetEndOnGenome
     // Get the index at which the alignment of our basevector onto its genome part ends.
     // Note that this gives the actual last position, not the position right after that.
     int GetEndOnGenome() const { return Pos_; }

  int StartOnTarget() const { return GetStartOnGenome(); }
  int EndOnTarget() const { return GetEndOnGenome(); }
  int QueryId() const { return read_id_; }
  int TargetId() const { return genome_id_; }

     // Method: IsRc
     // Tells whether it is actually the reverse complement of our basevector that aligns
     // to the genome at the specified point, not our basevector itself.
     bool IsRc() const { return rc_; }

     // Method: IsFw
     // Tells whether it is our basevector itself that aligns to the genome at the specified
     // point, rather than the reverse complement of our basevector.
     bool IsFw() const { return !rc_; }
     
     int GetCopyNumber() const { return nplaces_; }
     void SetCopyNumber( int n ) { nplaces_ = n; }

     private:

     int read_id_;
     int kmer_count_;
     int genome_id_;
     int pos_, Pos_;
     Bool rc_;
     int nplaces_;

};


// Semantic type: genome_placement_id_t
// Represents the id of a genome_placement
// in a vector of genome placements.
SemanticTypeStd( int, genome_placement_id_t );

#endif


// Synonyms: Various synonyms
//   genome placement - See <genome_placement>
//   placement - See <genome_placement>

