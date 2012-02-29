///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef ANNOTATEDCONTIG
#define ANNOTATEDCONTIG

#include <map>
#include "Alignment.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "system/Types.h"
#include "Vec.h"
#include "layout/ContigActualloc.h"
#include "Misc.h"

// A semiannotation of a contig gives a crude alignment of all or part of it with
// another "standard" contig.  The crude alignment data contains starting positions,
// orientation, and alignment length, but not an alignment itself.

// For both semiannotation's and annotation's (see below), if Rc is set,
// "This" is reversed relative to "Standard", and not the other way around.

class semiannotation {
  
 public:

  //
  // Constructors
  //
  semiannotation() :
    std_id_        ( -1    ),
    rc_            ( False ),
    start_on_this_ ( -1    ),
    start_on_std_  ( -1    ),
    implied_leftmost_on_std_( -1 ),
    length_        ( -1    ), 
    weight_        ( -1    ) {}
  

  semiannotation( int std_id, 
		  Bool rc, 
		  int start_on_this, 
		  int start_on_std, 
		  int implied_leftmost_on_std,
		  int length,
		  int weight ) :

    std_id_        ( std_id        ),
    rc_            ( rc            ),
    start_on_this_ ( start_on_this ),
    start_on_std_  ( start_on_std  ),
    implied_leftmost_on_std_ ( implied_leftmost_on_std ),
    length_        ( length        ),
    weight_        ( weight        ) {}
  //-----------------------------------------------------


  
  //
  // Const accessors
  //
  int StdID()             const { return std_id_;                  } // ID of the known genome it goes to (usually 0).
  Bool RC()               const { return rc_;                      } // Orientation of overlap to known genome
  
  int StartOnThis()       const { return start_on_this_;           } // Start of aligning segment, on this contig
  int StartOnStandard()   const { return start_on_std_;            } // Start of aligning segment, on the genome
  int Length()            const { return length_;                  } // Length of aligning segment
  int Weight()            const { return weight_;                  } // Number of reads "supporting" or "voting for" this aligning segment

  int ShiftFromStandard() const { return implied_leftmost_on_std_; } // Shift from standard (redundant, since we have start_on_std_).

  void Print( ostream &o ) const;                                    // Prints the semiannotation information

  const basevector& Standard( int i ) const { return (*standards_)[i]; } // The ith segment of the genome (if initialized).


  //--------------------------------------------------------------------

  int AddToWeight( int k ) { return weight_ += k; } // Increments the weight of this semiannotation. I.e. the number of reads supporting it.

  static void SetStandards( const vec< basevector > &standards ) { standards_ = &standards; } // Sets the known genome segments.
  

  //
  // Input/Output. If a semiannotation is printed with <<, then it is loadable with >>.
  // 
  friend ostream& operator<<( ostream &o, const semiannotation &a );
  friend istream& operator>>( istream &i,       semiannotation &a );
  //----------------------------------------------------------------

 private:
  
  int std_id_;                        // ID of known genome segment. Usually 0.
  Bool rc_;                           // Orientation of alignment to known genome
  int start_on_this_, start_on_std_;  // Start of alignment on this contig, and on the known genome.
  int implied_leftmost_on_std_;       /* A somewhat redundant, but convenient value.
					 If the semiannotation is fw, then implied start on standard
					 is start_on_std_ - start_on_this_;
					 But if the semiannotation is rc, and let L be the length of
					 the semiannotated contig, then the implied start on standard
					 is start_on_std_ - L + start_on_this + length_ */
  int length_;                        // Length of the alignment (semiannotation).
  int weight_;                        // Number of reads supporting this particular semiannotation
  static vec< basevector > const *standards_;

};

typedef vec< semiannotation >::iterator semiannot_itr;

struct semiannotation_loccomp {
  bool operator() ( const semiannotation &ca1, 
		    const semiannotation &ca2 ) 
       const ;
};


//
// An annotation of a contig gives an alignment of all or part of it with
// another "standard" contig.
//

class annotation {

 public:

  annotation() :
    std_id_    ( -1 ),
    rc_        ( True ) {}

  annotation( vec< basevector > &standards ) :
    std_id_    ( -1 ),
    rc_        ( True ) {}


  annotation( int std_id,
	      Bool rc,
	      vec< basevector > &standards ) :
    std_id_    ( std_id ),
    rc_        ( rc ) {}


  int StdID() const { return std_id_; }
  Bool Rc()   const { return rc_; }
  
  alignment&  Align()          { return align_; }
  const basevector& Standard( int i ) { return (*standards_)[i]; }

  static void SetStandards( const vec< basevector > &standards ) { standards_ = &standards; }
  
  friend ostream& operator<<( ostream &o, const annotation &a );
  friend istream& operator>>( istream &i, annotation &a );
 
 private:
  
  int std_id_;
  Bool rc_;
  alignment align_;
  static vec< basevector > const *standards_;

};

typedef vec< annotation >::iterator annot_itr;

//
// Merged contigs consist of arachne contigs. Class arachne_contig_placement
// holds this information. For example if a merged contig X of length 9000 consists
// of arachne contigs A, B of lengths 5000, with 1000 overlap between them, then
// this would give rise to two arachne_contig_placement objects:
// ( A, 0, 5000 ), and ( B, 4000, 9000 ).
//
class arachne_contig_placement {
 public:

  //
  // Constructors
  //
  arachne_contig_placement() :
    arachne_id_ ( "" ) {}

  arachne_contig_placement( String arachne_id,
			    int start,
			    int stop ) :
    arachne_id_ ( arachne_id ),
    start_      ( start      ),
    stop_       ( stop       ) {}
  //-----------------------------------------


  //
  // Const accessors
  //
  int ID    () const;                                        // Returns the Arachne layout ID
  int Start () const { return start_;                      } // Start of Arachne layout contig in final contig
  int Stop  () const { return stop_;                       } // Stop of Arachne layout contig in final contig
  int IsGap () const { return arachne_id_.Contains( "." ); } // Is this Arachne layout contig a gap contig?
  int Length() const { return stop_ - start_;              } // Length of Arachne layout contig
  
  String FullId () const { return arachne_id_; }
  //--------------------------------------------------------


  void Print( ostream &o ) const;
  
  friend ostream& operator<<( ostream &o, const arachne_contig_placement &a );
  friend istream& operator>>( istream &i, arachne_contig_placement &a );

 private:
  String arachne_id_;
  int start_;
  int stop_;
};

typedef vec< arachne_contig_placement >::iterator acp_itr;


// An annotated_contig is a contig, together with a list of semiannotations or
// annotations, or perhaps both, but the intention is that first one gets the
// semiannotations, and then converts them into annotations.

const int fudge_per_nucleotides = 100;
const int constant_fudge        = 100;
const int max_fudge             = 100; // Was infinity, until January 9, 2001.

class annotated_contig {

 public:
  //
  // Constructors
  // 
  annotated_contig() : 
    semiannotated_  ( False ),
    annotated_      ( False ) {}

  annotated_contig( int id,
		    const basevector                      &bases,
		    const qualvector                      &quals,
		    Bool                                  semiannotated,
		    Bool                                  annotated,
		    const vec< semiannotation >           &semiannotations,
		    const vec< annotation >               &annotations,
		    const vec< arachne_contig_placement > &arachne_ctgs ) :
    id_                     ( id              ),
    bases_                  ( bases           ),
    quals_                  ( quals           ),
    semiannotated_          ( semiannotated   ),
    annotated_              ( annotated       ),
    semiannotations_        ( semiannotations ),
    annotations_            ( annotations     ),
    arachne_ctg_placements_ ( arachne_ctgs    ) {}
  //------------------------------------------------------------------------


  //
  // Const accessors
  //
  Bool Semiannotated() const { return semiannotated_;                    } // Is this contig semiannotated (fragments of it matching the known genome)
  Bool Annotated()     const { return annotated_;                        } // Is it annotated (located in one or more potential places on the genome)

  int NumSemiannots()  const { return semiannotations_.size();           } // Number of different semiannotations
  int NumAnnots()      const { return annotations_.size();               } // Number of different annotations
  int NumACPs()        const { return arachne_ctg_placements_.size();    } // Number of different Arachne contigs comprising this final contig
  int Length()         const { return quals_.size();                     } // Number of bases in this final contig
  int ID()             const { return id_;                               } // ID of this final contig
  const basevector& Bases() const { return bases_; } // Basevector of this contig
  const qualvector& Quals() const { return quals_; } // Quality scores of this contig

void Print( ostream &out ) const;

  //---------------------------------------------------------------------

  // 
  // Accessors that are not const (yet)
  //
  const arachne_contig& ArachneContig( int i ) { return arachne_contigs_->find( i )->second; } // The ith Arachne layout contig in this final contig
  
  unsigned char  Base( int i )                 { return bases_[ i ];                         } // The ith base
  unsigned char& Qual( int i )                 { return quals_[ i ];                         } // The ith quality score
  const basevector& Standard( int i ) { return (*standards_)[i]; }
  //---------------------------------------------------------------------

  //
  // Setting values
  //
  Bool SetSemiannotated( Bool s ) { return ( semiannotated_ = s ); }
  Bool SetAnnotated    ( Bool a ) { return ( annotated_     = a ); }
  void SetBases( const basevector& b ) { bases_ = b; }
  void SetQuals( const qualvector& q ) { quals_ = q; }

  void SetID( int id ) { id_ = id; }
  
  static void SetStandards     ( const vec< basevector >          &standards       ) { standards_       = &standards;       }
  static void SetArachneContigs( const map< int, arachne_contig > &arachne_contigs ) { arachne_contigs_ = &arachne_contigs; }
  //---------------------------------------------------------------

  
  //
  // This function does all the work.
  // Starting from the arachne contig placements (i.e. the composition of this
  //   contig with Arachne layout contigs) and the actual locations (hypothesized
  //   during the Arachne layout phase) of the Arachne layout contigs,
  //   it computes all the plausible semiannotations of this final contig.
  // That is, a lot of bookkeeping.
  //
  void ComputeSemiannotations();
  //-----------------------------------------------------------------------------


  //
  // Input/Output operators that have the property that
  //   if the object is printed with <<, it can be retrieved with >>.
  //
  friend ostream& operator<<( ostream &o, const annotated_contig &a );
  friend istream& operator>>( istream &i,       annotated_contig &a );
  //------------------------------------------------------------------

  //
  // Iterators through the semiannotations, anotations, and arachne_contig_placements.
  //
  semiannot_itr FirstSemiannot() { return semiannotations_.begin();           }
  semiannot_itr EndSemiannot()   { return semiannotations_.end();             }
  
  annot_itr     FirstAnnot()     { return annotations_.begin();               }
  annot_itr     EndAnnot()       { return annotations_.end();                 }

  acp_itr       FirstACP()       { return arachne_ctg_placements_.begin();    }
  acp_itr       EndACP()         { return arachne_ctg_placements_.end();      }
  //---------------------------------------------------------------------------


  static vec< basevector >          const *standards_;
  static map< int, arachne_contig > const *arachne_contigs_;

 private:
  int                             id_;                      // Final contig id                     
  basevector                      bases_;                   // Bases of final contig
  qualvector                      quals_;                   // Qualitie scores of final contig
  Bool                            semiannotated_;           // Are the semiannotations computed?
  Bool                            annotated_;               // Are the annotations computed?
  vec< semiannotation >           semiannotations_;         // Semiannotations vector
  vec< annotation >               annotations_;             // Annotations vector
  vec< arachne_contig_placement > arachne_ctg_placements_;  // Vector of arachne contig placements

};

typedef vec< annotated_contig    >::iterator annotated_contig_itr;
typedef map< int, arachne_contig >::iterator arachne_ctg_map_itr;
typedef map< int, arachne_contig >::iterator const_arachne_ctg_map_itr;


//
// Same as annotated_contig, but this time it is a supercontig.
//
class annotated_supercontig {

 public:
  //
  // Constructors
  //
  annotated_supercontig() {}

  annotated_supercontig( const vec< annotated_contig > &contigs,
			 const vec< int >               gaps,
			 const vec< int >               gap_sds 
			 ) :
    contigs_   ( contigs ),
    gaps_      ( gaps    ),
    gap_sds_   ( gap_sds ) {
    
    Assert( contigs_.size() == gaps_.size() + 1 &&
	    gap_sds_.size() == gaps_.size() );
  }
  //--------------------------------------------------------------


  //
  // Const accessors
  //
  int NumContigs() const { return contigs_.size(); }    // Returns the number of final contigs in the supercontig
  int Length()     const;                               // Returns the sum of all contig lengths -- may need to modify this

  void Print( ostream &out ) const;
  //------------------------------------------------

  //
  // Non-const accessors
  //
  annotated_contig& Contig( int i )   { return contigs_[i]; } // Returns the ith annotated contig

  int& Gap           ( int i )        { return gaps_[i];    } // Returns the ith gap length
  int& GapStandardDev( int i )        { return gap_sds_[i]; } // Returns the ith gap length standard dev
  //---------------------------------------------------------


  //
  // Main function that does all the computations.
  // It computes the semiannotations. See also the similarly-named
  //   function in annotated_contig.
  //
  void ComputeSemiannotations();  
  //----------------------------


  void SetGap ( int i, 
		int g ) 
    { gaps_[i] = g; } // Sets the ith gap to have length g


  static void SetStandards( const vec< basevector > &standards ) { standards_ = &standards; }

  const basevector& Standard( int i ) { return (*standards_)[i]; }

  //
  // Iteration through the annotated_contigs
  //
  annotated_contig_itr FirstContig() { return contigs_.begin(); }
  annotated_contig_itr EndContig()   { return contigs_.end();   }
  //-------------------------------------------------------------


  // 
  // Input/Output operators that have the property that
  //   if the object is printed with <<, it can be retrieved with >>.
  //
  friend ostream& operator<<( ostream &o, const annotated_supercontig &a );
  friend istream& operator>>( istream &i,       annotated_supercontig &a );
  //-----------------------------------------------------------------------


 private:
  
  vec< annotated_contig > contigs_; // The list of annotated_contig objects 
  vec< int >              gaps_;    // The list of interleaving gap lengths
  vec< int >              gap_sds_; // The list of std deviations on interleaving gap lengths

  static vec< basevector > const *standards_;
};


//
// Compares the two supercontigs by comparing their lengths.
// The smallest one is the longest.
//
bool annot_sc_lengthcomp ( const annotated_supercontig &a,
			   const annotated_supercontig &b );


typedef vec< annotated_supercontig >::iterator annotated_supercontig_itr;


//
// Contains all the annotated_supercontig objects of the assembly.
// Also, it contains all the arachne contigs, and all the standard
//   (i.e. known) contigs, which are usually just one known genome
//   from which the simulated reads were taken.
//
class annotated_final_answer {

 public:
  annotated_final_answer() {}

  annotated_supercontig& SuperContig( int i ) { return supers_[i];    }
  basevector&            Standard   ( int i ) { return standards_[i]; }

  void PushSupercontig   ( const annotated_supercontig &s ) { supers_.push_back       ( s ); }
  void PushStandardContig( const basevector            &b ) { standards_.push_back    ( b ); }
  void PushArachneContig ( const arachne_contig        &a ) { 
    arachne_contigs_.insert ( pair< int, arachne_contig >( a.ID(), a ) );
  }
  arachne_contig& ArachneContig( int i );

  void ComputeSemiannotations();  
  
  vec< annotated_supercontig >&       RevealSupercontigs()    { return supers_;          }
  const vec< basevector >&            RevealStandardContigs() { return standards_;       }
  const map< int, arachne_contig >&   RevealArachneContigs()  { return arachne_contigs_; }

  annotated_supercontig_itr FirstSc() { return supers_.begin(); }
  annotated_supercontig_itr EndSc()   { return supers_.end();   }
  
  friend ostream& operator<<( ostream &o, const annotated_final_answer &a );
  friend istream& operator>>( istream &i, annotated_final_answer &a );

 private:
  
  map< int, arachne_contig >   arachne_contigs_;
  vec< annotated_supercontig > supers_;
  vec< basevector >            standards_;

};


ostream& operator<<( ostream &o, const semiannotation           &a );
ostream& operator<<( ostream &o, const annotation               &a );
ostream& operator<<( ostream &o, const arachne_contig_placement &a );
ostream& operator<<( ostream &o, const annotated_contig         &a );
ostream& operator<<( ostream &o, const annotated_supercontig    &a );
ostream& operator<<( ostream &o, const annotated_final_answer   &a );

istream& operator>>( istream &i,       semiannotation           &a );
istream& operator>>( istream &i,       annotation               &a );
istream& operator>>( istream &i,       arachne_contig_placement &a );
istream& operator>>( istream &i,       annotated_contig         &a );
istream& operator>>( istream &i,       annotated_supercontig    &a );
istream& operator>>( istream &i,       annotated_final_answer   &a );


#endif
