///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef CONTIG_ACTUALLOC_H
#define CONTIG_ACTUALLOC_H

#include <map>
#include <set>

#include "Vec.h"
#include "system/Crash.h"
#include "Misc.h"


//
// Class contig_actualloc
//
// It contains a potential actual location on the genome for a contig,
//   in the case of simulated data (primarily) or whenever plausible
//   locations of the reads on the genome are given to the layout
//   phase.
// 
// More specifically, it contains the position and orientation of the
//   contig on the genome, or on the "standard" contig of id StdID().
//   It also contains the weight of this position/orientation, which
//   is basically the number of reads "voting" for it.
//

class contig_actualloc {

 public:

  //
  // Constructors/ Trivial Destructor
  //

  contig_actualloc() {}

  contig_actualloc( int weight, 
		    int start, 
		    Bool rc ) :
    weight_ ( weight ),
    start_  ( start ),
    rc_     ( rc ) {}

  ~contig_actualloc() {}

  //---------------------------

  // 
  // Const Accessors
  // 
  int Weight() const { return weight_; } // Number of reads "voting" for this actualloc
  int Start()  const { return start_;  } // Shift from standard
  int RC()     const { return rc_;     } // Orientation wrt standard
  int StdID()  const { return 0;       } // Not implemented yet
  //------------------------------------


  //
  // Const output function
  //
  void Print( ostream &o ) const;
  //-----------------------------

  //
  // Adds to the weight of the actualloc
  //
  int AddToWeight( int w ) { return weight_ += w; }

  friend ostream& operator<<( ostream &o, const contig_actualloc &a );
  friend istream& operator>>( istream &i, contig_actualloc &a );
  
 private:
  int weight_; // Number of reads voting for it
  int start_;  // Shift from known genome
  Bool rc_;    // Orientation in known genome
};
    
bool operator<( const contig_actualloc &a1, const contig_actualloc &a2 );

typedef set< contig_actualloc >::iterator         ctg_actloc_itr;


class arachne_contig {
 public:
  arachne_contig() :
    id_     ( -1 ),
    length_ ( -1 ) {}

  arachne_contig( int id, int length, int sc_id, int sc_pos, const set< contig_actualloc> &actuallocs ) :
    id_         ( id         ),
    length_     ( length     ),
    sc_id_      ( sc_id      ),
    sc_pos_     ( sc_pos     ),
    actuallocs_ ( actuallocs ) {}

  int ID     () const { return id_;     }
  int Length () const { return length_; }
  int SC_ID  () const { return sc_id_;  }
  int SC_Pos () const { return sc_pos_; }
  void Print( ostream &o ) const;

  friend ostream& operator<<( ostream &o, const arachne_contig &a );
  friend istream& operator>>( istream &i, arachne_contig &a );

  const ctg_actloc_itr FirstActloc() const { return actuallocs_.begin(); }
  const ctg_actloc_itr EndActloc()   const { return actuallocs_.end();   }
  

 private:
  int id_;
  int length_;
  int sc_id_;
  int sc_pos_;
  set< contig_actualloc > actuallocs_;
};

ostream& operator<<( ostream &o, const contig_actualloc &a );
ostream& operator<<( ostream &o, const arachne_contig   &a );

istream& operator>>( istream &i, contig_actualloc &a );
istream& operator>>( istream &i, arachne_contig   &a );

typedef map< int, arachne_contig > arachne_contig_itr;

#endif
