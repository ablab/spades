/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PLACEMENT_H
#define PLACEMENT_H

#include "CommonSemanticTypes.h"
#include "SeqInterval.h"
#include "math/HoInterval.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

/**
   Class: placement

   Represents the placement of a unipath or a read on the <reference genome>.

   *NOTE*: Not to be confused with same-named class in <VAALTools.h> in polymorphism/ dir!
*/
class placement {

     public:

     placement( ) { }

     placement( genome_part_id_t genome_id, genome_part_pos_t pos, genome_part_pos_t Pos, Bool rc ) :
          genome_id_(genome_id), pos_(pos), Pos_(Pos), rc_(rc) { }

     genome_part_pos_t pos( ) const { return pos_; }
     genome_part_pos_t Pos( ) const { return Pos_; }
     int length( ) const { return Pos_ - pos_; }
     ho_interval Interval( ) const { return ho_interval( pos_, Pos_ ); }
     genome_part_id_t GenomeId( ) const { return genome_id_; }
     seq_interval SeqInterval() const { return seq_interval( -1, genome_id_, pos_, Pos_ ); }
     Bool Rc( ) const { return rc_; }
     Bool Fw( ) const { return !rc_; }
     orient_t Orient() const { return rc_ ? ORIENT_RC : ORIENT_FW; }

     friend Bool operator<( const placement& p1, const placement& p2 )
     {    if ( p1.genome_id_ < p2.genome_id_ ) return True;
          if ( p1.genome_id_ > p2.genome_id_ ) return False;
          if ( p1.pos_ < p2.pos_ ) return True;
          return False;    }

     friend Bool operator==( const placement& p1, const placement& p2 )
     {    return p1.genome_id_ == p2.genome_id_ && p1.pos_ == p2.pos_
               && p1.Pos_ == p2.Pos_ && p1.rc_ == p2.rc_;    }

     friend ostream& operator<<( ostream& out, const placement& p )
     {    return out << p.genome_id_ << "." << p.pos_ << "-" << p.Pos_
               << ( p.rc_ ? " (rc)" : " (fw)" );    }

     private:

     genome_part_id_t genome_id_;
     genome_part_pos_t pos_, Pos_;
     Bool rc_;
  int nplaces_; // not used

};

TRIVIALLY_SERIALIZABLE(placement);

typedef SerfVec<placement> PlacementVec;
typedef MasterVec<PlacementVec> VecPlacementVec;

#endif
