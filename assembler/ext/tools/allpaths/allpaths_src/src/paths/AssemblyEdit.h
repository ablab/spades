///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// An assembly_edit is a proposed change to a scaffold that replaces a given
// stretch of bases (the 'source') by another stretch of bases or list of 
// alternatives (the 'replacement').  There are two cases: (i) 'internal': either 
// the edit is entirely within a contig, or else (ii) 'gap-closer': it contains a 
// gap (in which case it must contain the entire gap, along with zero or more bases 
// from its flanks).

#ifndef ASSEMBLY_EDIT_H
#define ASSEMBLY_EDIT_H

#include "Basevector.h"
#include "CoreTools.h"
#include "efasta/EfastaTools.h"

class assembly_edit {

     public:

     enum edit_type { INTERNAL, GAP_CLOSER };

     assembly_edit( ) { }

     assembly_edit( const edit_type& etype, const int tig, const int start,
          const int stop, const vec<basevector>& reps )
          : etype_(etype), tig_(tig), start_(start), stop_(stop), reps_(reps)
     {    ForceAssert( etype == INTERNAL );
          ForceAssertLe( 0, tig );
          ForceAssertLe( 0, start );
          ForceAssertLe( start, stop );
          ForceAssertGe( reps.isize( ), 1 );    }

     assembly_edit( const edit_type& etype, const int tig1, const int start1,
          const int tig2, const int stop2, const vec<basevector>& reps )
          : etype_(etype), tig_(tig1), tign_(tig2), start_(start1), stop_(stop2), 
          reps_(reps)
     {    ForceAssert( etype == GAP_CLOSER );
          ForceAssertLe( 0, tig1 );
          ForceAssertLe( 0, tig2 );
          ForceAssertLe( 0, start1 );
          ForceAssertLe( 0, stop2 );
          ForceAssertGe( reps.isize( ), 1 );    }

     Bool Internal( ) const { return etype_ == INTERNAL; }
     Bool GapCloser( ) const { return etype_ == GAP_CLOSER; }

     int Tig( ) const 
     {    Assert( Internal( ) );
          return tig_;    }

     int Tig1( ) const
     {    Assert( GapCloser( ) );
          return tig_;    }

     int Tig2( ) const
     {    Assert( GapCloser( ) );
          return tign_;    }

     int Start( ) const
     {    Assert( Internal( ) );
          return start_;    }

     int Start1( ) const 
     {    Assert( GapCloser( ) );
          return start_;    }
     void SetStart1( int start1 ) { start_ = start1; }

     int Stop( ) const
     {    Assert( Internal( ) );
          return stop_;    }

     int Stop2( ) const 
     {    Assert( GapCloser( ) );
          return stop_;    }
     void SetStop2( int stop2 ) { stop_ = stop2; }

     int Nreps( ) const { return reps_.size( ); }
     const basevector& Rep( int i ) const { return reps_[i]; }
     basevector& Rep( int i ) { return reps_[i]; }
     const vec<basevector>& Reps( ) const { return reps_; }

     Bool PureDeletion( ) const
     {    return Internal( ) && Nreps( ) == 1 && Rep(0).size( ) == 0;    }

     friend Bool operator<( const assembly_edit& e1, const assembly_edit& e2 )
     {    Assert( e1.Internal( ) && e2.Internal( ) );
          if ( e1.Tig( ) < e2.Tig( ) ) return True;
          else if ( e1.Tig( ) > e2.Tig( ) ) return False;
          if ( e1.Start( ) < e2.Start( ) ) return True;
          else if ( e1.Start( ) > e2.Start( ) ) return False;
          if ( e1.Stop( ) < e2.Stop( ) ) return True;
          else if ( e1.Stop( ) > e2.Stop( ) ) return False;
          return e1.Reps( ) < e2.Reps( );    }

     friend ostream& operator<<( ostream& out, const assembly_edit& e )
     {    if ( e.Internal( ) )
          {    return out << "INTERNAL " << e.Tig( ) << "." << e.Start( ) << "-" 
                    << e.Stop( ) << " --> " << efasta( e.Reps( ) );    }
          else
          {    return out << "GAP_CLOSER " << e.Tig1( ) << "." << e.Start1( ) 
                    << " - " << e.Tig2( ) << "." << e.Stop2( )
                    << " --> " << efasta( e.Reps( ) );    }    }

     friend Bool Overlap( const assembly_edit& e1, const assembly_edit& e2 )
     {    if ( e1.Internal( ) && e2.Internal( ) )
          {    if ( e1.Tig( ) != e2.Tig( ) ) return False;
               return IntervalOverlap( e1.Start( ), e1.Stop( ),
                    e2.Start( ), e2.Stop( ) );    }
          else if ( e1.Internal( ) && e2.GapCloser( ) )
          {    if ( e1.Tig( ) == e2.Tig( ) ) return e1.Stop( ) > e2.Start1( );
               else if ( e1.Tig( ) == e2.Tig2( ) ) return e1.Start( ) < e2.Stop2( );
               else return False;    }
          else if ( e2.Internal( ) && e1.GapCloser( ) )
          {    if ( e2.Tig( ) == e1.Tig1( ) ) return e2.Stop( ) > e1.Start1( );
               else if ( e2.Tig( ) == e1.Tig2( ) ) return e2.Start( ) < e1.Stop2( );
               else return False;    }
          else
          {    if ( e1.Tig2( ) == e2.Tig1( ) ) return e1.Stop2( ) > e2.Start1( );
               if ( e2.Tig2( ) == e1.Tig1( ) ) return e2.Stop2( ) > e1.Start1( );
               return False;    }    }

     private:

     edit_type etype_;      // type of edit
     int tig_;              // contig id (the first contig in gap-closer case)
     int tign_;             // the second contig (gap-closer case only)
     int start_;            // start on contig (on first contig in gap-closer case)
     int stop_;             // stop on contig (on second contig in gap-closer case)
     vec<basevector> reps_; // the replacements

     public:

     size_t writeBinary(BinaryWriter& writer) const
     {    size_t len = 0;
          int et = etype_;
          len += writer.write(et);
          len += writer.write(tig_);
          len += writer.write(tign_);
          len += writer.write(start_);
          len += writer.write(stop_);
          len += writer.write(reps_);
          return len;    }

     void readBinary(BinaryReader& reader)
     {    int et;
          reader.read(&et);
          etype_ = (edit_type) et;
          reader.read(&tig_);
          reader.read(&tign_);
          reader.read(&start_);
          reader.read(&stop_);
          reader.read(&reps_);    }

      static size_t externalSizeof( ) { return 0; }
};

SELF_SERIALIZABLE(assembly_edit);

#endif
