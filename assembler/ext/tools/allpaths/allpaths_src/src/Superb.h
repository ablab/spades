///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// superb (pronounced "super b") is a supercontig structure which carries the
// following information:
// 1. contig ids
// 2. contig lengths
// 3. gap estimates
// 4. standard deviation estimates for each gap estimate.
// Thus this class carries slightly more information than the existing class
// super (pronounced "super a").

// superb_by_tig: encapsulates supercontig information, on a contig-by-contig basis.
// (Note that the class definition is likely to be expanded.)

#ifndef SUPERB_H
#define SUPERB_H

#include "Charvector.h"
#include "CoreTools.h"

// TODO: potentially dangerous truncation of index throughout class
class superb {

     private:

     // Class gapb is intended for internal use by class superb, but could be made 
     // public if we think of a good reason to do so:

     class gapb {

          public:

          gapb( ) { }
          gapb( int g, int d ) : gap_(g), dev_(d) { }

          int Gap( ) const { return gap_; }
          int Dev( ) const { return dev_; }

          void SetGap( int g ) { gap_ = g; }
          void SetDev( int d ) { dev_ = d; }

          private:
          
          int gap_, dev_;

     };

     public:

     // Public member functions and friends:

     superb( ) { }

     int Ntigs( ) const { return mtig_.size( ); }
     int Ngaps( ) const { return mtig_.size( ) - 1; }

     int Tig( int i ) const { return mtig_[i]; }

     int Len( int i ) const { return len_[i]; }
     int Gap( int i ) const { return gap_[i].Gap( ); }
     int Dev( int i ) const { return gap_[i].Dev( ); }

     void Clear( )
     {    mtig_.clear();
          len_.clear();
          gap_.clear();    }

     void SetNtigs( int n ) {
       if ( n == 0 ) { Clear( ); return; }
       mtig_.resize( n );
       len_.resize( n );
       gap_.resize( n - 1 );
     }

     void SetTig( int i, int m ) { mtig_[i] = m; }
     void SetLen( int i, int l ) { len_[i] = l; }
     void SetGap( int i, int g ) { gap_[i].SetGap(g); }
     void SetDev( int i, int d ) { gap_[i].SetDev(d); }

     int TigIndex( int m ) const
     {    for ( int i = 0; i < Ntigs( ); i++ )
               if ( Tig(i) == m ) return i;
          ForceAssert( 0 == 1 );
          return -1;    }

     int TrueBegin( ) const {
       int pos = 0;
       int beg = 0;
       for ( int i = 0; i < Ntigs( ); i++ ) {
	 pos += Len(i);
	 if ( i < Ntigs( ) - 1 ) pos += Gap(i);
	 beg = beg < pos ? beg : pos;
       }
       return beg;
     }
  
     int TrueEnd( ) const {
       int sum = 0;
       int end = 0;
       for ( int i = 0; i < Ntigs( ); i++ ) {
	 sum += Len(i);
	 end = end > sum ? end : sum;
	 if ( i < Ntigs( ) - 1 ) sum += Gap(i);
       }
       return end;
     }
  
     int TrueLength( ) const {
       return ( this->TrueEnd( ) - this->TrueBegin( ) );
     }
  
     int ReducedLength( ) const
     {    int sum = 0;
          for ( int i = 0; i < Ntigs( ); i++ )
               sum += Len(i);
          return sum;    }

     int FullLength( ) const
     {    int sum = 0;
          for ( int i = 0; i < Ntigs( ); i++ )
          {    sum += Len(i);
               if ( i < Ntigs( ) - 1 ) sum += Gap(i);    }
          return sum;    }

     int FullLengthNN( ) const
     {    int sum = 0;
          for ( int i = 0; i < Ntigs( ); i++ )
          {    sum += Len(i);
               if ( i < Ntigs( ) - 1 && Gap(i) >= 0 ) sum += Gap(i);    }
          return sum;    }

     int FullDevNN( ) const{    
       double sum = 0.0;
       for ( int i = 0; i < Ntigs( ) -1; i++ )    
	 sum += pow( Dev(i), 2.0 );     
       return round( sqrt(sum) );    
     }

     void Reverse( )
     {    mtig_.ReverseMe( );
          len_.ReverseMe( );
          gap_.ReverseMe( );    }

     void PlaceFirstTig( int m, int len )
     {    ForceAssertEq( mtig_.size( ), 0u );
          mtig_.push_back(m);
          len_.push_back(len);    }

     void AppendTig( int m, int len, int gap, int dev )
     {    mtig_.push_back(m);
          len_.push_back(len);
          gap_.push_back( gapb( gap, dev ) );    }

     void AppendTig( int m, int len )
     {    mtig_.push_back(m);
          len_.push_back(len);    }

     void AppendGap( int gap, int dev )
     {    gap_.push_back( gapb( gap, dev ) );    }

     void AppendSuper( const superb& s, int gap, int dev )
     {    for ( int i = 0; i < s.Ntigs( ); i++ )
          {    int g, d;
               if ( i == 0 )
               {    g = gap;
                    d = dev;    }
               else
               {    g = s.Gap(i-1);
                    d = s.Dev(i-1);    }
               AppendTig( s.Tig(i), s.Len(i), g, d );    }    }

     void InsertTig( int pos_of_new_tig,
                     int m, int len, 
                     int gap_before, int dev_before, 
                     int gap_after, int dev_after );

     // ReplaceTigBySuper: stuff a super into the place where a single contig was.

     void ReplaceTigBySuper( int tig_pos, const superb& s );

     // ReplaceGapBySuper: overloaded version.  
     // You have to specify the flanking gaps and deviations.

     void ReplaceTigBySuper( int tig_pos, const superb& s, int gap1, int dev1,
			     int gap2, int dev2 );
     
     // ReplaceTigsBySuper: stuff a super into the place where contigs were.

     void ReplaceTigsBySuper( int first_tig_pos, int last_tig_pos, const superb& s );

     // ReplaceGapBySuper: stuff a super into a gap together with possible 
     // additional edge gaps.  You have to specify the flanking gaps and deviations.

     void ReplaceGapBySuper( int gap_pos, const superb& s, int gap1, int dev1,
          int gap2, int dev2 );

     // RemoveTigByPos: remove a contig from a super.

     void RemoveTigByPos( int tig_pos );

     void RemoveTigByPosAlt( int tig_pos, int left_add, int right_add );

     // RemoveTigsByPos: remove a set of consecutive contigs from a super.

     void RemoveTigsByPos( int first_tig_pos, int last_tig_pos );
  
     // TigPosition: return the location in this super of the contig with this ID.
     // Return -1 if it's not found.
     int TigPosition( const int contig_ID ) const { return Position( mtig_, contig_ID ); }

     // SubSuperLength, SubSuperLengthDev: determine the estimated length and 
     // length s.d. for a given inclusive range of contigs.

     int SubSuperLength( int start, int stop ) const;
     int SubSuperLengthDev( int start, int stop ) const;

     // SubSuper: define a supercontig using only some of the contigs (defined
     // by "to_use") in a given supercontig.  We require that "to_use" is ordered.

     superb SubSuper( const vec<int>& to_use ) const;

     void Print( ostream& out, const String& super_name ) const;
     void Print( const Boolvector& s_rc, ostream& out, 
          const String& super_name ) const;
     // PrintRange: print scaffolding information between two contigs (in one line). 
     void PrintRange( ostream& out, const String& super_name, 
		      int tpos1, int tpos2) const;

     friend void WriteSuperbs( const String& fn, const vec<superb>& s );
     friend void ReadSuperbs( const String& fn, vec<superb>& s );

     // write the files  mergedcontigs.superb and mergedcontigs.summary
     friend void WriteSuperbsAndSummary( const String & dir,
					 const vec<superb> & s,
					 Bool index_only = False,
                                         const String& base = "mergedcontigs" );

     friend void WriteSummary( const String & fn, const vec<superb> & s, 
			       const vecBoolvector& s_rc = vecBoolvector() );

     // WriteSupercontigFiles: write all the supercontig assembly files.  These are
     // mergedcontigs.superb, mergedcontigs.summary, and mergedcontigs.answer.gz.
     // The last is only for backward compatibility.

     friend void WriteSupercontigFiles( const String& dir,
					const vec<superb>& s,
					Bool index_only = False );

     // Two supercontigs are declared abstractly equal if they have the same contig
     // sizes in the same or reverse order, with the same gaps and deviations.
     // We allow total variation of up to max_diff_percent, where the denominator 
     // is the total contig length.

     friend Bool AbstractlyEqual( const superb& s1, const superb& s2,
          float max_diff_percent );

     private:

     // Data definitions:

     vec<int> mtig_;
     vec<int> len_;
     vec<gapb> gap_;

};

class superb_by_tig {

     public:

     superb_by_tig( ) { }

     superb_by_tig( int super, int super_pos, int super_len )
          : super_(super), super_pos_(super_pos), super_len_(super_len) { }

     // Super: return numerical identifier of supercontig.

     int Super( ) const { return super_; }

     // SuperPos: return index of contig in the supercontig.

     int SuperPos( ) const { return super_pos_; }

     // SuperLen: return the number of contigs in the supercontig.

     int SuperLen( ) const { return super_len_; }

     private:

     int super_;
     int super_pos_;
     int super_len_;

};

// Test to see if mergedcontigs.superb_index seems to be up to date.

void TestSuperbIndex( const String& dir );

// Fetch superb index entries from disk.

void FetchSuperbIndexEntries( const String& dir, const vec<int>& tigs,
     vec<superb_by_tig>& indices );

#endif
