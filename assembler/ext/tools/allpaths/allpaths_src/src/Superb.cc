///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "AnnotatedContig.h"
#include "Basevector.h"
#include "Charvector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "Superb.h"
#include "math/Functions.h"
#include "system/file/FileReader.h"

#include <math.h>

void WriteSuperbs( const String& fn, const vec<superb>& s )
{    ForceAssert( !IsRegularFile( fn + ".gz" ) );
     int fd = Open( fn, O_WRONLY | O_CREAT );
     int n = s.size( );
     WriteBytes( fd, &n, sizeof(int) );
     for ( int i = 0; i < n; i++ )
     {    int ntigs = s[i].Ntigs( );
          if ( ntigs == 0 )
          {    cout << "In superb::Write: attempt to write super " << i
                    << ", which has zero contigs.  Abort.\n";
               TracebackThisProcess( );    }
          WriteBytes( fd, &ntigs, sizeof(int) );
          WriteBytes( fd, &s[i].mtig_[0], ntigs * sizeof(int) );
          WriteBytes( fd, &s[i].len_[0], ntigs * sizeof(int) );
          if ( ntigs > 1 )
               WriteBytes( fd, &s[i].gap_[0], (ntigs-1) * sizeof(superb::gapb) );   }
     close(fd);    }

void ReadSuperbs( const String& fn, vec<superb>& s )
{    FileReader fr(fn.c_str());
     int n;
     fr.read( &n, sizeof(int) );
     s.resize(n);
     for ( int i = 0; i < (int) s.size( ); i++ )
     {    int ntigs;
          fr.read( &ntigs, sizeof(int) );
          s[i].mtig_.resize(ntigs);
          s[i].len_.resize(ntigs);
          s[i].gap_.resize( ntigs - 1 );
          fr.read( &s[i].mtig_[0], ntigs * sizeof(int) );
          fr.read( &s[i].len_[0], ntigs * sizeof(int) );
          if ( ntigs > 1 ) 
               fr.read( &s[i].gap_[0], (ntigs-1) * sizeof(superb::gapb) ); }  }

void superb::InsertTig( int pos_of_new_tig,
                        int m, int len, 
                        int gap_before, int dev_before, 
                        int gap_after, int dev_after )
{
    ForceAssertLe( pos_of_new_tig, (int) mtig_.size() );
    
    // If the new contig won't be the last contig, insert the gap after it.
    if ( pos_of_new_tig < (int) mtig_.size() )
        gap_.insert( gap_.begin() + pos_of_new_tig, gapb( gap_after, dev_after ) );
    
    // If the new contig will be the last contig, resize the gap_ vector.
    if ( pos_of_new_tig == (int) mtig_.size() )
        gap_.resize( gap_.size() + 1 );
    
    // Set the gap before the new contig.
    if ( pos_of_new_tig > 0 )
        gap_[ pos_of_new_tig - 1 ] = gapb( gap_before, dev_before );
    
    mtig_.insert( mtig_.begin() + pos_of_new_tig, m );
    len_.insert( len_.begin() + pos_of_new_tig, len );
}

void superb::ReplaceTigBySuper( int tig_pos, const superb& s )
{    ForceAssertGe( tig_pos, 0 );
     ForceAssertLt( tig_pos, Ntigs( ) );
     if ( s.Ntigs( ) == 0 )
     {    RemoveTigByPos(tig_pos);
          return;    }
     int push = s.Ntigs( ) - 1;
     mtig_.resize( mtig_.size( ) + push );
     len_.resize( len_.size( ) + push );
     gap_.resize( gap_.size( ) + push );
     for ( int j = Ntigs( ) - 1; j >= tig_pos + 1 + push; j-- )
     {    SetTig( j, Tig( j - push ) );
          SetLen( j, Len( j - push ) );    }
     for ( int j = Ntigs( ) - 2; j >= tig_pos + push; j-- )
     {    SetGap( j, Gap( j - push ) );
          SetDev( j, Dev( j - push ) );    }
     for ( int j = 0; j < s.Ntigs( ); j++ )
     {    SetTig( tig_pos + j, s.Tig(j) );
          SetLen( tig_pos + j, s.Len(j) );
          if ( j < s.Ntigs( ) - 1 )
          {    SetGap( tig_pos + j, s.Gap(j) );
               SetDev( tig_pos + j, s.Dev(j) );    }    }    }

void superb::ReplaceTigBySuper( int tig_pos, const superb& s, 
				int gap1, int dev1,
				int gap2, int dev2){    
  ForceAssertGe( tig_pos, 0 );
  ForceAssertLt( tig_pos, Ntigs( ) );
  if ( s.Ntigs( ) == 0 ){    
     ForceAssertGe( tig_pos, 0 );
     ForceAssertLt( tig_pos, this->Ntigs() );
     if ( tig_pos != 0 &&  tig_pos != (int) gap_.size() ){
       gap_[ tig_pos - 1 ].SetGap( gap_[ tig_pos - 1 ].Gap() + 
				   len_[ tig_pos ] + 
				   gap_[ tig_pos ].Gap() );
       
       gap_[ tig_pos - 1 ].SetDev( gap_[ tig_pos - 1 ].Dev() +  
				   gap_[ tig_pos ].Dev() );
     }

     mtig_.erase( mtig_.begin() + tig_pos );
     len_.erase( len_.begin() + tig_pos );
     
     if ( tig_pos != (int) gap_.size() ) 
       gap_.erase( gap_.begin() + tig_pos );
     return;    
  }
  int NtigsOrig = Ntigs();
  int push      = s.Ntigs( ) - 1;
  mtig_.resize( mtig_.size( ) + push );
  len_.resize( len_.size( ) + push );
  gap_.resize( gap_.size( ) + push );
  for ( int j = Ntigs( ) - 1; j >= tig_pos + 1 + push; j-- ){    
    SetTig( j, Tig( j - push ) );
    SetLen( j, Len( j - push ) );    
  }
  for ( int j = Ntigs( ) - 2; j >= tig_pos + push; j-- ){    
    SetGap( j, Gap( j - push ) );
    SetDev( j, Dev( j - push ) );    
  }
  for ( int j = 0; j < s.Ntigs( ); j++ ){
    SetTig( tig_pos + j, s.Tig(j) );
    SetLen( tig_pos + j, s.Len(j) );
    if ( j == 0 && gap1 != 0 && tig_pos != 0 ){
      int gap = Gap( tig_pos - 1 );
      int dev = Dev( tig_pos - 1 );
      gap += gap1;
      // assuming independence
      dev  = sqrt( pow(dev, 2.0) + pow(dev1, 2.0) );
      
      SetGap( tig_pos - 1, gap );
      SetDev( tig_pos - 1, dev ); 
    }
    if ( j < s.Ntigs( ) - 1 ){    
      SetGap( tig_pos + j, s.Gap(j) );
      SetDev( tig_pos + j, s.Dev(j) );   
    }
    if ( j == s.Ntigs() -1 && tig_pos != NtigsOrig -1 ){
      int gap = Gap( tig_pos + j );
      int dev = Dev( tig_pos + j );
      if ( gap2 != 0 ){
	gap += gap2;
	// assuming independence
	dev = sqrt( pow(dev, 2.0) + pow(dev2, 2.0) );
      }
      SetGap( tig_pos + j, gap );
      SetDev( tig_pos + j, dev ); 
    }
  }    
}

void superb::ReplaceTigsBySuper( int first_tig_pos, int last_tig_pos, 
				 const superb& s ){    
  ForceAssertLe( first_tig_pos, last_tig_pos );
  ForceAssertGe( first_tig_pos, 0 );
  for ( int tig_pos = last_tig_pos; tig_pos > first_tig_pos; --tig_pos ){
    mtig_.erase( mtig_.begin() + tig_pos );
    len_.erase( len_.begin() + tig_pos );
    if ( tig_pos -1 != (int) gap_.size() ) 
      gap_.erase( gap_.begin() + tig_pos -1 );
  }
  
  this->ReplaceTigBySuper( first_tig_pos, s );
}


void superb::ReplaceGapBySuper( int gap_pos, const superb& s, int gap1, int dev1,
     int gap2, int dev2 )
{    int push = s.Ntigs( );
     mtig_.resize( mtig_.size( ) + push );
     len_.resize( len_.size( ) + push );
     gap_.resize( gap_.size( ) + push );
     for ( int j = Ntigs( ) - 1; j >= gap_pos + 1 + push; j-- )
     {    SetTig( j, Tig( j - push ) );
          SetLen( j, Len( j - push ) );    }
     for ( int j = Ntigs( ) - 2; j >= gap_pos + push + 1; j-- )
     {    SetGap( j, Gap( j - push ) );
          SetDev( j, Dev( j - push ) );    }
     for ( int j = 0; j < s.Ntigs( ); j++ )
     {    SetTig( gap_pos + j + 1, s.Tig(j) );
          SetLen( gap_pos + j + 1, s.Len(j) );
          if ( j < s.Ntigs( ) - 1 )
          {    SetGap( gap_pos + j + 1, s.Gap(j) );
               SetDev( gap_pos + j + 1, s.Dev(j) );    }    }
     SetGap( gap_pos, gap1 );
     SetDev( gap_pos, dev1 );
     SetGap( gap_pos + push, gap2 );
     SetDev( gap_pos + push, dev2 );    }

void superb::RemoveTigByPos( int tig_pos )
{
    ForceAssertGe( tig_pos, 0 );
    ForceAssertLt( tig_pos, this->Ntigs() );
    
    if ( tig_pos != 0 &&
         tig_pos != (int) gap_.size() )
    {
        gap_[ tig_pos - 1 ].SetGap( gap_[ tig_pos - 1 ].Gap() + 
                                    len_[ tig_pos ] + 
                                    gap_[ tig_pos ].Gap() );
	// assuming independence of gap size distributions
	gap_[ tig_pos - 1 ].SetDev( sqrt( pow(gap_[ tig_pos - 1 ].Dev(), 2.0) + 
					  pow(gap_[ tig_pos ].Dev(), 2.0)       )   );        
    }

    mtig_.erase( mtig_.begin() + tig_pos );
    len_.erase( len_.begin() + tig_pos );

    if ( tig_pos != (int) gap_.size() ) 
        gap_.erase( gap_.begin() + tig_pos );
}

void superb::RemoveTigByPosAlt( int tig_pos, int left_add, int right_add )
{    ForceAssertGe( tig_pos, 0 );
     ForceAssertLt( tig_pos, Ntigs( ) );
     if ( tig_pos-1 >= 0 && tig_pos != Ngaps( ) )
          SetGap( tig_pos-1, Gap(tig_pos) + right_add );
     if ( tig_pos-2 >= 0 ) SetGap( tig_pos-2, Gap(tig_pos-2) + left_add );
     mtig_.erase( mtig_.begin( ) + tig_pos );
     len_.erase( len_.begin( ) + tig_pos );
     if ( tig_pos != Ngaps( ) ) gap_.erase( gap_.begin( ) + tig_pos );    }

void superb::RemoveTigsByPos( int first_tig_pos, int last_tig_pos )
{
    ForceAssertLe( first_tig_pos, last_tig_pos );

    for ( int tig_pos = last_tig_pos; tig_pos >= first_tig_pos; --tig_pos )
        this->RemoveTigByPos( tig_pos );
}

void superb::Print( ostream& out, const String& super_name ) const
{    int n = Ntigs( );
     int total_super = 0;
     for ( int j = 0; j < n; j++ )
     {    total_super += Len(j);
          if ( j < n - 1 ) total_super += Gap(j);    }
     out << super_name << " (l = " << total_super << "):\n";
     for ( int j = 0; j < n; j++ )
     {    int m = Tig(j);
          out << m << " (l = " << Len(j) << ")\n";
          if ( j < n - 1 ) 
               out << " -- (" << Gap(j) << " +/- " << Dev(j) << ") --> ";    }
     out << "\n";    }

void superb::Print( const Boolvector& s_rc, ostream& out, 
     const String& super_name ) const
{    int n = Ntigs( );
     int total_super = 0;
     for ( int j = 0; j < n; j++ )
     {    total_super += Len(j);
          if ( j < n - 1 ) total_super += Gap(j);    }
     out << super_name << " (l = " << total_super << "):\n";
     for ( int j = 0; j < n; j++ )
     {    int m = Tig(j);
          out << m << "[" << ( s_rc[j] ? "-" : "+" ) << "]"
               << " (l = " << Len(j) << ")\n";
          if ( j < n - 1 ) 
               out << " -- (" << Gap(j) << " +/- " << Dev(j) << ") --> ";    }
     out << "\n";    }

void superb::PrintRange( ostream& out, const String& super_name,
			 int tpos1, int tpos2 ) const {    
  int n = Ntigs( );
  if ( tpos1 < 0 )    tpos1 = 0;
  if ( tpos2 > n -1 ) tpos2 = n -1;
 
  int total_super = 0;
  for ( int j = tpos1; j <= tpos2; j++ ){    
    total_super += Len(j);
    if ( j < n - 1 ) total_super += Gap(j);    
  }
  out << super_name << " between positions [" << tpos1 << "," << tpos2 << "] combined length = " << total_super << "\n";
  for ( int j = tpos1; j <= tpos2; j++ ){    
    int m = Tig(j);
    out << m << " (l = " << Len(j) << ") ";
    if ( j < tpos2 ) 
      out << " -- (" << Gap(j) << " +/- " << Dev(j) << ") --> ";    
  }
  out << "\n";    
}

superb superb::SubSuper( const vec<int>& to_use ) const
{    ForceAssert( to_use.UniqueOrdered( ) );
     ForceAssert( to_use.nonempty( ) );
     superb s;
     s.SetNtigs( to_use.size( ) );
     for ( size_t i = 0; i < to_use.size(); i++ )
     {    s.SetTig( i, Tig( to_use[i] ) );
          s.SetLen( i, Len( to_use[i] ) );    }
     for ( size_t i = 1; i < to_use.size(); i++ )
     {    int gap = 0, dev = 0;
          for ( int j = to_use[i-1] + 1; j <= to_use[i] - 1; j++ )
               gap += Len(j);
          for ( int j = to_use[i-1]; j <= to_use[i] - 1; j++ )
          {    gap += Gap(j);
               dev += Dev(j);    }
          s.SetGap( i - 1, gap );
          s.SetDev( i - 1, dev );    }
     return s;    }

void TestSuperbIndex( const String& dir )
{    FileReader fr( (dir + "/mergedcontigs.superb_index").c_str() );
     int entries;
     fr.read( &entries, sizeof(int) );
     int file_size = FileSize( dir + "/mergedcontigs.superb_index" );
     if ( file_size != (int) sizeof(int) + entries * (int) sizeof(superb_by_tig) )
     {    FatalErr( "The file mergedcontigs.superb_index in " << dir
               << " seems to be out of date." );    }    }

void FetchSuperbIndexEntries( const String& dir, const vec<int>& tigs,
     vec<superb_by_tig>& indices )
{    TestSuperbIndex(dir);
     BinaryReadSubset( dir + "/mergedcontigs.superb_index", tigs, indices );    }

void WriteSummary( const String & fn, const vec<superb> & s, const vecBoolvector& s_rc )
{    Ofstream( summary, fn );
     summary <<
     "This file contains a list of all the supercontigs in the final assembly.\n"
     "For each, we give the numbers and lengths of the constituent contigs, "
     "as well \n"
     "as the approximate gaps between them (if a positive number is shown in \n"
     "parentheses) or the approximate overlap between them (if a negative number \n"
     "is shown).\n\n";
     for ( unsigned int i = 0; i < s.size( ); i++ )
     {    if ( s_rc.size( ) == 0 ) s[i].Print( summary, "scaffold " + ToString(i) );
          else s[i].Print( s_rc[i], summary, "scaffold " + ToString(i) );    }    }

void WriteSuperbsAndSummary( const String & dir, const vec<superb> & s,
     Bool index_only, const String& base )
{
     // Write mergedcontigs.superb_index.

     int ntigs = 0;
     for ( size_t i = 0; i < s.size(); i++ )
          for ( int j = 0; j < s[i].Ntigs( ); j++ )
               ntigs = Max( ntigs, 1 + s[i].Tig(j) );
     vec<superb_by_tig> by_tig(ntigs);
     for ( size_t i = 0; i < s.size(); i++ )
          for ( int j = 0; j < s[i].Ntigs( ); j++ )
               by_tig[ s[i].Tig(j) ] = superb_by_tig( i, j, s[i].Ntigs( ) );
     BinaryWrite( dir + "/" + base + ".superb_index", by_tig );
     if (index_only) return;

     // Write mergedcontigs.superb and mergedcontigs.summary.

     WriteSuperbs( dir + "/" + base + ".superb", s );
     WriteSummary( dir + "/" + base + ".summary", s );    }

void WriteSupercontigFiles( const String& dir,
			    const vec<superb>& s, 
			    Bool index_only )
{    
     WriteSuperbsAndSummary( dir, s, index_only );

     // Backward-compatible code:

     PipeOstream( answer, dir + "/mergedcontigs.answer" );
     vec<annotated_supercontig> annotated_supers;
     for ( unsigned int i = 0; i < s.size( ); i++ )
     {    const superb& b = s[i];
          int n = b.Ntigs( );
          vec<annotated_contig> current_contigs;
          vec<int> current_gaps( n - 1 ), current_gaps_sd( n - 1, 0 );
          for ( int i = 0; i < n - 1; i++ )
               current_gaps[i] = b.Gap(i);
          for ( int j = 0; j < n; j++ )
          {    int m = b.Tig(j);
               vec<arachne_contig_placement> placements;
               current_contigs.push_back( annotated_contig( m,
                    basevector(0), qualvector(0),
                    False, False, vec<semiannotation>(0),
                    vec<annotation>(0), placements ) );    }
          annotated_supercontig a( current_contigs, current_gaps, current_gaps_sd );
          annotated_supers.push_back(a);    }
     answer << annotated_supers.size( ) << "\n";
     for ( unsigned int i = 0; i < annotated_supers.size( ); i++ )
          answer << annotated_supers[i] << "\n";    }

int superb::SubSuperLength( int start, int stop ) const
{    int sum = 0;
     for ( int i = start; i <= stop; i++ )
     {    sum += Len(i);
          if ( i < stop ) sum += Gap(i);    }
     return sum;    }

int superb::SubSuperLengthDev( int start, int stop ) const
{    if ( start == stop ) return 0;
     double var = 0.0;
     for ( int i = start; i < stop; i++ )
          var += (double) Dev(i) * (double) Dev(i);
     return int(rint(sqrt(var)));    }

Bool AbstractlyEqual( const superb& s1, const superb& s2, float max_diff_percent )
{    if ( s1.Ntigs( ) != s2.Ntigs( ) ) return False;
     int n = s1.Ntigs( );
     int delta = 0;
     int max_allowed_diff = int( floor( max_diff_percent/100.0 
          * float( Min( s1.ReducedLength( ), s2.ReducedLength( ) ) ) ) );
     for ( int i = 0; i < n; i++ )
     {    delta += Abs( s1.Len(i) - s2.Len(i) );
          if ( i < n - 1 ) delta += Abs( s1.Gap(i) - s2.Gap(i) );
          if ( i < n - 1 ) delta += Abs( s1.Dev(i) - s2.Dev(i) );    }
     if ( delta <= max_allowed_diff ) return True;
     delta = 0;
     for ( int i = 0; i < n; i++ )
     {    delta += Abs( s1.Len(i) - s2.Len(n-i-1) );
          if ( i < n - 1 ) delta += Abs( s1.Gap(i) - s2.Gap(n-i-2) );
          if ( i < n - 1 ) delta += Abs( s1.Dev(i) - s2.Dev(n-i-2) );    }
     if ( delta <= max_allowed_diff ) return True;
     return False;    }
