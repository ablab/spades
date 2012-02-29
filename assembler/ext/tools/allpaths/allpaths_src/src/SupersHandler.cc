///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "Superb.h"
#include "SupersHandler.h"
#include "math/Functions.h"

/*
 * shandler
 * Constructor
 */
shandler::shandler( int n_contigs ) :
  n_contigs_ ( n_contigs )
{ }

/*
 * shandler
 * Constructor
 */
shandler::shandler( int n_contigs, const String &supers_file ) :
  n_contigs_ ( n_contigs )
{
  this->LoadFromFile( supers_file );
}

/*
 * shandler
 * Constructor
 *
 * This constructor will make a copy of supers.
 */
shandler::shandler( const vec<superb> &supers, int *min_gap ) :
  n_contigs_ ( -1 ),
  supers_ ( supers )
{
  this->Setup( min_gap );
}

/*
 * shandler
 * LoadFromFile
 *
 * If n_contigs = -1 then the number of contigs is deduced from the
 * superb file itself, as the largest contig_id plus one.
 */
void shandler::LoadFromFile( const String &supers_file, int *min_gap )
{
  ReadSuperbs( supers_file, supers_ );
  this->Setup( min_gap );
}

/*
 * shandler
 * Setup
 *
 * If n_contigs_ = -1 then the number of contigs is deduced from the
 * superb file itself, as the largest contig_id plus one. If min_gap
 * is given, reset gaps to >= *min_gap (stdevs are not changed).
 */
void shandler::Setup( int *min_gap )
{
  if ( n_contigs_ < 0 ) {
    for (int ii=0; ii<(int)supers_.size( ); ii++)
      for (int jj=0; jj<supers_[ii].Ntigs( ); jj++)
	n_contigs_ = Max( n_contigs_, supers_[ii].Tig( jj ) );
    n_contigs_ += 1;
  }
  
  if ( min_gap ) {
    for (int ii=0; ii<(int)supers_.size( ); ii++) {
      for (int jj=0; jj<supers_[ii].Ntigs( )-1; jj++)
	if ( supers_[ii].Gap( jj ) < *min_gap )
	  supers_[ii].SetGap( jj, *min_gap );
    }
  }

  true_begins_.clear( );
  true_ends_.clear( );
  to_super_.clear( );
  start_on_super_.clear( );
  pos_on_super_.clear( );

  true_begins_.resize( supers_.size( ), 0 );
  true_ends_.resize( supers_.size( ), 0 );
  to_super_.resize( n_contigs_, -1 );
  start_on_super_.resize( n_contigs_, -1 );
  pos_on_super_.resize( n_contigs_, -1 );

  for (int ii=0; ii<(int)supers_.size( ); ii++) {
    true_begins_[ii] = supers_[ii].TrueBegin( );
    true_ends_[ii] = supers_[ii].TrueEnd( );
    int pos = 0;
    for (int jj=0; jj<(int)supers_[ii].Ntigs( ); jj++) {
      ForceAssertEq( to_super_[supers_[ii].Tig( jj )], -1 );
      to_super_[supers_[ii].Tig( jj )] = ii;
      start_on_super_[supers_[ii].Tig( jj )] = pos;
      pos_on_super_[supers_[ii].Tig( jj )] = jj;
      pos += supers_[ii].Len( jj );
      if ( jj < (int)supers_[ii].Ntigs( ) - 1 )
	pos += supers_[ii].Gap( jj );
    }
  }

}

