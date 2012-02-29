///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Intvector.h"
#include "util/SearchFastb2Core.h"
#include "paths/MapUnibasesToContigs.h"

/**
 * MapUnibasesToContigs
 */
void MapUnibasesToContigs( const int K,
			   const String &unibases_file,
			   const String &contigs_file,
			   const String &OUTBASE,
			   ostream *log )
{
  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;
  
  // File names.
  String kmers_file = OUTBASE + ".kmers.fastb";
  String map_file = OUTBASE + ".u2c";

  // Count kmers in contigs.
  out << Date( ) << ": counting kmers in contigs... " << flush;
  size_t nkmers = 0;
  vecbvec contigs( contigs_file );
  for (size_t ii=0; ii<contigs.size( ); ii++)
    if ( contigs[ii].isize( ) >= K )
      nkmers += contigs[ii].size( ) - K + 1;
  out << ToStringAddCommas( nkmers ) << " found" << endl;
  
  // Generate temp fastb files.
  out << Date( ) << ": generating and saving temp kmers file" << endl;
  vecbvec kmers;
  vec<size_t> to_contig;
  kmers.reserve( nkmers );
  to_contig.reserve( nkmers );
  for (size_t ii=0; ii<contigs.size( ); ii++) {
    for (int jj=0; jj<contigs[ii].isize( ) - K + 1; jj++) {
      kmers.push_back( bvec( contigs[ii], jj, K ) );
      to_contig.push_back( ii );
    }
  }
  kmers.WriteAll( kmers_file );
  
  // Align and sort (remove rc's).
  vec< triple<int64_t,int64_t,int> > aligns;
  out << Date( ) << ": running SearchFastb2... " << flush;
  SearchFastb2( kmers_file, unibases_file, K, &aligns, 0, -1, 0.9, False );
  out << ToStringAddCommas( aligns.size( ) ) << " aligns found" << endl;
  {
    out << Date( ) << ": removing rc aligns, and sorting" << endl;
    vec< triple<int64_t,int64_t,int> > fwaligns;
    fwaligns.reserve( aligns.size( ) / 2 );
    for (size_t ii=0; ii<aligns.size( ); ii++)
      if ( aligns[ii].third > -1 )
	fwaligns.push_back( aligns[ii] );
    sort( fwaligns.begin( ), fwaligns.end( ) );
    swap( fwaligns, aligns );
  }
  
  // Remove kmers file.
  out << Date( ) << ": removing temp kmers file" << endl;
  Remove( kmers_file );
  
  // Generate maps.
  out << Date( ) << ": generating contigs to unibases map" << endl;
  vec< pair<size_t,size_t> > c2u;
  for (size_t ii=0; ii<aligns.size( ); ii++) {
    size_t kmer_id = aligns[ii].first;
    size_t unibase_id = aligns[ii].second;
    size_t contig_id = to_contig[kmer_id];
    if ( c2u.size( ) < 1 ||
	 contig_id != c2u.back( ).first ||
	 unibase_id != c2u.back( ).second )
      c2u.push_back( make_pair( contig_id, unibase_id ) );
  }
  sort( c2u.begin( ), c2u.end( ) );
  c2u.erase( unique( c2u.begin( ), c2u.end( ) ), c2u.end( ) );

  out << Date( ) << ": generating unibases to contigs map" << endl;
  size_t n_unibases = MastervecFileObjectCount( unibases_file );
  UInt64VecVec u2c( n_unibases );
  for (size_t ii=0; ii<c2u.size( ); ii++)
    u2c[ c2u[ii].second ].push_back( c2u[ii].first );
  
  // Save map.
  out << Date( ) << ": saving unibases to contigs map" << endl;
  u2c.WriteAll( map_file );
  
  // Done.
  out << Date( ) << ": done" << endl;
  
}
