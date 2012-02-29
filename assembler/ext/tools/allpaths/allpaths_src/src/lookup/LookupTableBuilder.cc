/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// LookupTableBuilder. Builds a kmer lookup table for all the k-mers in a
/// given set of basevectors.

#include "MainTools.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "lookup/LookupTable.h"
#include "lookup/LookupTableBuilder.h"


/** \brief Creates a new lookup table with specified sequence, 
 *  Kmer length and chunk length/overlap.
 *
 *  Created lookup table will have a separate contig for each separate
 *  basevector in \c source. Contigs will be assigned default names 
 *  "contig_<n>", where <n> is the ordinal number of the contig.
 *
 *  @param source collection of contigs (each basevector = contig )
 *  @param table_name name of the disk file to write the lookup table to
 *  @param IGNORE_SHORT if \c true, the function will tolerate and append
 *   contigs that are too short ( < K ), otherwise
 *  @param K Kmer size for this lookup table. Default = 12
 *   the function will abort if such contig is found. Default = \c False
 *  @param CHUNK_SIZE size of the chunks to be created; default = 100M
 *  @param CHUNK_OVERLAP overlap between chunks (needed to make sure 
 *  alignments across a chunk boundary are not lost); default = 1M
 */
void LookupTableBuilder(const vecbasevector& source, 
			const String table_name, // file name to write to
			const Bool IGNORE_SHORT, 
			const unsigned int K,
			const unsigned int CHUNK_SIZE, 
			const unsigned int CHUNK_OVERLAP) {

  vecString cnames;

  // if no contig names were passed, create simple default ones...
  for ( size_t j = 0; j < source.size( ); j++ ) {
    cnames.push_back( "contig_" + ToString(j));
  }

  // ...and pass everything further to the function that actually does the job
  LookupTableBuilder(source, cnames, table_name, IGNORE_SHORT,
		     K, CHUNK_SIZE, CHUNK_OVERLAP);

}


/** \brief Creates a new lookup table with specified sequence, 
 *  Kmer length and chunk length/overlap.
 *
 *  Created lookup table will have a separate contig for each separate
 *  basevector in \c source. Contigs will be assigned names passed in the 
 *  \c contig_names vector.
 *
 *  @param source collection of contigs (each basevector = contig )
 *  @param contig_names names of the contigs (will be stored in the lookup 
 *   table)
 *  @param table_name name of the disk file to write the lookup table to
 *  @param IGNORE_SHORT if \c true, the function will tolerate and append
 *   contigs that are too short ( < K ), otherwise
 *  @param K Kmer size for this lookup table. Default = 12
 *   the function will abort if such contig is found. Default = \c False
 *  @param CHUNK_SIZE size of the chunks to be created; default = 100M
 *  @param CHUNK_OVERLAP overlap between chunks (needed to make sure 
 *  alignments across a chunk boundary are not lost); default = 1M
 */
void LookupTableBuilder(const vecbasevector& source, // contigs to build the ref from
			const vecString &contig_names,
			const String table_name, // file name to write to
			const Bool IGNORE_SHORT, 
			const unsigned int K,
			const unsigned int CHUNK_SIZE, 
			const unsigned int CHUNK_OVERLAP) 
{
  
  // Impose requirement on K.
  if ( K > 15 ) {
    cout << "The maximum allowed K value is 15.  Abort.\n";
    exit(1);
  }
  
  // make sure that numbers of contigs and contig names are the same
  if ( source.size() != contig_names.size() ) {
    cout << "LookupTableBuilder: counts of passed contigs and ";
    cout << "of contig names are not the same." << endl;
    exit(1);
  }

  // Create Lookup file
  Remove( table_name );
  int fd = Open( table_name, O_WRONLY | O_CREAT );

  // Count number of bases.
  const longlong max_bases(4000000000u);
  longlong base_count = 0;
  for ( size_t i = 0; i < source.size( ); i++ )
    base_count += source[i].size( );    
  if ( base_count >= max_bases ) {
    cout << "There are more than four billion bases in "
	 << "your fasta files.  Abort.\n";
    exit(1);
  }

  // Set up lookup table.
  lookup_table look(fd);
  look.SetK(K);
  look.SetChunkParams( CHUNK_SIZE, CHUNK_OVERLAP );

  // Examine source contigs
  vec<char> bases(base_count);
  base_count = 0;

  for ( size_t j = 0; j < source.size( ); j++ ) {
    unsigned int n = source[j].size( );
    look.AddContigName( contig_names[j], "", j );
    look.AddContigStart(base_count);
    look.AddContigSize(n);
    if ( ( n < K && !IGNORE_SHORT ) || n > 2000000000u ) {    
      cout << "Size of contig \"" 
	   << look.LastContigName( ) << "\" is " << n;
      if ( n < K ) 
	cout << ", which is too small (min value = K).\n";
      else 
	cout << ", which is too large "
	     << "(max value = 2,000,000,000).\n";
      exit(1);    
    }
    for ( unsigned int k = 0; k < n; k++ )
      bases[ base_count++ ] = as_base( source[j][k] );    
  }

  // Build the table
  BuildTableFromContigs(look, bases, K, CHUNK_SIZE, CHUNK_OVERLAP);
}


void BuildTableFromContigs(lookup_table& look, 
			   const vec<char>& bases, 
			   const unsigned int K,
			   const unsigned int CHUNK_SIZE,
			   const unsigned int CHUNK_OVERLAP ) 
{
  // index_loc is a list of (kmer number, offset in reference) pairs
  vec< pair<unsigned int, unsigned int> > index_loc;
  index_loc.reserve(CHUNK_SIZE);
  look.WriteHeader( );

  unsigned int base_count_in_chunk = 0;
  for ( unsigned int i = 0; i < look.NContigs( ); i++ ) {
    if ( look.ContigSize(i) < K )
      continue;
    for ( unsigned int j = 0; j < look.ContigSize(i) - K + 1; j++ ) {
      unsigned int pos = look.ContigStart(i) + j;
      if (base_count_in_chunk > 0 && base_count_in_chunk % CHUNK_SIZE == 0) {
	// Dump current chunk of data to disk and reset
	look.DumpChunk( index_loc, bases );
	if ( j < CHUNK_OVERLAP )
	  j = 0;
	else
	  j = j - CHUNK_OVERLAP;
	pos = look.ContigStart(i) + j;
	base_count_in_chunk = 0;
      }
      ++base_count_in_chunk;
      
      Bool ambiguous = False;
      for ( unsigned int l = 0; l < K; l++ ) {
	if ( AmbiguousBase( bases[pos + l] ) ) {
	  ambiguous = True;
	  break;
	}
      }
      if (ambiguous)
	continue;
      
      // Compute index of kmer, put entry in table.
      unsigned int index = Index( bases, pos, K );
      index_loc.push_back( make_pair( index, pos ) );
    }
  }
  // Dump final chunk of data to disk and reset
  look.DumpChunk( index_loc, bases );
  look.WriteHeader( );
}

