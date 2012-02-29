/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Bitvector.h"
#include "lookup/LookupTable.h"

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
void LookupTableBuilder(const vecbasevector& source, // contigs to put into the ref
			const String table_name,
			const Bool IGNORE_SHORT = False,
			const unsigned int K = 12, 
			const unsigned int CHUNK_SIZE = 100 * 1000 * 1000,
			const unsigned int CHUNK_OVERLAP = 1000 * 1000 );


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
void LookupTableBuilder(const vecbasevector& source,
			const vecString & contig_names,
			const String table_name,
			const Bool IGNORE_SHORT = False,
			const unsigned int K = 12, 
			const unsigned int CHUNK_SIZE = 100 * 1000 * 1000,
			const unsigned int CHUNK_OVERLAP = 1000 * 1000 );


/** BuildTableFromContigs. Completes construction of the lookup table.
    Called by LookupTableBuilder and MakeLookupTable.cc
*/
void BuildTableFromContigs(lookup_table& look, const vec<char>& bases, 
			   const unsigned int K,
			   const unsigned int CHUNK_SIZE,
			   const unsigned int CHUNK_OVERLAP );
