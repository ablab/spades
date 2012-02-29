///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"

#include "lookup/LookupTable.h"

/** \brief Sets the size of the lookup table's K-mer to \c k.
 *
 *  This method sets the K-mer size and automatically performs
 *  all required synchronization (internally maintained constants
 *  \c four_to_K, \c Kmask_ are updated; \c lookup and \c freq tables are
 *  resized - but \i not initialized, there is no data yet). Use this method
 *  only when creating new lookup tables and only prior to adding
 *  sequence data.
 */
void lookup_table::SetK( unsigned int K )
{    K_ = K;
     four_to_K_ = 1;
     Kmask_ = KmerBitmask(K);
     for ( unsigned int i = 0; i < K; i++ )
          four_to_K_ *= 4;
     lookup_.resize( four_to_K_ );
     freq_.resize( four_to_K_, 0 );
     control_[0] = K;
}

/** \brief Adds a new conting name to the lookup table.
 *
 *  This is a dangerous "direct access" method that leaves lookup table
 *  in an inconsistent state. Use it \i only along with other setter methods
 *  when building a table.
 *
 *  This method adds the specified name for the new contig to the list of
 *  contig names maintained by the lookup table (name must be shorter than
 *  252 symbols, otherwise the name is truncated at 251 characters with the
 *  last 3 symbols, [248...250] replaced with '...'). Simultaneously,
 *  an alternative name of the form "<index_in_file>:<file_name>" is added
 *  for the same contig (must be shorter than 256 characters). The counter
 *  of the number of contigs in this lookup table is incremented. All other
 *  data structures (tables of contig start positions and sizes,
 *  list of all bases for this lookup table) stay \i unchanged.
 *
 *  @see AddContigStart()
 *  @see AddContigSize()
 *  @see SetChunkParams()
 *  @see DumpChunk()
 *  @see LookupTableBuilder::BuildTableFromContigs()
 */
void lookup_table::AddContigName( String contig_name, const String&
     file_name, int index_in_file )
{    if ( contig_name.size( ) >= 252 )
          contig_name = contig_name.substr(0, 247) + " ...";
     contig_name_.push_back(contig_name);
     String contig_name_alt = ToString(index_in_file) + ":" + file_name;
     ForceAssertLt( contig_name_alt.size( ), 256u );
     contig_name_alt_.push_back(contig_name_alt);
     ++control_[2];
}

/** \brief Reads in the lookup table header from disk and initializes
 *  associated data structures (K, number, names and sizes of contigs,
 *  number of chunks, etc).
 *
 *  This method is intended for use in constructor only. The data stored in
 *  the file header are read in, and additional data maintained by the lookup
 *  table are initialized (such as chunk start offsets in the file,
 *  absolute positions of the contig starts in the full concatenated reference
 *  sequence, etc). All accessor
 *  methods defined in this class should return correct values after this
 *  method is executed.
 */
void lookup_table::ReadHeader( )
{
     ForceAssert(-1!=fd_);
     off_t bytes_in_header_ = 0;
     bytes_in_header_ += read( fd_, &control_[0], 1024 );
     if ( bytes_in_header_ != 1024 ) FatalErr( "Lookup table read failed." );
     bytes_in_header_ += read( fd_, &control2_[0], 1024 );
     SetK( control_[0] );
     nchunks_ = control_[1];
     unsigned int& ncontigs = control_[2];
     given_chunk_size_ = control_[3];
     given_chunk_overlap_ = control_[4];
     chunk_sizes_.resize(nchunks_);

     for ( int i = 0; i < (int) nchunks_; i++ )
       chunk_sizes_[i] = control_[ 5 + i ]; // this is size of loc_ of chunk i
     // note: chunk_sizes_ actually holds the sizes of locs_ for each chunk!!

     unsigned int max_chunk_size = Max(chunk_sizes_);
     bytes_in_header_ += read( fd_, &freq_[0], FourToK( ) * 4 );

     contig_name_.resize(ncontigs);
     contig_name_alt_.resize(ncontigs);
     contig_sizes_.resize(ncontigs);

  for ( int i = 0; i < (int) ncontigs; i++ )
    {
      bytes_in_header_ += read( fd_, &contig_sizes_[i], 4 );
      char cname[256];
      bytes_in_header_ += read( fd_, &cname[0], 252 );
      contig_name_[i] = String( &cname[0] );
      char cname_alt[256];
      bytes_in_header_ += read( fd_, &cname_alt[0], 256 );
      contig_name_alt_[i] = String( &cname_alt[0] );
    }

     // compute contig start positions for all contigs:
     // we have all contig sizes, so that we can calculate absolute
     // start positions with respect to the full (concatenated) reference
     // sequence.
     contig_start_.resize(ncontigs);
     unsigned int pos = 0;
     for ( int i = 0; i < (int) ncontigs; i++ )
     {    contig_start_[i] = pos;
          pos += contig_sizes_[i];    }
     locs_.reserve(max_chunk_size);

     // Initialize chunk_start_ vector:
     // we have enough info to compute start positions of all chunks right away
     chunk_start_.resize(nchunks_);
     off_t offset = bytes_in_header_;
     for ( unsigned int i=0; i < nchunks_; ++i ) {
       chunk_start_[i] = offset;
       offset += off_t(FourToK())*4;       // size of the lookup table
       offset += off_t(chunk_sizes_[i])*4; // size of the locs_ table
       offset += off_t( ( NBasesInChunk(i) + 15 ) / 16 )*4; // size of bases
     }
}

/** \brief Writes complete header of the lookup table file to disk.
 *
 *  This method (re)writes complete lookup table header on disk: control block
 *  1, control block 2, frequency table, and size/name/alt_name triplets
 *  for all contigs present in this lookup table (see the file format
 *  description). After first execution of this method, the
 *  header part of the file is "allocated" on disk and
 *  information that remains to be appended is
 *  chunks (i.e. lookup/location/bases triplets
 *  for each chunk).
 *
 *  Use this method only when building a new lookup table.
 *  This method forces (re)writing the current state of the lookup table
 *  header onto the disk and it \i does \i not perform any consistency checks
 *  (such as e.g.check for consistency between the lengths of the contig name
 *  and contig size lists and the number of contigs directly stored in the
 *  control block - the method may crash but leaving inconsistent data in the
 *  file is possible).
 *
 *  This method can be used only \i after variable-length
 *  part of the header was defined (i.e.the number of contigs is finalized
 *  via "adding" them all to the lookup table); otherwise the size of the
 *  allocated header on disk is wrong and appending chunks later will
 *  destroy the data. It is safe to re-write the header provided the
 *  number of contigs is set in advance and not changed between calls to
 *  this method.
 *
 * Currently, this is the only method that writes the header information
 * onto disk, in all its entirety (no partial updates of the header fields
 * on disk are available),
 * so actually this method \i has to be called multiple times: first time to
 * allocate the header in the file, and at least once after all chunks
 * are appended (writing chunks also updates chunk-related information
 * in the header). It's ok to call this method every time after a chunk
 * is appended for extra safety.
 *
 *  @see AddContigName()
 *  @see AddContigSize()
 *  @see AddContigStart()
 *  @see WriteChunk()
 */
void lookup_table::WriteHeader( ) const
{    lseek( fd_, 0, SEEK_SET );
     WriteBytes( fd_, &control_[0], 1024 );
     WriteBytes( fd_, &control2_[0], 1024 );
     WriteBytes( fd_, &freq_[0], FourToK( ) * 4 );
     vec<char> cname;
     for ( int i = 0; i < (int) contig_name_.size( ); i++ )
     {    WriteBytes( fd_, &contig_sizes_[i], 4 );
          cname.resize_and_set( 252, 0 );
          for ( int j = 0; j < (int) contig_name_[i].size( ); j++ )
               cname[j] = contig_name_[i][j];
          WriteBytes( fd_, &cname[0], 252 );
          cname.resize_and_set( 256, 0 );
          for ( int j = 0; j < (int) contig_name_alt_[i].size( ); j++ )
               cname[j] = contig_name_alt_[i][j];
          WriteBytes( fd_, &cname[0], 256 );
     }
}

/** \brief Appends chunk data to the file and updates the corresponding
 *  fields in the header (\i does \i not rewrite the header!!)
 *
 *  This is a dangerous unchecked direct access method, use with care.
 *  This method assumes that the lookup and location tables for this chunk
 *  are already initialized (in-memory) in the lookup table under
 *  construction. These tables are written to disk
 *  (appended to the file
 *  associated with this lookup table), followed by the bases that belong to
 *  this chunk (taken from the sequence passed as argument). The \c bases
 *  vector passed to this method should contain the \i complete, concatenated
 *  reference sequence, from which the chunks are being built; otherwise
 *  the chunk start positions stored in the file will be incorrect.
 *  The actual location of the chunk on the reference sequence \c bases is
 *  determined by the content of the pre-computed locations table used by this
 *  method. This method updates the chunk size- and
 *  position-related fields in the header
 *  (size of locations table for this chunk and counter of the number of
 *  chunks in control block 1; start of this chunk's basevector/start
 *  of chunk's bases in the basevector/number of bases in the chunk's
 *  basevector in control block 2 - see file format description), but
 *  does \i not write header info into the file on disk. Frequency table
 *  is \i not updated by this method. Use WriteHeader()
 *  to rewrite the updated header info onto the disk.
 *
 *  @see WriteHeader()
 *  @see MakeChunk()
 */

void lookup_table::WriteChunk( const vec<char>& bases )
{    WriteBytes( fd_, &lookup_[0], FourToK( ) * 4 );
     WriteBytes( fd_, &locs_[0], ( (longlong) NLocs( ) ) * LLCONST(4) );
     ForceAssert( control_[1] + 5 < 256 );
     control_[ control_[1] + 5 ] = locs_.size( );

     unsigned int start = Min(locs_);
     vec<char>::const_iterator itr(bases.begin()+start);
     vec<char>::const_iterator end(bases.begin()+Max(locs_)+K());
     basevector b(end-itr);
     for ( unsigned int idx = 0; itr != end; ++idx, ++itr )
     {
         char base = *itr;
         if ( !Base::isBase(base) ) base = 'A';
         b.Set(idx,Base::char2Val(base));
     }

     int byte_count = ( ( b.size( ) + 15 ) / 16 ) * 4;
     off_t fd_pos = lseek( fd_, 0, SEEK_CUR );
     int c2start = control_[1] * 4;
     ForceAssertLt( c2start, 252 ); // will fail if > 62 chunks
     control2_[ c2start ] = *( (unsigned int*) (&fd_pos) );
     control2_[ c2start + 1 ] = *( (unsigned int*) (&fd_pos) + 1 );
     control2_[ c2start + 2 ] = start;
     control2_[ c2start + 3 ] = b.size();
     BinaryWriteContent(fd_,b);
     ++control_[1];
}

/** \brief Organizes, in memory, lookup and locations data for the new chunk;
 *  does \i not update other chunk-related data, many accessors will not work
 *  yet.
 *
 *  This method takes pre-computed vector of (Kmer, location) pairs as
 *  its argument. Here 'location' is
 *  the position of the corresponding Kmer in the original base sequence;
 *  multiple pairs (e.g. multiple locations) for the same Kmer are allowed.
 *  As the result of this method's execution, 'lookup' and 'locations'
 *  data are (re)initialized according to the passed data, and the
 *  Kmer frequency data are updated (frequency of each Kmer is
 *  incremented by the number of its occurences in the chunk under
 *  construction). After execution of this method, the data are ready to use
 *  or to write the chunk onto disk (WriteChunk(), then WriteHeader()).
 *  NOTE: make sure that the same sequence
 *  that was used to calculate location offsets in index_loc is passed to
 *  WriteChunk() (do \i not truncate the sequence to the actual chunk size
 *  after offsets are calculated, WriteChunk() will take care of that;
 *  otherwise data written on disk will be inconsistent).
 *
 *  @see WriteChunk()
 *  @see WriteHeader()
 */
void
lookup_table::MakeChunk( vec< pair<unsigned int, unsigned int> >& index_loc )
{
  if ( index_loc.size( ) == 0 ) return;

  // Sort (kmer, offset) pairs by kmer then by offset
  sort( index_loc.begin( ), index_loc.end( ) );

  // Put the offsets in the locs_ vector.
  locs_.reserve( index_loc.size() );
  locs_.clear();
  for ( unsigned int r = 0; r < index_loc.size( ); r++ )
    locs_.push_back( index_loc[r].second );

  // Note we don't clear freq_ -- it contains frequency information
  // across the whole reference, not just the current chunk!  For
  // single-chunk lookup tables, we will construct lookup_ to maintain
  // the invariant that freq_[k] == lookup_[k+1]-lookup_[k].  For
  // multiple-chunk tables the expression on the right is the
  // frequency in the current chunk.

  // Initialize the lookup_ vector.
  const unsigned int undefined = 4000000000u;
  for ( unsigned int i = 0; i < four_to_K_; i++ )
    lookup_[i] = undefined;

  // Now fill in lookup_, walking backwards through index_loc so that
  // the last value into lookup[k] is the smallest one that has that
  // kmer number.  Also increment freq_ accordingly.
  for ( unsigned int r = index_loc.size( ) - 1; ; r-- ) {
    lookup_[index_loc[r].first] = r;
    ++freq_[index_loc[r].first];
    if ( r == 0 ) break;
  }
  // Fix up the lookup table for kmers not present in the index_loc.
  // Propagate lookup values backwards from the end so that the
  // invariant above, freq_[k] == lookup_[k+1]-lookup_[k], becomes
  // true for all k < 4^K.
  if ( undefined == lookup_[four_to_K_ - 1] )
    lookup_[four_to_K_ - 1] = locs_.size();
  for ( unsigned int i = four_to_K_ - 2; ; i-- ) {
    if ( undefined==lookup_[i] )
      lookup_[i] = lookup_[i+1];
    if ( i == 0 ) break;
  }
}

/** \brief Prepares chunk and writes it to disk.
 *
 *  Prepares all the data required to create a new chunk, writes chunk
 *  to disk and performs all the required updates in the lookup table
 *  header (but does \i not update the header on disk!): control blocks,
 *  frequency table.
 *
 *  The passed \c bases sequence can be longer than the actual chunk
 *  to be written: the actual content of the chunk is controlled
 *  by the first argument. Only the part of the sequence spanned
 *  by the offests present in \c index_loc will be written into the chunk.
 *
 *  NOTE: both arguments will be cleared after this method executes.
 *
 *  @param index_loc vector of (Kmer,offset) pairs, where offset is a
 *  position of the corresponding Kmer in the \c bases sequence (multiple
 *  pairs with the same Kmer, i.e. multiple positions are allowed).
 *  @param bases original sequence, in which the offsets of each Kmer
 *  (first argument) were calculated.
 *
 *  @see WriteHeader()
 */
void lookup_table::DumpChunk(
		       vec< pair<unsigned int, unsigned int> >& index_loc,
		       const vec<char>& bases )
{
  MakeChunk( index_loc );
  WriteChunk(bases);
  locs_.clear();
  index_loc.clear();
}


void lookup_table::FetchBasesFromDisk( unsigned int start, unsigned int stop,
     basevector& b )
{   ForceAssert(-1!=fd_);
    basevector b0;
    for ( unsigned int i = 0; i < NChunks( ); i++ )
     {    if ( start < StartBaseInChunk(i) ) continue;
          if ( stop > StopBaseInChunk(i) ) continue;
          unsigned int bases_offset = start - StartBaseInChunk(i);
          off_t d = DiskStartOfChunkBases(i) + off_t( bases_offset/4 );
          unsigned int before_bases = bases_offset - 4 * (bases_offset/4);
          b0.resize( stop - start + before_bases );

          // Note: the following was originally done with pread, but it did not
          // perform correctly on ia32 machines.

          off_t current = lseek( fd_, 0, SEEK_CUR );
          lseek( fd_, d, SEEK_SET );
          BinaryReadContent(fd_,b0);
          lseek( fd_, current, SEEK_SET );

          b.resize( stop - start );
          for ( unsigned int j = start; j < stop; j++ )
               b.Set( j - start, b0[ j - start + before_bases ] );
          return;    }
     FatalErr( "Internal error.  Requested start = " << start << " and stop = "
          << stop << " do not exist in a base chunk on disk." );
}

/** \breaf Reads in i-th chunk from disk.
 *
 *  This method reads in and makes available the i-th chunk data: its lookup,
 *  location offsets, and bases. Additionally, the chunk-dependent
 *  information is updated (first contig in chunk, stop positions of contigs
 *  in the chunk), so that all the accessor methods work correctly.
 */
void lookup_table::ReadChunk( int i )
{
  if (i==chunk_)
    return;  // nothing to do!
  ForceAssert(-1!=fd_);

  // For some reason the values of chunk_start_[i] are screwed up.  Therefore
  // we compute them from scratch here.

  int64_t chunk_start_i = chunk_start_[0];
  for ( int j = 0; j < i; j++ )
  {    chunk_start_i += FourToK( ) * 4;
       chunk_start_i += ( (longlong) chunk_sizes_[j] ) * LLCONST(4);
       chunk_start_i += BaseVec::physicalSize( NBasesInChunk(j) );    }

  if (i!=1+chunk_) {
    // Seek to right place
    lseek( fd_, chunk_start_i, SEEK_SET);
  }
  // Read lookup and locs
  read( fd_, &lookup_[0], FourToK( ) * 4 );
  // bad name. chunk_sizes_ actually holds sizes of locs_ for each chunk,
  // see ReadHeader(), WriteChunk() and file format specification...
  locs_.resize( chunk_sizes_[i] );
  read( fd_, &locs_[0], ( (longlong) chunk_sizes_[i] ) * LLCONST(4) );

  // Read bases for this chunk:
  b_start_ = StartBaseInChunk(i);
  unsigned int base_count = NBasesInChunk(i);
  //PRINT2(b_start_, base_count);
  int byte_count = ( ( base_count + 15 ) / 16 ) * 4;
  b_.Setsize(base_count);
  BinaryReadContent(fd_,b_);

  chunk_ = i; // remember what chunk is currently loaded

  // Set up ContigsInChunk info
  unsigned int pos, num_contigs_in_chunk;

  // translate absolute offset of chunk start and end ->
  // into indexes of contigs chunk start/end fall into and offsets of
  // start/end within those contigs:
  GetContigPos(b_start_, first_contig_in_chunk, pos);
  //PRINT2(first_contig_in_chunk, pos);

  GetContigPos(b_start_ + base_count - 1, num_contigs_in_chunk, pos);
  //PRINT2(num_contigs_in_chunk, pos);
  num_contigs_in_chunk -= first_contig_in_chunk;
  ++num_contigs_in_chunk;
  //PRINT2(first_contig_in_chunk, num_contigs_in_chunk);

  stop_base_of_contig_in_chunk.resize(num_contigs_in_chunk);
  for (unsigned int j=0; j<num_contigs_in_chunk-1; ++j) {
    stop_base_of_contig_in_chunk[j] = ContigStop(j+first_contig_in_chunk);
    //PRINT2(j, stop_base_of_contig_in_chunk[j]);
  }
  stop_base_of_contig_in_chunk[num_contigs_in_chunk-1] = b_start_ + base_count;
  //PRINT2(num_contigs_in_chunk-1, b_start_ + base_count);
  //cout << "Done reading chunk " << i << endl;
}

/** Builds a "basic" contig name of the form <name>[<id>], where
 *  <id> is the same as used in alternative name (<id>:<path_string>),
 *  and <name> is the <path_string> with directory path and file extension
 *  stripped out.
 */
String lookup_table::ContigNameBasic(int i)
{
     String alt = ContigNameAlt(i);
     String id = alt.Before( ":" ), rest = alt.After( ":" );
     if ( rest.Contains( "/" ) )
     {    int index;
          for ( index = (int) rest.size( ) - 1; index >= 0; index-- )
               if ( rest[index] == '/' ) break;
          rest = rest.substr( index + 1, rest.size( ) );
     }
     if ( rest.Contains( "." ) ) rest = rest.Before( "." );
     return rest + "[" + id + "]";
}


lookup_table::lookup_table( unsigned int K, const basevector &b ) :
  fd_(-1), // Don't ever try to read from or close fd
  control_(256,0), control2_(256,0),
  nchunks_(1), // Only one chunk in this style of table
  chunk_sizes_(1, b.size()),
  chunk_(0), // Chunk 0 is already loaded
  chunk_start_(0), first_contig_in_chunk(0), stop_base_of_contig_in_chunk(1,b.size())
{
  SetK(K);
  SetChunkParams(0, 0);
  InitializeFromContig("contig", b);
}

void lookup_table::InitializeFromContig( const String &name, const basevector &b )
{
  AddContigName(name, "", 0);
  AddContigStart(0);
  AddContigSize(b.size());
  b_ = b;
  b_start_ = 0;

  const unsigned int N = b.size();

  // index_loc is a list of (kmer number, offset in reference) pairs.
  // Fill it in from the bases.
  vec< pair<unsigned int, unsigned int> > index_loc;
  index_loc.reserve(N);

  KmerIndexSeq next(K_);
  next.Reset(b);
  const unsigned int last = b.size()-(K_-1);
  for (unsigned int i=0; i<last; ++i) {
    index_loc.push_back(make_pair(next(b), i));
  }

  // Make a lookup table chunk from the pairs
  MakeChunk( index_loc );
}


/// Implementation of high-level query interface
//////////////////////////////////////////////////////////////////////
/// Transform a vecbasevector into queries.  In this method,
/// a \i linear array of queries corresponding to all Kmers in all
/// the reads is created. The queries can be later distinguished from each other
/// and re-attributed to specific reads by their positions: for a query
/// returned by this method, its position is set to (read_position +
/// within_the_read_position), where read_position is the position of the
/// read this query corresponds to on a single base sequence made by
/// concatenation of all reads, and the within_the_read_position is the
/// offset within the read of the Kmer this query represents.
/// [NOTE: in this representation the total length of concatenated
/// \c bases should be less than 2^31 ].
///
/// If maxFreq is positive, don't put kmers with higher frequency
/// (or frequency 0) into the query set.
///
/// HACK: if maxFreq is negative, interpret it as -blockSize.  For
/// each block of the read, starting at the beginning, take the
/// lowest (nonzero) frequency kmer in that chunk, and put only
/// those kmers into the query set.
/// HACK2 : if <maxFreqDiscardRead> is <true>, there will be no queries
/// at all generated for a read that has at least one kmer with frequency
/// > maxFreq (i.e. the read will be effectively discarded)

void lookup_table::BasesToQueries(const vecbasevector &bases, vec<Query> &queries,
		      int maxFreq, unsigned int npasses,
				  int firstRead, int lastRead, bool maxFreqDiscardRead)
  {
    const unsigned int B = (maxFreq<0) ? -maxFreq : 0; // only used if maxFreq<0
    const unsigned int MF = (maxFreq>0) ? maxFreq : 0; // only used if maxFreq>0
    const unsigned int infinite=numeric_limits<unsigned int>::max();
                     // only used if maxFreq<0

    int v;
    unsigned int q=0;
    if (-1 == lastRead) lastRead = bases.size();
    for (v=firstRead; v < lastRead; ++v)
      q += max(0, bases[v].isize() - int(K_-1));
    queries.reserve(q);

    KmerIndexSeq next(K_);
    unsigned int pos = 0;
    basevector b;
    for (v=firstRead; v<lastRead; ++v) {
      if (bases[v].size()>=K_) {
	b = bases[v];
	for (unsigned int pass = 0; pass < npasses; ++pass) {
	  if (1==pass)
	    b.ReverseComplement();
	  const unsigned int last = b.size()-K_+1;
	  next.Reset(b);
	  if (maxFreq==0) {
	    // No need to test frequencies
	    for (unsigned int i=0; i<last; ++i) {
	      queries.push_back(Query(pos+i, next(b), pass>0));
	    }
	  } else if (maxFreq>0) {
	    if ( maxFreqDiscardRead ) {

	      unsigned int index, freq;
	      for (unsigned int i=0; i<last; ++i) {
		index = next(b);
		freq = freq_[index];
		if ( freq > MF ) {
		  vec<Query>::iterator from = queries.end();
		  from -= i;
		  queries.erase(from, queries.end());
		  break; // erase all queries accumulated for this read so far and abort: go get next read
		}
		queries.push_back(Query(pos+i, index, pass>0));
	      }

	    } else {

	      unsigned int index, freq;
	      for (unsigned int i=0; i<last; ++i) {
		index = next(b);
		freq = freq_[index];
		if (0 < freq && freq <= MF)
		  queries.push_back(Query(pos+i, index, pass>0));
	      }

	    }

	  } else { //HACK: find lowest nonzero freqs in blocks
	    // For pass 0 start at beginning, pass 1 start at end of read
	    int blockFirst = (0==pass) ? 0 : (max(B,last) - B);
	    int blockLast = (0==pass) ? min(B, last) : last;
	    unsigned int index, freq;
	    while ( 0 <= blockFirst && blockFirst < blockLast ) {
	      unsigned int minFreq = infinite, bestIndex=infinite, bestI=infinite;
	      for (int i=blockFirst; i<blockLast; ++i) {
		index = Index(b, i, K_);
		freq = freq_[index];
		if (0 < freq && freq < minFreq) {
		  minFreq = freq;
		  bestIndex = index;
		  bestI = i;
		}
	      }
	      if (minFreq<infinite)
		queries.push_back(Query(pos+bestI, bestIndex, pass>0));
	      if (0==pass) {
		blockFirst = blockLast;
		blockLast = min(B+blockLast, last);
	      } else {
		blockLast = blockFirst;
		blockFirst -= B;
	      }
	    }
	  }
	}
      }
      pos += bases[v].size(); // advance pos even if skipping over the read
    }
  }



