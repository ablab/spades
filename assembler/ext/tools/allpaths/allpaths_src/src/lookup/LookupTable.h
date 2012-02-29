///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOOKUP_H
#define LOOKUP_H

#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "lookup/KmerIndex.h"
#include "lookup/Hit.h"
#include "TaskTimer.h"

#include <set>

typedef pair<unsigned int, unsigned int> LocSeq;

///
/// Class: lookup_table
///
/// Maintains a K-mer lookup table on disk.  The lookup
/// table consists of four parts (but see note on chunks, below):
///
///  - A frequency table freq_ of length 4^K, which contains the number of
/// times each kmer occurs in the reference.
///  - A list of indices lookup_ of length 4^k, which contains the index
/// into the offsets table for the first occurence of each kmer.
///  - A list of offsets locs_ of variable length, which is the
/// concatenation of the lists of offsets into the reference for
/// each kmer.
///  - The reference genome b_, concatenated into a single sequence, in
/// basevector format.
///
/// A reference position is represented by an *offset*, which is the
/// index in the concatenation of all the reference contigs.  This
/// allows positions in any reference of size < 2^32 bases to be
/// represented with a single unsigned int value (32 bits).
///
/// Note on chunks: Because a long reference makes a very large lookup
/// table (the human genome translates to a 13 Gb table), large
/// references are broken into overlapping chunks, and a separate
/// lookup_, locs_, and b_ are stored for each one.
///
/// Note: The freq_ table is for the entire reference, *not* the
/// current chunk.  If there is only one chunk then the freq_ table is
/// redundant, as you have the invariant freq_[k] ==
/// lookup_[k+1]-lookup_[k] for all k<4^K.  Memory locality will be
/// better if you use lookup_ rather than the freq_ to determine the
/// kmer frequencies in the current chunk.
///
/// Structure of lookup table on disk:
///
/// (begin example)
///
///     k (4 bytes)                                          |
///     number of chunks (4 bytes)                           |
///     number of contigs (4 bytes)                          | This
///     given chunk size (4 bytes)                           | is
///     given chunk overlap (4 bytes)                        | control
///     size of location chunk 1 in 4-byte words (4 bytes)   | block
///     size of location chunk 2 in 4-byte words (4 bytes)   | 1.
///     ...                                                  |
///     (pad to 1024 bytes)                                  |
///
///     start on disk of bases chunk 1 (8 bytes)             | This
///     start of bases in bases chunk 1 (4 bytes)            | is
///     number of bases in bases chunk 1 (4 bytes)           | control
///     ...                                                  | block
///     (pad to 1024 bytes)                                  | 2.
///
///     frequency table (4 * 4^k bytes)
///
///     contig 1 size (4 bytes)
///     contig 1 name (252 bytes)
///     contig 1 alt name (256 bytes)
///     ...
///
///     lookup to location chunk 1 (4 * 4^k bytes)
///     location chunk 1
///     bases chunk 1
///     ...
/// (end example)
///
/// Note that in the current design, there can be at most 64 chunks.
/// Otherwise control block 2 overflows.
class lookup_table {

public:

  //////////////////////////////////////////////////////////////////////
  ////// Interface for reading an existing lookup table
  //////////////////////////////////////////////////////////////////////

  //// The core interface for searching for hits in a lookup table.
  //// In overview, you construct the table from the file and then
  //// for each chunk i do the following:
  ////
  //// ReadChunk(i); then for each kmer k=Index(b,i,K) in the query
  //// you process the list of indices StartLocs(k) <= l <
  //// StopLocs(k).  These indices l map to offsets into the
  //// concatenated reference, Locs(l).  You may optionally test
  //// the overall frequency of the kmer by checking Freq(k), for
  //// instance to avoid kmers that are too common in the
  //// reference.
  //////////////////////////////////////////////////////////////////////

  /// This constructor is used for reading an existing lookup table
  /// from disk.  The file is closed by the destructor.
  lookup_table( const String & lookup_file ) :
    fd_(OpenForRead(lookup_file)),
    control_(256, 0),
    control2_( 256, 0 ),
    chunk_(-1)
  { ReadHeader(); }


  /// Copy constructor - takes an existing table and duplicates it, but opens
  /// a new file descriptor for completely independent access. Resets the current
  /// chunk to -1, so you must next load in a chunk to process
  /// The file is closed by the destructor.
  lookup_table( const lookup_table & table,  const String & lookup_file ) :
    control_(table.control_), control2_(table.control2_),
    K_(table.K_), four_to_K_(table.four_to_K_), Kmask_(table.Kmask_),
    nchunks_(table.nchunks_),
    given_chunk_size_(table.given_chunk_size_), given_chunk_overlap_(table.given_chunk_overlap_),
    chunk_sizes_(table.chunk_sizes_),
    contig_name_(table.contig_name_), contig_name_alt_(table.contig_name_alt_),
    contig_sizes_(table.contig_sizes_), contig_start_(table.contig_start_),
    freq_(table.freq_),  lookup_(table.lookup_),  locs_(table.locs_),
    b_(table.b_),  b_start_(table.b_start_),
    chunk_(table.chunk_),  chunk_start_(table.chunk_start_),
    first_contig_in_chunk(table.first_contig_in_chunk),
    stop_base_of_contig_in_chunk(table.stop_base_of_contig_in_chunk),
    ambiguous_reads_(table.ambiguous_reads_),
    start_read_(table.start_read_)
  {
    ForceAssertNe(table.fd_, -1);
    fd_ = OpenForRead(lookup_file);
    chunk_= -1;
    lseek( fd_, chunk_start_[0], SEEK_SET);
  }

  /// The K value.  You need this to compute kmer numbers using
  /// Index, the nonmember function in KmerIndex.h
  unsigned int K( ) const { return K_; }

  /// Number of chunks in lookup table.  The reference is broken
  /// into chunks to limit the total memory used.  Note chunk
  /// boundaries have nothing to do with contig boundaries.
  unsigned int NChunks( ) const { return chunk_sizes_.size( ); }

  /** \breaf Reads in i-th chunk from disk.
   *
   *  This method reads in and makes available the i-th chunk data: its lookup,
   *  location offsets, and bases. Additionally, the chunk-dependent
   *  information is updated (first contig in chunk, stop positions of contigs
   *  in the chunk), so that all the accessor methods work correctly.
   */
  void ReadChunk( int i );

  /// The number of times the kmer numbered index occurs in the reference
  /// (this is a global counter across all chunks!)
  inline unsigned int Freq( const unsigned int index ) const
    { return freq_[index]; }

  /// First index into the offsets table locs_ coorresponding to kmer l
  /// (for the currently loaded chunk only!).
  inline unsigned int Lookup( unsigned int l ) const { return lookup_[l]; }

  /// First index into the offsets table locs_ coorresponding to the
  /// specified kmer index (for currently loaded chunk only!).
  inline unsigned int StartLocs( unsigned int index ) const { return Lookup(index); }

  /// One-past-last index into the offsets table locs_ coorresponding to
  /// the specified kmer index (for currently loaded chunk only!).
  unsigned int StopLocs( unsigned int index ) const
  {
    unsigned int stop;
    if ( index < FourToK( ) - 1 ) stop = Lookup(index+1);
    else stop = NLocs( );
    return stop;
  }

  /// iterator into the sequence of locations (offsets) on the reference genome
  typedef vec<unsigned int>::const_iterator locs_iterator;

  /// Returns iterator pointing to the first location offset for the specified
  /// Kmer index (only offsets in the currently loaded chunk are accessible).
  /// NOTE: this is \e not an iterator into the reference sequence (bases)
  /// itself, the returned iterator steps over offset \e values (if any).
  locs_iterator StartLocsIterator(unsigned int index) const
    { return locs_.begin()+Lookup(index); }

  /// Returns iterator pointing to the last location offset value for the
  /// specified
  /// Kmer index (only offsets in the currently loaded chunk are accessible).
  /// NOTE: this is \e not an iterator into the reference sequence itself
  /// (bases), but rather into the collection of offset \e values.
  locs_iterator StopLocsIterator(unsigned int index) const {
    if ( index < FourToK()-1 ) return locs_.begin()+Lookup(index+1);
    return locs_.end();
  }

  /// The (start, stop) pair of indices into the offsets table locs_
  /// that coorrespond to kmer index.
  LocSeq LookupSeq( unsigned int index )
  { return LocSeq(StartLocs(index), StopLocs(index)); }

  /// Returns an offset into the concatenated reference.
  unsigned int Locs( unsigned int index ) const { return locs_[index]; }

  /** \brief Translate absolute base offset in the reference sequence into a
   *   contig number and base offset from the start of that contig.
   *
   *  @param [in] pos absolute offset in the concatenated reference sequence
   *  @param [out] c contig, which \c pos falls into
   *  @param [out] cpos offset within the conting \c c (contig start=0)
   *  corresponding to the absolute position \c pos.
   */
  void GetContigPos( unsigned int pos, unsigned int& c, unsigned int& cpos )
  {
    // Easier to understand but slower:
    //
    // for ( unsigned int i = 0; i < NContigs( ); i++ )
    // {    if ( pos >= ContigStart(i) && pos < ContigStop(i) )
    //      {    c = i;
    //           cpos = pos - ContigStart(i);
    //           return;    }    }

    ForceAssertGe( pos, contig_start_.front( ) );
    c = upper_bound( contig_start_.begin( ), contig_start_.end( ), pos )
      - contig_start_.begin( ) - 1;
    cpos = pos - ContigStart(c);
    if ( !( cpos < ContigSize(c) ) ) {
      PRINT3( pos, c, ContigStart(c) );
      PRINT2( BasesStart( ), BasesStop( ) );
      ForceAssertLt( cpos, ContigSize(c) );
    }
  }

  void GetContigPos( unsigned int pos, int & c, unsigned int& cpos )
  {
    ForceAssertGe( pos, contig_start_.front( ) );
    c = (int)(upper_bound( contig_start_.begin( ), contig_start_.end( ), pos )
	      - contig_start_.begin( ) - 1);
    cpos = pos - ContigStart(c);
    if ( !( cpos < ContigSize(c) ) ) {
      PRINT3( pos, c, ContigStart(c) );
      PRINT2( BasesStart( ), BasesStop( ) );
      ForceAssertLt( cpos, ContigSize(c) );
    }
  }


  /// Some higher-level code for working with the lookup table.
  //////////////////////////////////////////////////////////////////////

  /// Find alignments of the queries, parameterized by an operator
  /// that processes the hits that are found, allowing filtering
  /// in almost any way desired.
  ///
  /// Interface of HitReceiver: its operator() accepts hits, defined
  /// as a particular implied start position of a query on a target,
  /// that are specified as
  ///
  /// (unsigned int queryNumber, unsigned int queryPos, bool isRc,
  ///  unsigned int targetOffset, unsigned int targetContig,
  ///  unsigned int numberHits)
  ///
  /// and may do anything it wants with them.  The hits will be
  /// transmitted in clumps by query, in increasing order by query.
  /// The number of clumps per query is at most the number of chunks
  /// in the lookup table.  Within a clump for a single query the hits
  /// willl be sorted by targetOffset.  A particular queryPos =>
  /// targetOffset combination will not be repeated, except for hits
  /// occurring in the overlaps between chunks of a multi-chunk lookup
  /// table.  Each increasing-order run of hits will be terminated by
  /// a call to HitReceiver's ChunkDone() method.
  ///
  /// For meaning of maxFreq and maxFreqDiscardRead options see BasesToQueries,
  /// which they are forwarded to.
  template<typename HitReceiver>
  void FindHits(const vecbasevector &queries, HitReceiver &transmit,
		int maxFreq = 0, unsigned int npasses = 2,
		int firstRead=0, int lastRead = -1, bool maxFreqDiscardRead = false)
  {
    // TODO: Write an implementation of FindHits that avoids making
    // the q vector for single-chunk lookup tables.
    // TODO; potentially dangerous truncation of index by firstRead/lastRead
    vec<Query> q;
    BasesToQueries(queries, q, maxFreq, npasses, firstRead, lastRead, maxFreqDiscardRead);
    //if (queries.size()==1)
    //  copy(q.begin(), q.end(), ostream_iterator<Query>(cout, "\n"));
    FindHits(queries, q, transmit, firstRead, lastRead);
  }

  template<typename HitReceiver>
   void FindAmbiguousHits(const vecbasevector &queries, HitReceiver &transmit,
			 int maxFreq = 0, unsigned int npasses = 2,
			 int firstRead=0, int lastRead = -1)
  {
    vec<Query> q;
    BasesToQueries(queries, q, maxFreq, npasses, firstRead, lastRead,false);
    //if (queries.size()==1)
    //  copy(q.begin(), q.end(), ostream_iterator<Query>(cout, "\n"));
    FindAmbiguousHits(queries, q, transmit, firstRead, lastRead);
  }

  ////////////////////////////////////////////////////////////////////////
  //// Information about the reference genome, also stored in lookup table
  ////////////////////////////////////////////////////////////////////////

  //// Contig information...

  /// The number of contigs in table.
  unsigned int NContigs( ) const { return contig_name_.size( ); }

  /// Returns size (in bases) of the i-th contig
  unsigned int ContigSize( unsigned int i ) const { return contig_sizes_[i]; }

  /// Returns absolute start position (with respect to the full
  /// concatenated sequence of the reference genome) of the i-th contig
  unsigned int ContigStart( unsigned int i ) const { return contig_start_[i]; }

  /// Returns absolute offset (with respect to the full concatenated
  /// sequence of the reference genome) of the next base after the
  /// last base of the i-th contig
  unsigned int ContigStop( unsigned int i ) const
  {    return ContigStart(i) + ContigSize(i);    }

  /// Returns start position (absolute, with respect to the full
  /// concatenated sequence of the reference genome) of the last contig
  unsigned int LastContigStart( ) const { return contig_start_.back( ); }

  /// Returns name of the i-th contig
  String ContigName( int i ) const { return contig_name_[i]; }

  /// Returns alternative name of the i-th contig
  String ContigNameAlt( int i ) const { return contig_name_alt_[i]; }

  /** Builds a "basic" contig name of the form <name>[<id>], where
   *  <id> is the same as used in alternative name (<id>:<path_string>),
   *  and <name> is the <path_string> with directory path and file extension
   *  stripped out.
   */
  String ContigNameBasic( int i );

  /// Returns name of the last contig in this reference genome
  String LastContigName( ) const { return contig_name_.back( ); }

  /// Which contigs and bases are in the chunk...

  /// The number of contigs in currently-loaded chunk.
  unsigned int ContigsInChunk() { return stop_base_of_contig_in_chunk.size(); }

  /// The overall index of the first contig in currently-loaded chunk.
  unsigned int FirstContigInChunk() { return first_contig_in_chunk; }

  /// The stop position of contigOffset (that is, relative to first
  /// contig in chunk) in currently-loaded chunk.
  unsigned int StopBaseOfContigInChunk(unsigned int contigOffset)
  { return stop_base_of_contig_in_chunk[contigOffset]; }

  /// Whether this range of offsets falls in a single chunk.
  Bool CanFetchBasesFromDisk( unsigned int start, unsigned int stop )
  {
    for ( unsigned int i = 0; i < NChunks( ); i++ ) {
      if ( start < StartBaseInChunk(i) ) continue;
      if ( stop > StopBaseInChunk(i) ) continue;
      return True;
    }
    return False;
  }

  void FetchBasesFromDisk( unsigned int start, unsigned int stop,
			   basevector& b );

  /// The currently-available chunk of genome, all bases concatenated.
  const basevector& Bases( ) const { return b_; }

  /// The absolute start offset of Bases() in complete concatenated
  /// reference sequence.
  unsigned int BasesStart( ) const { return b_start_; }

  /// The absolute offset of one past the last base in Bases()
  /// (with respect to the complete concatenated reference sequence).
  unsigned int BasesStop( ) const { return b_start_ + b_.size( ); }

  /// Whether this offset is within the current chunk
  Bool BaseInMemory( unsigned int m )
  { return BasesStart( ) <= m && m < BasesStop( ); }

  /** \brief Translates an absolute offset (with respect to complete
   *  concatenated reference sequence)into the base (i.e. ACTG). \e UNCHECKED
   *  \e METHOD: if the specified offset is not in the currently loaded chunk,
   *  "array index out of bound" error will occur.
   */
  unsigned char Base( unsigned int m ) const
    { return b_[ m - BasesStart( ) ]; }

  //// Less vital header information
  //////////////////////////////////////////////////////////////////////

  /// The number of entries in the frequencies table -- mostly internal use
  unsigned int FourToK( ) const { return four_to_K_; }

  /// The number of kmers in the currently loaded chunk of reference genome.
  unsigned int NLocs( ) const { return locs_.size( ); }

  /// I don't think these accessors are used by any code, but the
  /// information is present in the header.
  unsigned int GivenChunkSize( ) { return given_chunk_size_; }
  unsigned int GivenChunkOverlap( ) { return given_chunk_overlap_; }

  /// The absolute offset (in the full concatenated reference sequence)
  /// of the first base in chunk c
  unsigned int StartBaseInChunk( int c ) { return control2_[ (c * 4) + 2 ]; }

  /// The number of bases stored in chunk c
  unsigned int NBasesInChunk( int c ) { return control2_[ (c * 4) + 3 ]; }

  /// Absolute offset (with respect to the complete concatenated reference sequence)
  /// of the one past the last base in this chunk
  unsigned int StopBaseInChunk( int c )
  {    return StartBaseInChunk(c) + NBasesInChunk(c);    }


  //////////////////////////////////////////////////////////////////////
  //// Interface for constructing a single-chunk lookup table in-memory
  //////////////////////////////////////////////////////////////////////

  /// Create a new lookup table from a basevector
  lookup_table( unsigned int K, const basevector &b );

  //////////////////////////////////////////////////////////////////////
  //// Interface for constructing a lookup table on disk
  //////////////////////////////////////////////////////////////////////

  /// This constructor is used for creating a new lookup table.
  lookup_table( int fd ) : chunk_(-1)
  {
    STATIC_ASSERT_M(sizeof(unsigned int) == 4 , bad_uint_size );
    STATIC_ASSERT_M( sizeof(off_t) == 8, bad_off_t_size );
    fd_ = fd;
    control_.resize( 256, 0 ), control2_.resize( 256, 0 );
  }

 /** \brief Sets the size of the lookup table's K-mer to \c k.
  *
  *  This method sets the K-mer size and automatically performs
  *  all required synchronization (internally maintained constants
  *  \c four_to_K, \c Kmask_ are updated; \c lookup and \c freq tables are
  *  resized - but \i not initialized, there is no data yet). Use this method
  *  only when creating new lookup tables and only prior to adding
  *  sequence data.
  */
  void SetK( unsigned int K );

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
  void AddContigName( String contig_name, const String& file_name,
		      int index_in_file );

  /** \brief Adds a new contig start position to the table of start positions.
   *
   *  This is a dangerous "direct access" method that leaves lookup table
   *  in an inconsistent state. Use it \i only along with other setter methods
   *  when building a table.
   *
   *  This method \i only adds a new start position, leaving all the
   *  other data structures in the lookup table unchanged. Use this method
   *  only along with other setter methods to ensure consistent state of
   *  the lookup table.
   *
   *  @see AddContigName()
   *  @see AddContigSize()
   *  @see SetChunkParams()
   *  @see DumpChunk()
   *  @see LookupTableBuilder::BuildTableFromContigs()
   */
  void AddContigStart( unsigned int start ) { contig_start_.push_back(start); }

  /** \brief Adds a new contig size to the table of sizes.
   *
   *  This is a dangerous "direct access" method that leaves lookup table
   *  in an inconsistent state. Use it \i only along with other setter methods
   *  when building a table.
   *
   *  This method \i only adds a new size, leaving all the
   *  other data structures in the lookup table unchanged. Use this method
   *  only along with other setter methods to ensure consistent state of
   *  the lookup table.
   *
   *  @see AddContigName()
   *  @see AddContigStart()
   *  @see SetChunkParams()
   *  @see DumpChunk()
   *  @see LookupTableBuilder::BuildTableFromContigs()
   */
  void AddContigSize( unsigned int count ) { contig_sizes_.push_back(count); }

  /** \brief Writes complete header of the lookup table file to disk.
   *
   *  This method (re)writes complete lookup table header on disk: control
   *  block 1, control block 2, frequency table, and size/name/alt_name
   *  triplets for all contigs present in this lookup table (see the file
   *  format description). After first execution of this method, the
   *  header part of the file is "allocated" on disk and
   *  information that remains to be appended is
   *  chunks (i.e. lookup/location/bases triplets
   *  for each chunk).
   *
   *  Use this method only when building a new lookup table.
   *  This method forces (re)writing the current state of the lookup table
   *  header onto the disk and it \i does \i not perform any consistency
   *  checks (such as e.g.check for consistency between the lengths of the
   *  contig name
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
   *  @see SetChunkParams()
   *  @see AddContigStart()
   *  @see WriteChunk()
   */
  void WriteHeader( ) const;

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
   *  determined by the content of the pre-computed locations table used by
   *  this method. This method updates the chunk size- and
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
  void WriteChunk( const vec<char>& bases );

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
  void DumpChunk( vec< pair<unsigned int, unsigned int> >& index_loc,
		  const vec<char>& bases );

  /** \brief Sets chunk size and overlap between chunks (both measured in
   *  bases) for ths lookup table.
   */
  void SetChunkParams( unsigned int chunk_size, unsigned int chunk_overlap )
  { control_[3] = chunk_size; control_[4] = chunk_overlap; }

  /////////////////////////////////////////////////////////////////////////////
  /// Destructor closes the open file fd_.
  /////////////////////////////////////////////////////////////////////////////
  ~lookup_table() { if (-1!=fd_) Close(fd_); }

  /// Other odd stuff
  /////////////////////////////////////////////////////////////////////////

  /** \brief Reads in the lookup table header from disk and initializes
   *  associated data structures (K, number, names and sizes of contigs,
   *  number of chunks, etc).
   *
   *  This method is intended for use in constructor only. The data stored in
   *  the file header are read in, and additional data maintained by the lookup
   *  table are initialized (such as chunk start offsets in the file,
   *  absolute positions of the contig starts in the full concatenated
   *  reference sequence, etc). All accessor
   *  methods defined in this class should return correct values after this
   *  method is executed.
   */
  void ReadHeader( );

  /// Make sure the lookup table is not holding any significant amount of memory
  void DestroyStuff( )
  {
    Destroy(freq_);
    Destroy(lookup_);
    Destroy(locs_);
    b_.Reinitialize( );
  }


  // we don't align ambiguous reads
  void SetAmbiguousReads(vec<int> &amb) { ambiguous_reads_ = amb; }
  void SetStartRead(int s) { start_read_=s; }

  //////////////////////////////////////////////////////////////////////
  /// Member variables
  //////////////////////////////////////////////////////////////////////
private:


  /// Control information
  int fd_;
  vec<unsigned int> control_, control2_;
  unsigned int K_, four_to_K_, Kmask_;
  unsigned int nchunks_;
  unsigned int given_chunk_size_, given_chunk_overlap_;
  vec<unsigned int> chunk_sizes_;
  vec<String> contig_name_, contig_name_alt_;
  vec<unsigned int> contig_sizes_, contig_start_;

  /// The key lookup-table structures
  vec<unsigned int> freq_;
  vec<unsigned int> lookup_;
  vec<unsigned int> locs_;

  /// Target bases
  basevector b_;
  unsigned int b_start_;

  /// Which chunk is currently in memory, which implies the file
  /// position is right for the next chunk.  Initialized to -1
  /// meaning no chunk in memory but ready to read chunk 0.
  int chunk_;
  vector<off_t> chunk_start_; ///< Where to seek to, for each chunk

  /// The first contig that is in the current chunk
  unsigned int first_contig_in_chunk;
  /// The stop offsets of the contigs in current chunk
  vec<unsigned int> stop_base_of_contig_in_chunk;

  vec<int> ambiguous_reads_;
  int start_read_;


  //////////////////////////////////////////////////////////////////////
  //// Private methods
  //////////////////////////////////////////////////////////////////////

  /// Turn a single basevector into a lookup table
  void InitializeFromContig( const String &name, const basevector &b );

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
   *  construction). After execution of this method, the chunk data are ready
   *  to use or to write the chunk onto disk (WriteChunk(),
   *  then WriteHeader()).
   *  NOTE: make sure that the same sequence
   *  that was used to calculate location offsets in index_loc is passed to
   *  WriteChunk() (do \i not truncate the sequence to the actual chunk size
   *  after offsets are calculated, WriteChunk() will take care of that;
   *  otherwise data written on disk will be inconsistent).
   *
   *  @see WriteChunk()
   *  @see WriteHeader()
   */
  void MakeChunk( vec< pair<unsigned int, unsigned int> >& index_loc );

  /// Internal use: where the chunk bases start on disk
  off_t DiskStartOfChunkBases( int c )
  {
    off_t answer;
    unsigned int* answerp = (unsigned int*) &answer;
    answerp[0] = control2_[ c * 4 ];
    answerp[1] = control2_[ (c * 4) + 1 ];
    return answer;
  }

  /// Implementation of high-level query interface
  //////////////////////////////////////////////////////////////////////
  /// Transform a vecbasevector into queries.  In this method,
  /// a \i linear array of queries corresponding to Kmers in all
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
  void BasesToQueries(const vecbasevector &bases, vec<Query> &queries,
		      int maxFreq, unsigned int npasses,
		      int firstRead=0, int lastRead=-1, bool maxFreqDiscardRead = false);

  /// Find all hits in the current chunk for the interval of query
  /// kmers.  Returned hits, expressed as a vector of RawHits, are in
  /// order of query kmers.  posOffset is subtracted from the
  /// QueryPos's indicated in the query kmers, supporting
  /// concatenation of multiple queries.
  void FindHits(vec<Query>::const_iterator query,
		vec<Query>::const_iterator lastQuery,
		int posOffset,
		vec<RawHit> &hits) {
    hits.clear();
    unsigned int j;
    LocSeq locseq;
    for ( ; query != lastQuery; ++query) {
      locseq = LookupSeq(query->Kmer());
      for (j=locseq.first; j!=locseq.second; ++j) {
	hits.push_back(RawHit(locs_[j], query->QueryPos()-posOffset, query->IsRc()));
      }
    }
  }

  /// Find alignments of the queries, parameterized by an operator
  /// that processes the hits that are found, allowing filtering in
  /// almost any way desired. Interface of HitReceiver described in
  /// FindHits above.
  template<typename HitReceiver>
  void FindHits(const vecbasevector &bases, ///< Used for lengths only!
		const vec<Query> &queries, ///< Sorted by query position.
		HitReceiver &transmit,
		int firstRead=0, int lastRead=-1)
  {
    if (-1 == lastRead) lastRead = bases.size();
    CompareQueriesByPos compareQueriesByPos;
    CompareRawHitsByOffset compareRawHitsByOffset;
    CompareRawHitsByQueryStart compareRawHitsByQueryStart;
    CompareRawHitsByQueryStartOffset compareRawHitsByQueryStartOffset;
    vec<RawHit> hits;
    vec<Query>::const_iterator firstQuery, lastQuery;
    vec<RawHit>::iterator firstRcHit, firstHit, lastHit;
    vec<RawHit>::iterator it, last, runLast;
    unsigned int j, startpos, stoppos, contig, pos, firstcontig, ncontigs;
    longlong startOnTarget;
    for (unsigned int chunk=0; chunk<nchunks_; ++chunk) {
      ReadChunk(chunk);
      ncontigs = ContigsInChunk();
      firstcontig = FirstContigInChunk();
      startpos=0;
      firstQuery = queries.begin();
      for (int i=firstRead; i< lastRead; ++i, startpos=stoppos) {

	// Determine which of the queries are from this sequence
	stoppos = startpos + bases[i].size();
	lastQuery = lower_bound(firstQuery, queries.end(), Query(stoppos), compareQueriesByPos);

	// don't bother with reads declared ambiguous
	if ( !ambiguous_reads_.empty() &&
	     BinMember(ambiguous_reads_,(i+start_read_)) ) {
	  firstQuery=lastQuery;
	  continue;
	}

	FindHits(firstQuery, lastQuery, startpos, hits); // get the RawHits

	// separate into fw, rc and process each independently as [firstHit, lastHit)
	firstRcHit = partition(hits.begin(), hits.end(), RawHitIsFw());
	firstHit = hits.begin();
	lastHit = firstRcHit;
	while (firstHit != hits.end()) {
	  sort(firstHit, lastHit, compareRawHitsByQueryStartOffset);
	  it = firstHit;
	  for (j=0; j<ncontigs; ++j) {
	    // Separate hits by target contig, creating subseq [it, last) of hits.
	    if (j+1 < ncontigs) { // More than one contig remains: find junction
	      // As an optimization, first locate where the junction
	      // falls in the hits sorted by start position.  The
	      // start position is no larger than the offset, so there
	      // can't be any hits for this contig to the right of
	      // that point.  As an optimization, skip to next contig
	      // if there isn't at least one hit for this one.
	      RawHit contigEnd(StopBaseOfContigInChunk(j), 0, true);
	      if (it->QueryStartOnTarget() > contigEnd.Offset()) continue;
	      last = lower_bound(it, lastHit, contigEnd,
				 compareRawHitsByQueryStartOffset);
	      // Then partition those by their actual offsets.
	      last = stable_partition(it, last, RawHitOffsetIsBefore(contigEnd));
	    } else { // Only one contig left, so all hits belong to it
	      last = lastHit;
	    }
	    // Now transmit each run of hits with same query start in this contig group.
	    for ( ; it!=last; it = runLast) {
	      // one past end of the *it run
	      runLast = it;
	      startOnTarget = it->QueryStartOnTarget();
	      do
		++runLast;
	      while (runLast!=last && runLast->QueryStartOnTarget()==startOnTarget);
	      transmit(i, it->QueryPos(), it->IsRc(), it->Offset(), firstcontig+j,
		       distance(it, runLast));
	    }
	  }
	  // Now advance to rc part of hits, if any
	  firstHit = lastHit;
	  lastHit = hits.end();
	}
	firstQuery = lastQuery;
      }
      transmit.ChunkDone();
    }
  }


template<typename HitReceiver>
void FindAmbiguousHits(const vecbasevector &bases, ///< Used for lengths only!
		  const vec<Query> &queries, ///< Sorted by query position.
		  HitReceiver &transmit,
		  int firstRead=0, int lastRead=-1)
{
  if (-1 == lastRead) lastRead = bases.size();
   CompareQueriesByPos compareQueriesByPos;
  vec<Query>::const_iterator firstQuery, lastQuery;
  unsigned int j, startpos, stoppos, contig, pos, firstcontig, ncontigs;
  longlong startOnTarget;


  vec<bool> done_ovlp(NChunks(),false);

  for (unsigned int chunk=0; chunk<nchunks_; ++chunk)
  {
    ReadChunk(chunk);
    ncontigs = ContigsInChunk();
    firstcontig = FirstContigInChunk();
    startpos=0;
    firstQuery = queries.begin();

    for (int i=firstRead; i< lastRead; ++i, startpos=stoppos)
    {
       // Determine which of the queries are from this sequence
      stoppos = startpos + bases[i].size();
      lastQuery = lower_bound(firstQuery, queries.end(), Query(stoppos), compareQueriesByPos);

      if ( transmit.IsAmbiguous(i) ) {
	firstQuery=lastQuery;
	continue;
      }

      set<longlong> startPos;
      for ( ; firstQuery != lastQuery; ++firstQuery )
      {
	LocSeq locseq = LookupSeq(firstQuery->Kmer());
	unsigned int k(0);

// 	cout << "\t\t#matches this query: " << distance(queries.begin(),firstQuery) <<" "
// 	     <<  locseq.second-locseq.first << " "<<startPos.size() << endl;

	for ( k=locseq.first; k != locseq.second; ++k )
	{
	  if ( transmit.IsAmbiguous(i) )
	    break;

	  RawHit h(locs_[k], firstQuery->QueryPos()-startpos, firstQuery->IsRc());

	  // only transmit each start position one time
	  longlong min_pos(std::max(longlong(0),h.QueryStartOnTarget()-transmit.Bandwidth()));
	  longlong max_pos(h.QueryStartOnTarget()+transmit.Bandwidth());
	  for ( ; min_pos <= max_pos; ++min_pos )
	    if ( ! startPos.insert(min_pos).second )
	      break;

	  if ( min_pos != max_pos+1 )
	    continue;

	  // find this contig

	  unsigned int j(0);
	  while ( j<ncontigs &&
		  h.Offset() >= StopBaseOfContigInChunk(j) ) {
	    ++j;
	  }

// 	  cout << "transmitting: "
// 	       << i <<" "
// 	       << h.QueryStartOnTarget() <<" "
// 	       << firstcontig+j <<" "
// 	       << h.QueryPos() << " "
// 	       << h.Offset() << endl;

	  if ( chunk>0 &&
	       locs_[k] < (StartBaseInChunk(chunk)+GivenChunkOverlap()) &&
	       done_ovlp[chunk-1] ) {
	    continue;
	  }

	  // Now transmit each run of hits with same query start in this contig group.
	  transmit(i, h.QueryPos(), h.IsRc(), h.Offset(), firstcontig+j,1);
	}

	if ( transmit.IsAmbiguous(i) ) {
	  firstQuery=lastQuery;
	  break;
	}
     }  //queries
    } // reads

    done_ovlp[chunk]=true;
  } // chunks

}



  /// Deliberately PRIVATE UNIMPLEMENTED copy constructor.  Don't
  /// copy lookup tables.
  lookup_table(const lookup_table &);

  /// Deliberately PRIVATE UNIMPLEMENTED assignment operator.
  /// Don't assign lookup tables.
  void operator=(const lookup_table &);
};


inline Bool AmbiguousBase( char base )
{ return !Base::isBase(base); }


#endif

// Synonyms: Various synonyms
//   lookup table - See <lookup_table>
