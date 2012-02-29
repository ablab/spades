///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "FastaFilestream.h"
#include "STLExtensions.h"
#include "system/System.h"

#include <strstream>

template <typename vecT, typename seqT, typename convT, typename verifT>
FastaFilestream<vecT,seqT,convT,verifT>::FastaFilestream( const String& filename,
                                                          FastaNameParser* name_parser )
  : filename_(filename), 
    preview_ptr_(0), 
    parsed_(false), 
    needs_pipe_(false), 
    name_parser_(name_parser)
{
  RequireRegularFile(filename_);
  String filename_plus( filename_ );
  if ( filename_plus.Contains( ".gz", -1 )  ||
       filename_plus.Contains( ".Z", -1 ) )
    needs_pipe_ = true;
  // it's compressed, so we'll need to open a pipe
}

template <typename vecT, typename seqT, typename convT, typename verifT>
FastaFilestream<vecT,seqT,convT,verifT>::FastaFilestream( const FastaFilestream& original )
  : filename_(original.filename_), 
    preview_ptr_(0), 
    parsed_(original.parsed_),
    needs_pipe_(original.needs_pipe_), 
    name_parser_(original.name_parser_)
{
  if ( original.preview_ptr_ != 0 )
    preview_ptr_ = new FastaFilestreamPreview( *(original.preview_ptr_ ) );
}


template <typename vecT, typename seqT, typename convT, typename verifT>
FastaFilestream<vecT,seqT,convT,verifT>& 
FastaFilestream<vecT,seqT,convT,verifT>::operator= ( const FastaFilestream& original )
{
  filename_ = original.filename_;
  preview_ptr_ = 0;
  parsed_ = original.parsed_;
  needs_pipe_ = original.needs_pipe_;
  name_parser_ = original.name_parser_;

  if ( original.preview_ptr_ != 0 )
    preview_ptr_ = new FastaFilestreamPreview( *(original.preview_ptr_ ) );
  
  return *this;
}  
  
template <typename vecT, typename seqT, typename convT, typename verifT>
int 
FastaFilestream<vecT,seqT,convT,verifT>::estimatedSize()
{
  if ( ! preview_ptr_ ) {
    istream* istream_ptr = getIstream_();
    preview_( istream_ptr );
    closeIstream_( istream_ptr );
  }
  return preview_ptr_->getSequenceSizes().size();
}


template <typename vecT, typename seqT, typename convT, typename verifT>
longlong
FastaFilestream<vecT,seqT,convT,verifT>::estimatedData()
{
  if ( ! preview_ptr_ ) {
    istream* istream_ptr = getIstream_();
    preview_( istream_ptr );
    closeIstream_( istream_ptr );
  }
  vec<streampos> seqSizes = preview_ptr_->getSequenceSizes();
  longlong total = accumulate( seqSizes.begin(), seqSizes.end(), (longlong) 0 );
  return total;
}


template <typename vecT, typename seqT, typename convT, typename verifT>
bool
FastaFilestream<vecT,seqT,convT,verifT>::verify()
{
  istream* istream_ptr = getIstream_();

  bool problem = false;

  longlong line_number = 1;

  const int buffer_size = 1024+1;
  char* buffer = new char[buffer_size];

  verifT verifier;

  while ( *istream_ptr && ! problem )
  {
    istream_ptr->getline( buffer, buffer_size );
    if ( ! verifier.verifyLine( buffer ) )
      problem = true;

    while ( istream_ptr->fail() && ! istream_ptr->eof() && ! problem)
    {
      istream_ptr->clear();
      istream_ptr->getline( buffer, buffer_size );
      if ( ! verifier.verifyRestOfLine( buffer ) )
	problem = true;
    }

    if ( ! problem )
      line_number++;
  }

  if ( problem )
    cout << "There was an illegal character on line " << line_number
	 << " of the file: " << endl << filename_ << endl;

  delete [] buffer;

  closeIstream_( istream_ptr );

  return ! problem;
}

template <typename vecT, typename seqT, typename convT, typename verifT>
void
FastaFilestream<vecT,seqT,convT,verifT>::parse_( vecString &names, 
                                                 vecT *p_sequences, 
                                                 const vec<int> *p_indices )
{
  bool verbose = false;

  if ( parsed_ )
    return;
  
  istream* istream_ptr = getIstream_();
  
  if ( ! preview_ptr_ ) {
    preview_( istream_ptr );
    resetIstream_( istream_ptr );
  }

  vec<streampos> & sequence_sizes = preview_ptr_->getSequenceSizes();
  long max_sequence_size = preview_ptr_->getMaxSequenceSize();
  streampos start_offset = preview_ptr_->getStartOffset();

  // Here we derive the start positions of each sequence in the file as 
  // an offset in bytes from the beginning.
  vec<streampos> sequence_positions;

  // The first sequence starts at start_offset.
  sequence_positions.push_back( start_offset );

  // We derive the positions of the remaining sequences by performing a
  // partial sum on their sizes and pushing the results onto the back
  // of sequence_positions.  

  // We need to add start_offset to the first sequence size for the
  // sum to come out right.
  if ( sequence_sizes.size() )
  {
      sequence_sizes[0] += start_offset;

      partial_sum( sequence_sizes.begin(), sequence_sizes.end(),
                   back_inserter( sequence_positions ) );

      // We take it back off because we use sequence_sizes later.
      sequence_sizes[0] -= start_offset;
  }

  char* buffer = new char[ max_sequence_size+1 ];

  convT converter( name_parser_ );

  String basename = getBasename_( filename_ );

  if ( verbose )
    cout << "Reading in data from " << basename << " (one . per 100K reads) " << flush;

  streampos curr_position = 0;

  vec<int> seq_indices;
  if ( p_indices )
    seq_indices = *p_indices;
  else
  {
    seq_indices.resize( sequence_sizes.size() );
    iota( seq_indices.begin(), seq_indices.end(), 0 );
  }

  const int chunk_size = 8192;
  char dumb_buffer[ chunk_size ];

  String name;
  seqT datum;

  // Go through the vector of the indices of the desired reads.
  for (unsigned int i = 0; i < seq_indices.size(); ++i ) {
    int seq_index = seq_indices[ i ];

    streampos offset_to_seq = sequence_positions[ seq_index ] - curr_position;

    // I realize this is tremendously ugly.  Unfortunately, a
    // pipe-based stream will apparently fail if you use either
    // tellg() or seekg(), so we have to advance through the file
    // by read()ing.  If someone figures out a way around this,
    // by all means, fix it.

    while ( offset_to_seq > 0 && istream_ptr->good() )
    {
      streampos bytes_to_be_read = min( (streampos) chunk_size, offset_to_seq );
      istream_ptr->read( dumb_buffer, bytes_to_be_read );
      int bytes_read = istream_ptr->gcount();
      offset_to_seq -= bytes_read;
      curr_position += bytes_read;
    }

    if ( ! istream_ptr->good() )
    {
      cout << "Read failed on " << filename_ << "." << endl;
      cout << "Exiting..." << endl;
      ForceAssert( istream_ptr->good() == true );
    }
    ForceAssertLe( sequence_sizes[seq_index], max_sequence_size + 1 );
    istream_ptr->read( buffer, sequence_sizes[ seq_index ] );

    int this_buffer_size = istream_ptr->gcount();
    curr_position += this_buffer_size;
    buffer[ this_buffer_size ] = 0;
    converter.extractAllFromBuffer( buffer, name, datum ); 

    if ( p_sequences )
      p_sequences->push_back( datum );

    //We use push_back_reserve here because the amount of space reserved
    //for the names was arbitrarily set at 50* number_of_names, and in some
    //cases names are longer than 50 on average.
    names.push_back_reserve( name );

    if ( verbose && (i+1 & 100000) == 0)
      cout << "." << flush;
  }

  if ( verbose )
    cout << ". done." << endl;

  delete [] buffer;

  parsed_ = 1;

  closeIstream_( istream_ptr );
}


template <typename vecT, typename seqT, typename convT, typename verifT>
const String 
FastaFilestream<vecT,seqT,convT,verifT>::getBasename_( const String& filename ) const
{
  int i;
  for (i = filename_.size() - 1; i >= 0; i--)
    if ( filename_[ i ] == '/' )
      break;
  
  return filename_.substr( i+1, filename_.size() - i );
}

template <typename vecT, typename seqT, typename convT, typename verifT>
void
FastaFilestream<vecT,seqT,convT,verifT>::preview_( istream* fasta_istream )
{
  String basename = getBasename_( filename_ );
  // cout << "Scanning file " << basename << "... " << flush;
  
  preview_ptr_ = new FastaFilestreamPreview( *fasta_istream );

  // cout << "found " << preview_ptr_->getSequenceSizes().size() << " reads." << endl;
}

template <typename vecT, typename seqT, typename convT, typename verifT>
istream* 
FastaFilestream<vecT,seqT,convT,verifT>::getIstream_() const
{
  istream* istream_ptr;

  if ( needs_pipe_ ) {
    // it's compressed, so we need to open a pipe
    string command = "gzip -dc " + filename_;
  
    procbuf* zcat_pipe = new procbuf( command.c_str(), ios::in );
    istream_ptr = new istream( zcat_pipe );
  }
  else {
    istream_ptr = new ifstream( filename_.c_str() );
  }
  
  //  istream_ptr->rdbuf()->allocate();
  
  return istream_ptr;
}

template <typename vecT, typename seqT, typename convT, typename verifT>
void
FastaFilestream<vecT,seqT,convT,verifT>::resetIstream_( istream*& fasta_istream ) const
{
  if ( needs_pipe_ ) {
    closeIstream_( fasta_istream );
    fasta_istream = getIstream_();
  }
  else {
    fasta_istream->clear();
    fasta_istream->seekg( 0 );
    fasta_istream->clear();
  }
}

template <typename vecT, typename seqT, typename convT, typename verifT>
void
FastaFilestream<vecT,seqT,convT,verifT>::closeIstream_( istream* fasta_istream ) const
{
  if ( needs_pipe_ ) {
    delete fasta_istream->rdbuf();
  }
  delete fasta_istream;
}

#define INSTANTIATE(aT,bT,cT,dT) \
template FastaFilestream<aT,bT,cT,dT>::FastaFilestream( const String&, FastaNameParser* ); \
template FastaFilestream<aT,bT,cT,dT>::FastaFilestream( const FastaFilestream& ); \
template FastaFilestream<aT,bT,cT,dT>& FastaFilestream<aT,bT,cT,dT>::operator= ( const FastaFilestream& ); \
template int FastaFilestream<aT,bT,cT,dT>::estimatedSize(void); \
template longlong FastaFilestream<aT,bT,cT,dT>::estimatedData(void); \
template void FastaFilestream<aT,bT,cT,dT>::parse_( vecString&, aT*, const vec<int>* ); \
template bool FastaFilestream<aT,bT,cT,dT>::verify();

INSTANTIATE( veccompseq, CompressedSequence, FastaSequenceConverter, FastaSequenceVerifier )

INSTANTIATE( veccharvector, charvector, FastaNullConverter, FastaNullVerifier )

INSTANTIATE( vecqualvector, qualvector, FastaQualityConverter, FastaQualityVerifier )

