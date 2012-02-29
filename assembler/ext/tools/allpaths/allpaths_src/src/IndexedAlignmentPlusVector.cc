///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "IndexedAlignmentPlusVector.h"
#include "VecAlignmentPlus.h"

// These strings are at the head of the vector and index files for verification purposes.

const String VecAlignmentPlusHeaderIO::mStrVectorHeader( "alignment_plus vector V.2\n" );
const String VecAlignmentPlusHeaderIO::mStrIndexHeader ( "alignment_plus index V.2\n" );

// A vector file contains the first header followed by the number of
// alignments in raw integer form, followed by the alignments
// themselves in a format described below, sorted by id1:

// An index file contains the second header, followed by an entry for
// each id1 ranging from 0 to the maximum id1 of any alignment in the
// vector.  Each entry contains the offset in the vector file of the
// start of the first alignment with that id1, followed by the number
// of alignments with that id1.

inline
int SizeOfIndexEntry() { return sizeof(longlong) + sizeof(int); }

// New alignment format:
//
// A header byte contains a set of flags indicating whether certain
// portions of the alignment have been compressed.
// 
// header byte:
//   bit 0 (0x01): use previous id1
//   bit 1 (0x02): second sequence is rc
//   bit 2 (0x04): pos1, pos2, errors are compressed
//   bit 3 (0x08): block 0 is compressed
//   bit 4 (0x10): block 1 is compressed
//   bit 5 (0x20): block 2 is compressed
//   bit 6 (0x40): block 3 is compressed
//   bit 7 (0x80): block 4 is compressed

const unsigned char HB_use_previous_id1  = 0x01;
const unsigned char HB_rc2               = 0x02;
const unsigned char HB_ppe_is_packed     = 0x04;
const unsigned char HB_block_is_packed[] = { 0x08, 0x10, 0x20, 0x40, 0x80 };
const unsigned char HB_all_blocks_packed = 0xF8;
const unsigned char HB_all_flags_set     = 0xFF;

const int HB_blocks_in_header = 5;

// If the "use previous id1" flag is set, then id1 is the id1 of the
// previous align and id2 is the next int.  Otherwise, id1 is the next
// int and id2 is the int after that.
//
// If the "pos1, pos2, errors are compressed flag is set, then the
// next int contains these values at these offsets:

const int PPE_pos1_packable_bits         = 0x000007FF; // 11 bits
const int PPE_pos1_packable_bits_short   =     0x07FF; // 11 bits
const int PPE_pos1_unpackable_bits       = ~PPE_pos1_packable_bits;
const int PPE_pos1_shift                 = 0;

const int PPE_pos2_packable_bits         = 0x000007FF; // 11 bits
const int PPE_pos2_packable_bits_short   =     0x07FF; // 11 bits
const int PPE_pos2_unpackable_bits       = ~PPE_pos2_packable_bits;
const int PPE_pos2_shift                 = 11;

const int PPE_errors_packable_bits       = 0x000003FF; // 10 bits
const int PPE_errors_packable_bits_short =     0x03FF; // 10 bits
const int PPE_errors_unpackable_bits     = ~PPE_errors_packable_bits;
const int PPE_errors_shift               = 22;

// otherwise, the next int contains pos1, the next pos2, and the next errors.
//
// This is followed by the score of the align in float form.
//
// Next comes the number of blocks in int form.
// 
// Each block can be either uncompressed:
//
//   gap    (signed int)
//   length (unsigned int)
//
// or compressed:

const int   BLOCK_len_packable_bits       = 0x000007FF; // 11 bits
const short BLOCK_len_packable_bits_short =     0x07FF;
const int   BLOCK_len_unpackable_bits     = ~BLOCK_len_packable_bits;
const int   BLOCK_len_shift               = 0;

const int BLOCK_gap_packable_bits         = 0x0000001F; // 5 bits (4 digits, 1 sign)
const int BLOCK_gap_sign_bit              = 0x00000010; // sign bit
const int BLOCK_gap_packable_mag_bits     = 0x0000000F; // magitude bits
const int BLOCK_gap_sign_bits             = ~BLOCK_gap_packable_mag_bits;
const int BLOCK_gap_shift                 = 11;

const int BLOCK_max_size                  = 2 * sizeof(int);

// The flags for the first five blocks are found in the header block.
// If there are more than five blocks, the flags for each set of eight
// blocks following the first five will be in another byte of the
// form:
//
// block compression byte:
//   bit 0 (0x01): block n+1 is compressed
//   bit 1 (0x02): block n+2 is compressed
//   bit 2 (0x04): block n+3 is compressed
//   bit 3 (0x08): block n+4 is compressed
//   bit 4 (0x10): block n+5 is compressed
//   bit 5 (0x20): block n+6 is compressed
//   bit 6 (0x40): block n+7 is compressed
//   bit 7 (0x80): block n+8 is compressed
//
// (where n was the last block read).  This byte will precede the next
// (up to) 8 alignments.

const unsigned char BCB_block_is_packed[] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };
const unsigned char BCB_all_blocks_packed = 0xFF;

const int maxHeaderSize = 1 +                  // header byte
                          2 * sizeof( int ) +  // id1, id2
                          3 * sizeof( int ) +  // pos1, pos2, errors
                          sizeof( float ) +    // score
                          sizeof( int );       // numBlocks

const int maxOutputChunkSize = max( maxHeaderSize,
                                    1 +                   // block compression byte
                                    8 * BLOCK_max_size ); // 8 big blocks


const int gPageSize = 8192;


ReadBuffer::ReadBuffer( const int bufferSize )
  : mBufferSize( max(bufferSize,2*gPageSize) ),
    mpBufferStart( new unsigned char[ mBufferSize ] ),
    mpBufferEnd( mpBufferStart + mBufferSize ),
    mpBufferCursor( mpBufferEnd ),
    mBufferIsStale( true )
{
}  

ReadBuffer::~ReadBuffer()
{
  delete [] mpBufferStart;
}

inline
void
ReadBuffer::CheckBuffer( const int amountToBuffer )
{
  int bytesLeft = mpBufferEnd - mpBufferCursor;

  if ( bytesLeft < amountToBuffer )
  {
    if ( bytesLeft < 0 )
    {
      cout << "Unexpected end-of-file." << endl;
      TracebackThisProcess( );
    }

    if ( mpIstream->good() )
    {
      if ( amountToBuffer > mBufferSize )
      {
        unsigned char *pNewBufferStart = new unsigned char[amountToBuffer];
     
        if ( bytesLeft > 0 )
          memmove( pNewBufferStart, mpBufferCursor, bytesLeft );
     
        delete [] mpBufferStart;
        mpBufferStart = pNewBufferStart;
      }
      else
      {
        if ( bytesLeft > 0 )
          memmove( mpBufferStart, mpBufferCursor, bytesLeft );
      }

      mpBufferCursor = mpBufferStart;
      mpBufferEnd = mpBufferStart + bytesLeft;

      streamoff bytesToRead = mBufferSize - ( mpBufferEnd - mpBufferStart );
      streampos pageEnd = mpIstream->tellg() + bytesToRead;
      pageEnd = pageEnd / gPageSize * gPageSize;  // always end on a block edge
      bytesToRead = pageEnd - mpIstream->tellg();

      mpIstream->read( (char *) mpBufferEnd, bytesToRead );
      mpBufferEnd += mpIstream->gcount();

      mBufferIsStale = false;
    }
  }
}

inline
void
ReadBuffer::JumpTo( streampos filePos )
{
  // if we've hit the end of the file, we need to reset the eof flag
  // before getting the file position
  if ( mpIstream->eof() )
    mpIstream->clear();

  streampos currPos = mpIstream->tellg();
  
  // if the desired position is somewhere in the current buffer and
  // the buffer is not stale (i.e. the current position of the
  // filestream actually corresponds to the end of the buffer).
  if ( ! mBufferIsStale && 
       filePos < currPos &&
       currPos - filePos <= mpBufferEnd - mpBufferStart )
  {
    mpBufferCursor = mpBufferEnd - ( currPos - filePos );
  }
  // otherwise, jump to the right spot and mark the buffer as stale
  else
  {
    // empty the buffer
    mpBufferCursor = mpBufferEnd;
    // if we're at the end of the file, we need to reset the status of
    // the stream so that it seeks properly
    bool atEndOfFile = mpIstream->eof();
    if ( atEndOfFile )
      mpIstream->clear();
    // seek to the correct position
    mpIstream->seekg( filePos, ios::beg );
    // if we were at the end of the file, the eof bit will be
    // spuriously set by the seekg, so we clear it
    if ( atEndOfFile )
      mpIstream->clear();
    // mark the buffer as stale
    mBufferIsStale = true;
  }
}




AlignmentPlusReader::AlignmentPlusReader( ReadBuffer *pReadBuffer )
  : mpReadBuffer( pReadBuffer ),
    mLastId1( -1 ),
    numBlocks( 0 )
{
  gaps.Setsize( 65535 );
  lens.Setsize( 65535 );
}

VecAlignmentPlusReader::VecAlignmentPlusReader( const String &strVectorFilename,
                                                const int bufferSize )
  : mStrVectorFilename( strVectorFilename ),
    mStrIndexFilename( strVectorFilename + "_index" ),
    mpDataIn( 0 ),
    mpIndexIn( 0 ),
    mDataPtrIsShared( false ),
    mNumAligns( -1 ),
    mRequestedId1Start( 0 ),
    mNumAlignsWithRequestedId1( 0 ),
    mReadBuffer( bufferSize ),
    mReader( &mReadBuffer )
{
}

VecAlignmentPlusReader::VecAlignmentPlusReader( istream *pVectorIstream,
                                                const int bufferSize )
  : mStrVectorFilename( "" ),
    mStrIndexFilename( "" ),
    mpDataIn( pVectorIstream ),
    mpIndexIn( 0 ),
    mDataPtrIsShared( true ),
    mNumAligns( -1 ),
    mRequestedId1Start( 0 ),
    mNumAlignsWithRequestedId1( 0 ),
    mReadBuffer( bufferSize ),
    mReader( &mReadBuffer )
{
  mReadBuffer.Attach( mpDataIn );

  if ( ! mHeaderIO.VerifyDataHeader( mpDataIn ) )
  {
    cout << "File may be in incorrect format." << endl;
    TracebackThisProcess();
  }
}

VecAlignmentPlusWriter::VecAlignmentPlusWriter( const String &strVectorFilename,
                                                const bool createIndex,
                                                const int bufferSize )
  : mStrVectorFilename( strVectorFilename ),
    mStrIndexFilename( createIndex ? strVectorFilename + "_index" : "" ),
    mpDataOut( 0 ),
    mpIndexOut( 0 ),
    mNumAligns( 0 ),
    mLastId1( -1 ),
    mLastId1Start( 0 ),
    mNumAlignsWithLastId1( 0 ),
    numBlocks( 0 ),
    mBufferSize( max(bufferSize,2*gPageSize ) ),
    mpBufferStart( new unsigned char[ mBufferSize ] ),
    mpBufferEnd( mpBufferStart + mBufferSize ),
    mpBufferCursor( mpBufferStart )
{
  gaps.Setsize( 65535 );
  lens.Setsize( 65535 );
}

VecAlignmentPlusWriter::VecAlignmentPlusWriter( ostream *pVectorOstream,
                                                const int bufferSize )
  : mStrVectorFilename( "" ),
    mStrIndexFilename( "" ),
    mpDataOut( pVectorOstream ),
    mpIndexOut( 0 ),
    mDataPtrIsShared( true ),
    mNumAligns( 0 ),
    mLastId1( -1 ),
    mLastId1Start( 0 ),
    mNumAlignsWithLastId1( 0 ),
    numBlocks( 0 ),
    mBufferSize( max(bufferSize,2*gPageSize ) ),
    mpBufferStart( new unsigned char[ mBufferSize ] ),
    mpBufferEnd( mpBufferStart + mBufferSize ),
    mpBufferCursor( mpBufferStart )
{
  mHeaderIO.WriteDataHeader( mpDataOut );

  gaps.Setsize( 65535 );
  lens.Setsize( 65535 );
}

VecAlignmentPlusReader::~VecAlignmentPlusReader()
{
  if ( ! mDataPtrIsShared )
    delete mpDataIn;

  delete mpIndexIn;
}

VecAlignmentPlusWriter::~VecAlignmentPlusWriter()
{
  delete [] mpBufferStart;
  if ( ! mDataPtrIsShared )
    delete mpDataOut;
  delete mpIndexOut;
}

inline
void
VecAlignmentPlusWriter::OpenDataFileForWriting()
{
  if ( ! mpDataOut )
  {
    // This open mode truncates an existing file.
    mpDataOut = new ofstream( mStrVectorFilename.c_str(), ios::out | ios::binary );

    if ( ! mpDataOut->good() )
    {
      cout << "Unable to open alignment vector file " 
           << mStrVectorFilename << " for writing." << endl;
      TracebackThisProcess();
    }

    mHeaderIO.WriteDataHeader( mpDataOut );

    // Empty the write buffer.
    mpBufferCursor = mpBufferStart;
  }
}

inline
void
VecAlignmentPlusWriter::OpenIndexFileForWriting()
{
  if ( ! mpIndexOut )
  {
    // This open mode truncates an existing file.
    mpIndexOut = new ofstream( mStrIndexFilename.c_str(), ios::out | ios::binary );

    if ( ! mpIndexOut->good() )
    {
      cout << "Unable to open alignment index file " 
           << mStrIndexFilename << " for writing." << endl;
      TracebackThisProcess();
    }

    mHeaderIO.WriteIndexHeader( mpIndexOut );
  }
}


inline
void
VecAlignmentPlusWriter::OpenDataFileForAppending()
{
  if ( ! mpDataOut )
  {
    if ( IsRegularFile( mStrVectorFilename.c_str() ) )
    {
    // This open mode only allows appending to an existing file.
      mpDataOut = new ofstream( mStrVectorFilename.c_str(), ios::out | ios::app | ios::binary );

      if ( ! mpDataOut->good() )
      {
        cout << "Unable to open alignment vector file " 
             << mStrVectorFilename << " for writing." << endl;
        TracebackThisProcess();
      }

      // Empty the write buffer.
      mpBufferCursor = mpBufferStart;
    }
    else
      OpenDataFileForWriting();
  }
}

inline
void
VecAlignmentPlusWriter::OpenIndexFileForAppending()
{
  if ( ! mpIndexOut )
  {
    if ( IsRegularFile( mStrIndexFilename.c_str() ) )
    {
      // This open mode allows repositioning and subsequent writing
      // within an existing file.  We may need to overwrite the last
      // entry in the index, so we need this ability.
      mpIndexOut = new ofstream( mStrIndexFilename.c_str(), ios::in | ios::out | ios::binary );

      if ( ! mpIndexOut->good() )
      {
        cout << "Unable to open alignment index file " 
             << mStrIndexFilename << " for writing." << endl;
        TracebackThisProcess();
      }
    }
    else
      OpenIndexFileForWriting();
  }
}

void
VecAlignmentPlusWriter::LoadLastId1FromIndex()
{
  if ( IsRegularFile( mStrIndexFilename.c_str() ) )
  {
    ifstream indexIn( mStrIndexFilename.c_str(), ios::in | ios::binary );
    
    if ( ! indexIn.good() ||
         ! mHeaderIO.VerifyIndexHeader( &indexIn ) )
    {
      cout << "Unrecognized alignment index file format in: " << endl;
      cout << mStrIndexFilename << endl;
      cout << "Unable to open alignment index file " 
           << mStrIndexFilename << " for appending." << endl;
      TracebackThisProcess();
    }
 
    indexIn.seekg( 0, ios::end );
    int indexSize = indexIn.tellg();
    
    indexSize -= mHeaderIO.GetIndexHeaderSize();
    mLastId1 = indexSize / SizeOfIndexEntry() - 1;

    if ( mLastId1 >= 0 )
    {
      indexIn.seekg( -SizeOfIndexEntry(), ios::end );
      indexIn.read( (char *) &mLastId1Start, sizeof( mLastId1Start ) );
      indexIn.read( (char *) &mNumAlignsWithLastId1, sizeof( mNumAlignsWithLastId1 ) );
    }

    indexIn.close();
  }
  else
    mLastId1 = -1;
}


void
VecAlignmentPlusWriter::LoadNumAlignsFromData()
{
  if ( IsRegularFile( mStrVectorFilename.c_str() ) )
  {
    ifstream dataIn( mStrVectorFilename.c_str(), ios::in | ios::binary );
    
    if ( ! dataIn.good() ||
         ! mHeaderIO.VerifyDataHeader( &dataIn ) )
    {
      cout << "Unrecognized alignment vector file format in: " << endl;
      cout << mStrVectorFilename << endl;
      cout << "Unable to open alignment vector file " 
           << mStrVectorFilename << " for appending." << endl;
      TracebackThisProcess();
    }

    mHeaderIO.ReadNumAlignments( &dataIn, mNumAligns );

    dataIn.close();
  }
  else
    mNumAligns = 0;
}
      

inline
void
VecAlignmentPlusReader::OpenDataFileForReading()
{
  if ( ! mpDataIn )
  {
    mpDataIn = new ifstream( mStrVectorFilename.c_str(), ios::in | ios::binary );
    
    if ( ! mpDataIn->good() )
    {
      cout << "Unable to open alignment vector file " 
           << mStrVectorFilename << " for reading." << endl;
      TracebackThisProcess();
    }

    if ( ! mHeaderIO.VerifyDataHeader( mpDataIn ) )
    {
      cout << "Unrecognized alignment vector file format in: " << endl;
      cout << mStrVectorFilename << endl;
      TracebackThisProcess();
    }

    // Attach the newly created stream to the read buffer.
    mReadBuffer.Attach( mpDataIn );
  }
}


inline
void
VecAlignmentPlusReader::OpenIndexFileForReading()
{
  if ( ! mpIndexIn )
  {
    if ( mStrIndexFilename.empty() )
    {
      cout << "VecAlignmentPlusReader was initialized without a filename." << endl;
      cout << "Alignment index file was not specified." << endl;
      TracebackThisProcess();
    }

    mpIndexIn = new ifstream( mStrIndexFilename.c_str(), ios::in | ios::binary );
    
    if ( ! mpIndexIn->good() )
    {
      cout << "Unable to open alignment index file " 
           << mStrIndexFilename << " for reading." << endl;
      TracebackThisProcess();
    }
    
    if ( ! mHeaderIO.VerifyIndexHeader( mpIndexIn ) )
    {
      cout << "File " << mStrIndexFilename << " may be in incorrect format." << endl;
      TracebackThisProcess();
    }

    mpIndexIn->seekg( 0, ios::end );
    mIndexSize = mpIndexIn->tellg();
    mpIndexIn->seekg( 0, ios::beg );
  }
}


void
VecAlignmentPlusHeaderIO::WriteDataHeader( ostream *pDataOut )
{
  pDataOut->seekp( 0, ios::beg );
  *pDataOut << mStrVectorHeader << flush;
  // Seek past space for numAligns.
  pDataOut->seekp( sizeof(int), ios::cur );
}

void
VecAlignmentPlusHeaderIO::WriteIndexHeader( ostream *pIndexOut )
{
  pIndexOut->seekp( 0, ios::beg );
  *pIndexOut << mStrIndexHeader << flush;
}

bool
VecAlignmentPlusHeaderIO::VerifyDataHeader( istream *pDataIn )
{
  pDataIn->seekg( 0, ios::beg );

  char *pHeaderData = new char[ mStrVectorHeader.size() + 1 ];
  pDataIn->read( pHeaderData, mStrVectorHeader.size() );

  pHeaderData[ mStrVectorHeader.size() ] = 0;

  bool verified = ( 0 == strcmp( pHeaderData, mStrVectorHeader.c_str() ) );

  delete [] pHeaderData;

  return verified;
}


bool
VecAlignmentPlusHeaderIO::VerifyIndexHeader( istream *pIndexIn )
{
  pIndexIn->seekg( 0, ios::beg );

  char *pHeaderData = new char[ mStrIndexHeader.size() + 1 ];
  pIndexIn->read( pHeaderData, mStrIndexHeader.size() );

  pHeaderData[ mStrIndexHeader.size() ] = 0;

  bool verified = ( 0 == strcmp( pHeaderData, mStrIndexHeader.c_str() ) );

  delete [] pHeaderData;
  
  return verified;
}


inline
void
VecAlignmentPlusHeaderIO::WriteNumAlignments( const String &strDataFilename, const int numAligns )
{
  ofstream dataOut( strDataFilename.c_str(), ios::in | ios::out | ios::binary );
  
  if ( ! dataOut.good() )
  {
    cout << "Unable to open alignment data file " 
         << strDataFilename << " for writing." << endl;
    TracebackThisProcess();
  }
  
  dataOut.seekp( mStrVectorHeader.size(), ios::beg );
  dataOut.write( (char *) &numAligns, sizeof( numAligns ) );
}

inline
void
VecAlignmentPlusHeaderIO::WriteNumAlignments( ostream *pDataOut, const int numAligns )
{
  pDataOut->seekp( mStrVectorHeader.size(), ios::beg );
  pDataOut->write( (char *) &numAligns, sizeof( numAligns ) );
}

inline
void
VecAlignmentPlusHeaderIO::ReadNumAlignments( istream *pDataIn, int &numAligns )
{
  pDataIn->seekg( mStrVectorHeader.size(), ios::beg );
  pDataIn->read( (char *) &numAligns, sizeof( numAligns ) );
}

inline
int
VecAlignmentPlusHeaderIO::GetDataHeaderSize()
{
  return mStrVectorHeader.size() + sizeof( int );
}

inline
int
VecAlignmentPlusHeaderIO::GetIndexHeaderSize()
{
  return mStrIndexHeader.size();
}


int
VecAlignmentPlusReader::GetSize()
{
  if ( mNumAligns < 0 )
  {
    OpenDataFileForReading();
    
    mHeaderIO.ReadNumAlignments( mpDataIn, mNumAligns );
  }

  return mNumAligns;
}

inline
void
VecAlignmentPlusWriter::UnpackAlignment( const alignment_plus &anAlign )
{
  id1 = anAlign.Id1();
  id2 = anAlign.Id2();
  rc2 = anAlign.Rc2();
  score = anAlign.score;

  // We actually call a version of packalign::Unpack() instead of
  // alignment::Unpack() on anAlign.a so we can specify a non-negative
  // numBlocks, which prevents unnecessary resizes of the gaps and
  // lens avectors.  Note that the length member of gaps and lens
  // should therefore be ignored and numBlocks used instead.
  const packalign *packalignPtr = &anAlign.a;
  packalignPtr->Unpack( pos1, pos2, gaps, lens, numBlocks );
  errors = anAlign.a.Errors();
}

inline
void
AlignmentPlusReader::PackAlignment( alignment_plus &anAlign )
{
  anAlign.SetId1( id1 );
  anAlign.SetId2( id2 );
  anAlign.SetRc2( rc2 );
  anAlign.score = score;

  anAlign.a.Set( pos1, pos2, errors, gaps, lens, numBlocks );
}

inline
void
VecAlignmentPlusWriter::WriteIndexEntry()
{
  // If mLastId1 is valid, write out its entry.  An index entry is
  // composed of the position in the data file of the start of the
  // first alignment with id1 == mLastId1 (8 bytes) and the number of
  // alignments where id1 == mLastId1 (4 bytes).
  if ( mNumAlignsWithLastId1 > 0 )
  {
    mpIndexOut->write( (char *) &mLastId1Start, sizeof( mLastId1Start ) );
    mpIndexOut->write( (char *) &mNumAlignsWithLastId1, sizeof( mNumAlignsWithLastId1 ) );
  }

  // Reset the values.  The next start is the offset in the file after
  // the last buffer flush plus the offset in the buffer after the last
  // value copy.
  mLastId1Start = mpDataOut->tellp() + streamoff( mpBufferCursor - mpBufferStart );
  mNumAlignsWithLastId1 = 0;

  // Write out trivial values for all the id1s between the last one and this one.
  int missingId1 = mLastId1;
  for ( ++missingId1; missingId1 < id1; ++missingId1 )
  {
    mpIndexOut->write( (char *) &mLastId1Start, sizeof( mLastId1Start ) );
    mpIndexOut->write( (char *) &mNumAlignsWithLastId1, sizeof( mNumAlignsWithLastId1 ) );
  }
}

inline
void
VecAlignmentPlusReader::ReadIndexEntry( const int sequenceId )
{
  streampos indexPos = 
    mHeaderIO.GetIndexHeaderSize() + sequenceId * SizeOfIndexEntry();

  if ( indexPos >= mIndexSize )
  {
    mRequestedId1Start = mHeaderIO.GetDataHeaderSize();
    mNumAlignsWithRequestedId1 = 0;
  }
  else
  {
    mpIndexIn->seekg( indexPos, ios::beg );

    mpIndexIn->read( (char *) &mRequestedId1Start, sizeof( mRequestedId1Start ) );
    mpIndexIn->read( (char *) &mNumAlignsWithRequestedId1, sizeof( mNumAlignsWithRequestedId1 ) );
  }
}


inline
bool
VecAlignmentPlusWriter::BlockIsPackable( const int gap, const int len )
{
  // For a block to be packable, there can be no bits in the length
  // set outside the storable region, and the bits in the gap outside
  // the magnitude region must either all be 1 or all be 0.
  mReverseMaskedGap = gap & BLOCK_gap_sign_bits;
  return ! ( ( len & BLOCK_len_unpackable_bits ) ||
             ( mReverseMaskedGap       ) &&
             ( mReverseMaskedGap       ) != BLOCK_gap_sign_bits );
}


inline
void
VecAlignmentPlusWriter::WriteAlignmentHeader()
{
  mHeaderByte = HB_all_flags_set;

  if ( id1 != mLastId1 )
    mHeaderByte &= ~HB_use_previous_id1;
  
  if ( ! rc2 )
    mHeaderByte &= ~HB_rc2;

  // If any of the bits are set in pos1, pos2, or errors outside of
  // their respectable storable regions, we can't pack them together.
  if ( ( pos1   & PPE_pos1_unpackable_bits   ) ||
       ( pos2   & PPE_pos2_unpackable_bits   ) ||
       ( errors & PPE_errors_unpackable_bits ) )
    mHeaderByte &= ~HB_ppe_is_packed;

  int blockLimit = min( numBlocks, HB_blocks_in_header );
  for ( int blockIdx = 0; blockIdx < blockLimit; ++blockIdx )
    if ( ! BlockIsPackable( gaps( blockIdx ), lens( blockIdx ) ) )
      mHeaderByte &= ~HB_block_is_packed[ blockIdx ];

  CopyValueToBuffer( mHeaderByte );
}

inline
void
AlignmentPlusReader::ReadAlignmentHeader()
{
  mpReadBuffer->CopyValueFromBuffer( mHeaderByte );

  rc2 = mHeaderByte & HB_rc2;
}


inline
void
VecAlignmentPlusWriter::WriteIds()
{
  if ( ! ( mHeaderByte & HB_use_previous_id1 ) )
    CopyValueToBuffer( id1 );

  CopyValueToBuffer( id2 );
}

inline
void
AlignmentPlusReader::ReadIds()
{
  if ( ! ( mHeaderByte & HB_use_previous_id1 ) )
  {
    mpReadBuffer->CopyValueFromBuffer( id1 );
    mLastId1 = id1;
  }
  else
    id1 = mLastId1;

  mpReadBuffer->CopyValueFromBuffer( id2 );
}


inline
void
VecAlignmentPlusWriter::WritePos1Pos2Errors()
{
  if ( mHeaderByte & HB_ppe_is_packed )
  {
    mPackedPPE = 
      ( static_cast<unsigned int>( pos1 )   << PPE_pos1_shift ) | 
      ( static_cast<unsigned int>( pos2 )   << PPE_pos2_shift ) |
      ( static_cast<unsigned int>( errors ) << PPE_errors_shift );
    
    CopyValueToBuffer( mPackedPPE );
  }
  else
  {
    CopyValueToBuffer( pos1 );
    CopyValueToBuffer( pos2 );
    CopyValueToBuffer( errors );
  }
}
    
inline
void
AlignmentPlusReader::ReadPos1Pos2Errors()
{
  if ( mHeaderByte & HB_ppe_is_packed )
  {
    mpReadBuffer->CopyValueFromBuffer( mPackedPPE );

    pos1   = ( mPackedPPE >> PPE_pos1_shift   ) & PPE_pos1_packable_bits_short;
    pos2   = ( mPackedPPE >> PPE_pos2_shift   ) & PPE_pos2_packable_bits_short;
    errors = ( mPackedPPE >> PPE_errors_shift ) & PPE_errors_packable_bits_short;
  }
  else
  {
    mpReadBuffer->CopyValueFromBuffer( pos1 );
    mpReadBuffer->CopyValueFromBuffer( pos2 );
    mpReadBuffer->CopyValueFromBuffer( errors );
  }
}
    

inline
void
VecAlignmentPlusWriter::WriteScore()
{
  CopyValueToBuffer( score );
}

inline
void
AlignmentPlusReader::ReadScore()
{
  mpReadBuffer->CopyValueFromBuffer( score );
}


inline
void
VecAlignmentPlusWriter::WriteNumBlocks()
{
  mShortNumBlocks = static_cast<unsigned short>( numBlocks );
  CopyValueToBuffer( mShortNumBlocks );
}

inline
void
AlignmentPlusReader::ReadNumBlocks()
{
  mpReadBuffer->CopyValueFromBuffer( mShortNumBlocks );

  numBlocks = mShortNumBlocks;
}


inline
void
VecAlignmentPlusWriter::WriteBuffer( bool cached )
{
  if ( ! cached || 
       mpBufferEnd - mpBufferCursor < maxOutputChunkSize )
  {
    mpDataOut->write( (char *) mpBufferStart, mpBufferCursor - mpBufferStart );
    mpBufferCursor = mpBufferStart;
  }
}

inline
void
VecAlignmentPlusWriter::PackBlock( const int gap, const int len )
{
  mPackedBlock = ( static_cast<unsigned short>( gap ) << BLOCK_gap_shift |
                   static_cast<unsigned short>( len ) << BLOCK_len_shift );
}

inline
void
AlignmentPlusReader::UnpackBlock( int &gap, int &len )
{
  gap = mPackedBlock >> BLOCK_gap_shift;

  // If the sign bit is set on the gap, we need to restore the rest of
  // the sign bits.
  if ( gap & BLOCK_gap_sign_bit )
    gap |= BLOCK_gap_sign_bits;

  len = mPackedBlock & BLOCK_len_packable_bits_short;
}


inline
void
VecAlignmentPlusWriter::WriteBlocks()
{
  int blockLimit = min( numBlocks, HB_blocks_in_header );

  int blockIdx;
  for ( blockIdx = 0; blockIdx < blockLimit; ++blockIdx )
  {
    if ( mHeaderByte & HB_block_is_packed[ blockIdx ] )
    {
      PackBlock( gaps(blockIdx), lens(blockIdx) );
      CopyValueToBuffer( mPackedBlock );
    }
    else
    {
      CopyValueToBuffer( gaps( blockIdx ) );
      CopyValueToBuffer( lens( blockIdx ) );
    }
  }

  WriteBuffer();
  
  while ( blockIdx < numBlocks )
  {
    mpBlockCompressionByte = mpBufferCursor;
    ++mpBufferCursor;

    *mpBlockCompressionByte = BCB_all_blocks_packed; // assume all are compressed
    
    blockLimit = min( numBlocks, blockIdx + 8 );

    for ( int bcbIdx = 0; blockIdx < blockLimit; ++blockIdx, ++bcbIdx )
    {
      mGap = gaps( blockIdx );
      mLen = lens( blockIdx );
      mReverseMaskedGap = mGap & BLOCK_gap_sign_bits;

      if ( ! BlockIsPackable( mGap, mLen ) )
      {
        // turn appropriate bit off
        *mpBlockCompressionByte &= ~( BCB_block_is_packed[ bcbIdx ] );
        
        CopyValueToBuffer( mGap );
        CopyValueToBuffer( mLen );
      }
      else
      {
        PackBlock( gaps(blockIdx), lens(blockIdx) );
        CopyValueToBuffer( mPackedBlock );
      }
    }

    WriteBuffer();
  }
}


inline
void
AlignmentPlusReader::ReadBlocks()
{
  mpReadBuffer->CheckBuffer( numBlocks * BLOCK_max_size );

  int blockLimit = min( numBlocks, HB_blocks_in_header );

  int blockIdx = 0;
  
  // Unpack blocks whose compression bits are stored in the header block.

  // special case for all blocks packed
  if ( ( mHeaderByte & HB_all_blocks_packed ) == HB_all_blocks_packed )
  {
    for ( ; blockIdx < blockLimit; ++blockIdx )
    {
      mpReadBuffer->CopyValueFromBuffer( mPackedBlock );
      UnpackBlock( gaps( blockIdx ), lens( blockIdx ) );
    }
  }

  // general case
  else
  {
    for ( ; blockIdx < blockLimit; ++blockIdx )
    {
      if ( mHeaderByte & HB_block_is_packed[ blockIdx ] )
      {
        mpReadBuffer->CopyValueFromBuffer( mPackedBlock );
        UnpackBlock( gaps( blockIdx ), lens( blockIdx ) );
      }
      else
      {
        mpReadBuffer->CopyValueFromBuffer( gaps( blockIdx ) );
        mpReadBuffer->CopyValueFromBuffer( lens( blockIdx ) );
      }
    }
  }

  // Unpack remaining blocks, if any.

  while ( blockIdx < numBlocks )
  {
    blockLimit = min( numBlocks, blockIdx + 8 );

    // special case for all 8 blocks packed
    if ( mpReadBuffer->AtCursor() == BCB_all_blocks_packed )
    {
      mpReadBuffer->Skip(1);

      for ( ; blockIdx < blockLimit; ++blockIdx )
      {
        mpReadBuffer->CopyValueFromBuffer( mPackedBlock );
        UnpackBlock( gaps( blockIdx ), lens( blockIdx ) );
      }
    }

    // general case
    else
    {
      mpReadBuffer->CopyValueFromBuffer( mBlockCompressionByteCopy );

      for ( int bcbIdx = 0; blockIdx < blockLimit; ++blockIdx, ++bcbIdx )
      {
        if ( mBlockCompressionByteCopy & BCB_block_is_packed[ bcbIdx ] )
        {
          mpReadBuffer->CopyValueFromBuffer( mPackedBlock );
          UnpackBlock( gaps( blockIdx ), lens( blockIdx ) );
        }
        else
        {
          mpReadBuffer->CopyValueFromBuffer( gaps( blockIdx ) );
          mpReadBuffer->CopyValueFromBuffer( lens( blockIdx ) );
        }
      }
    }
  }
}


inline
void
VecAlignmentPlusWriter::WriteUnpackedAlignment()
{
  WriteAlignmentHeader();
  WriteIds();
  WritePos1Pos2Errors();
  WriteScore();
  WriteNumBlocks();

  WriteBuffer();

  WriteBlocks();
}

inline
void
AlignmentPlusReader::ReadUnpackedAlignment()
{
  mpReadBuffer->CheckBuffer( maxHeaderSize );

  ReadAlignmentHeader();
  ReadIds();
  ReadPos1Pos2Errors();
  ReadScore();
  ReadNumBlocks();

  ReadBlocks();
}


inline
bool
VecAlignmentPlusWriter::IndexingIsOn()
{
  return ! mStrIndexFilename.empty();
}


void
VecAlignmentPlusWriter::Write( const vec<alignment_plus> &vecAligns )
{
  OpenDataFileForWriting();
  if ( IndexingIsOn() )
    OpenIndexFileForWriting();

  if ( mStrVectorFilename.empty() )
    mHeaderIO.WriteNumAlignments( mpDataOut, vecAligns.size() );
  else
    mHeaderIO.WriteNumAlignments( mStrVectorFilename, vecAligns.size() );

  id1 = -1;
  mLastId1 = -1;
  mNumAlignsWithLastId1 = 0;
  for ( unsigned int alignIdx = 0; alignIdx < vecAligns.size(); ++alignIdx )
  {
    const alignment_plus &anAlign = vecAligns[ alignIdx ];

    UnpackAlignment( anAlign );

    if ( IndexingIsOn() )
    {
      if ( mLastId1 > id1 )
      {
        cout << "Alignments not sorted.  Unable to index." << endl;
        cout << "Id1 at position " << alignIdx-1 << " of vector to write: " << mLastId1 << endl;
        cout << "Id1 at position " << alignIdx   << " of vector to write: " << id1 << endl;

        TracebackThisProcess();
      }
      
      if ( mLastId1 != id1 )
        WriteIndexEntry();
    }

    WriteUnpackedAlignment();

    mLastId1 = id1;
    ++mNumAlignsWithLastId1;
  }

  if ( IndexingIsOn() )
    WriteIndexEntry();

  WriteBuffer( false );
 
  mpDataOut->flush();
  if ( IndexingIsOn() )
    mpIndexOut->flush();
}


void
VecAlignmentPlusWriter::Write( const vec_alignment_plus &vecAlignPlus )
{
  OpenDataFileForWriting();
  if ( IndexingIsOn() )
    OpenIndexFileForWriting();

  if ( mStrVectorFilename.empty() )
    mHeaderIO.WriteNumAlignments( mpDataOut, vecAlignPlus.GetNumberAlignments() );
  else
    mHeaderIO.WriteNumAlignments( mStrVectorFilename, vecAlignPlus.GetNumberAlignments() );

  alignment_plus anAlign;
  id1 = -1;
  mLastId1 = -1;
  mNumAlignsWithLastId1 = 0;
  for ( unsigned int alignIdx = 0; alignIdx < vecAlignPlus.GetNumberAlignments(); ++alignIdx )
  {
    vecAlignPlus.GetAlignment( anAlign, alignIdx );

    UnpackAlignment( anAlign );

    if ( IndexingIsOn() )
    {
      if ( mLastId1 > id1 )
      {
        cout << "Alignments not sorted.  Unable to index." << endl;
        cout << "Id1 at position " << alignIdx-1 << " of vector to write: " << mLastId1 << endl;
        cout << "Id1 at position " << alignIdx   << " of vector to write: " << id1 << endl;

        TracebackThisProcess();
      }
      
      if ( mLastId1 != id1 )
        WriteIndexEntry();
    }

    WriteUnpackedAlignment();

    mLastId1 = id1;
    ++mNumAlignsWithLastId1;
  }

  if ( IndexingIsOn() )
    WriteIndexEntry();

  WriteBuffer( false );
 
  mpDataOut->flush();
  if ( IndexingIsOn() )
    mpIndexOut->flush();
}

void
VecAlignmentPlusWriter::Append( const vec<alignment_plus> &vecAligns )
{
  LoadNumAlignsFromData();
  LoadLastId1FromIndex();

  OpenDataFileForAppending();
  if ( IndexingIsOn() )
    OpenIndexFileForAppending();

  mHeaderIO.WriteNumAlignments( mStrVectorFilename, mNumAligns + vecAligns.size() );

  mpDataOut->seekp( 0, ios::end );

  if ( IndexingIsOn() )
    if ( mLastId1 < 0 )
      mpIndexOut->seekp( 0, ios::end );
    else
      mpIndexOut->seekp( -SizeOfIndexEntry(), ios::end );

  for ( int alignIdx = 0; alignIdx < vecAligns.isize(); ++alignIdx )
  {
    const alignment_plus &anAlign = vecAligns[ alignIdx ];

    UnpackAlignment( anAlign );

    if ( IndexingIsOn() )
    {
      if ( mLastId1 > id1 )
      {
        cout << "Alignments not sorted.  Unable to index." << endl;

        if ( alignIdx == 0 )
          cout << "Id1 of last alignment in file: " << mLastId1 << endl;
        else
          cout << "Id1 at position " << alignIdx-1 << " of vector to write: " << mLastId1 << endl;

        cout << "Id1 at position " << alignIdx   << " of vector to write: " << id1 << endl;

        TracebackThisProcess();
      }
      
      if ( mLastId1 != id1 )
        WriteIndexEntry();
    }

    WriteUnpackedAlignment();

    mLastId1 = id1;
    ++mNumAlignsWithLastId1;
  }

  if ( IndexingIsOn() )
    WriteIndexEntry();

  WriteBuffer( false );

  mpDataOut->flush();
  if ( IndexingIsOn() )
    mpIndexOut->flush();
}


void
VecAlignmentPlusReader::ReadAll( vec<alignment_plus> &vecAligns )
{
  OpenDataFileForReading();

  int numAligns;
  mHeaderIO.ReadNumAlignments( mpDataIn, numAligns );
  
  vecAligns.resize( numAligns );
  
  mReader.SetLastId1( -1 );
  for ( unsigned int alignIdx = 0; alignIdx < vecAligns.size(); ++alignIdx )
  {
    mReader.ReadUnpackedAlignment();
    mReader.PackAlignment( vecAligns[ alignIdx ] );
  }
}

void
VecAlignmentPlusReader::ReadHalf( vec<alignment_plus> &vecAligns,
                                  vec<Bool> &vecAlignIsFlipped,
                                  vec<int>  &vecAlignIds )
{
  OpenDataFileForReading();
  OpenIndexFileForReading();

  int numAligns;
  mHeaderIO.ReadNumAlignments( mpDataIn, numAligns );

  vecAligns.clear();
  vecAligns.reserve( numAligns/2 );

  vecAlignIsFlipped.clear();
  vecAlignIsFlipped.reserve( numAligns );

  vecAlignIds.clear();
  vecAlignIds.reserve( numAligns );

  // If we cache the alignment ids, the runtime is cut by two thirds
  // at a memory cost of about 20%.

  vec<int> vecAlignId1;
  vecAlignId1.reserve( numAligns/2 );

  vec<int> vecAlignId2;
  vecAlignId2.reserve( numAligns/2 );

  int indexDataSize = mIndexSize - mHeaderIO.GetIndexHeaderSize();
  int biggestId1 = mIndexSize/SizeOfIndexEntry() - 1;

  vec<int> vecFirstAlignWithId( biggestId1 + 1, -1 );


  int id1, id2;

  mReader.SetLastId1( -1 );
  for ( int alignIdx = 0; alignIdx < numAligns; ++alignIdx )
  {
    mReader.ReadUnpackedAlignment();

    id1 = mReader.GetId1();
    id2 = mReader.GetId2();
    if ( id1 < id2 ) 
    {
      vecAlignIds.push_back( vecAligns.size() );
      vecAlignIsFlipped.push_back( False );
      
      vecAlignId1.push_back( id1 );
      vecAlignId2.push_back( id2 );

      vecAligns.push_back( alignment_plus() );
      mReader.PackAlignment( vecAligns.back() );

      if ( vecFirstAlignWithId[ id1 ] < 0 )
        vecFirstAlignWithId[ id1 ] = vecAlignIds.back();
    }
    else
    {
      int existingAlignIdx = vecFirstAlignWithId[ id2 ];
      
      if ( existingAlignIdx == -1 )
      {
        cout << mStrVectorFilename << " is not symmetrical." << endl;
        cout << "Found alignment with id1=" << id1 << " and id2=" << id2 
             << " but no alignment with id1=" << id2 << " and id2=" << id1 << "." << endl;
        TracebackThisProcess();
      }

      for ( ; existingAlignIdx < (int)vecAligns.size(); ++existingAlignIdx )
      {
        if ( vecAlignId2[ existingAlignIdx ] == id1 )
          break;
        if ( vecAlignId1[ existingAlignIdx ] != id2 )
          break;
      }

      if ( vecAlignId1[ existingAlignIdx ] != id2 )
      {
        cout << mStrVectorFilename << " is not symmetrical." << endl;
        cout << "Found alignment with id1=" << id1 << " and id2=" << id2 
             << " but no alignment with id1=" << id2 << " and id2=" << id1 << "." << endl;
        TracebackThisProcess();
      }
        
      vecAlignIsFlipped.push_back( True );
      vecAlignIds.push_back( existingAlignIdx );
    }
  }

  ForceAssertEq( vecAlignIsFlipped.size(), vecAlignIds.size() );
}


int
VecAlignmentPlusReader::GetNumAlignmentsForId( const int sequenceId )
{
  OpenIndexFileForReading();

  ReadIndexEntry( sequenceId );

  return mNumAlignsWithRequestedId1;
}


void
VecAlignmentPlusReader::ReadSubset( vec<alignment_plus> &vecAligns,
                                    const int sequenceId,
                                    bool append )
{
  if ( ! append )
    vecAligns.clear();

  OpenDataFileForReading();
  OpenIndexFileForReading();

  ReadIndexEntry( sequenceId );

  if ( mNumAlignsWithRequestedId1 == 0 )
    return;

  mReadBuffer.JumpTo( mRequestedId1Start );

  unsigned int alignIdx = vecAligns.size();

  vecAligns.resize( vecAligns.size() + mNumAlignsWithRequestedId1 );
  
  mReader.SetLastId1( -1 );
  for ( ; alignIdx < vecAligns.size(); ++alignIdx )
  {
    mReader.ReadUnpackedAlignment();

    if ( mReader.GetId1() != sequenceId )
    {
      cout << "Alignment vector not indexed properly." << endl;
      cout << "Expected id1: " << sequenceId << ", "
           << "found id1: " << mReader.GetId1() << endl;
      TracebackThisProcess();
    }

    mReader.PackAlignment( vecAligns[ alignIdx ] );
  }
}  


void
VecAlignmentPlusReader::ReadSubset( vec<alignment_plus> &vecAligns,
                                    const vec<int> &vecSequenceIds,
                                    bool append )
{
  OpenIndexFileForReading();

  int numAlignsToRead = 0;
  for ( vec<int>::const_iterator sequenceIdIter = vecSequenceIds.begin();
        sequenceIdIter != vecSequenceIds.end(); ++sequenceIdIter )
  {
    ReadIndexEntry( *sequenceIdIter );
    numAlignsToRead += mNumAlignsWithRequestedId1;
  }

  if ( ! append )
    vecAligns.clear();

  vecAligns.reserve( vecAligns.size() + numAlignsToRead );

  for ( vec<int>::const_iterator sequenceIdIter = vecSequenceIds.begin();
        sequenceIdIter != vecSequenceIds.end(); ++sequenceIdIter )
    ReadSubset( vecAligns, *sequenceIdIter, true );
}


void
VecAlignmentPlusReader::ReadNext( alignment_plus &nextAlign )
{
  OpenIndexFileForReading();

  mReader.ReadUnpackedAlignment();
  mReader.PackAlignment( nextAlign );
}

  
void 
OldEfficientRead( istream& in,
                  alignment_plus& ap,
                  int &p,
                  avector<unsigned int>& x,
                  int previous_ap_id1 );

void
VecAlignmentPlusWriter::Convert( const String &strOldVectorFilename,
                                 const int limit )
{
  ForceAssert( IsRegularFile( strOldVectorFilename ) );

  istream *in_ptr = 0;
  procbuf *pb_ptr = 0;
  if ( strOldVectorFilename.Contains( ".gz", -1 ) )
  {
    String cmd = "gzip -dc ";
    cmd += strOldVectorFilename;
    pb_ptr = new procbuf( cmd.c_str(), ios::in | ios::binary );
    in_ptr = new istream( pb_ptr );
  }
  else
  {
    in_ptr = new ifstream( strOldVectorFilename.c_str(), ios::in | ios::binary );
  }
  
  istream &in = *in_ptr;
  
  OpenDataFileForWriting();
  if ( IndexingIsOn() )
    OpenIndexFileForWriting();

  int numAligns;
  in >> numAligns;
  char c;
  in.get(c);
  
  int realLimit = ( limit == 0 ? numAligns : limit );

  // cout << "Converting " << realLimit << " alignments:" << endl;

  avector<unsigned int> x(BufSize);
  in.read( (char*) &(x(0)), x.length * sizeof(int) );
  int p = 0;
  alignment_plus anAlign;

  if ( mStrVectorFilename.empty() )
    mHeaderIO.WriteNumAlignments( mpDataOut, realLimit );
  else
    mHeaderIO.WriteNumAlignments( mStrVectorFilename, realLimit );

  mLastId1 = -1;
  for ( int alignIdx = 0; alignIdx < realLimit; ++alignIdx )
  {
    OldEfficientRead( in, anAlign, p, x, mLastId1 );

    // if ( (alignIdx + 1) % 1000000 == 0 )
    //   Dot( cout, alignIdx / 1000000 );

    UnpackAlignment( anAlign );

    if ( IndexingIsOn() )
    {
      if ( mLastId1 > id1 )
      {
        cout << "Alignments not sorted.  Unable to index." << endl;
        cout << "Id1 at position " << alignIdx-1 << " of vector to write: " << mLastId1 << endl;
        cout << "Id1 at position " << alignIdx   << " of vector to write: " << id1 << endl;

        TracebackThisProcess();
      }
      
      if ( id1 != mLastId1 )
        WriteIndexEntry();
    }

    WriteUnpackedAlignment();

    mLastId1 = id1;
    ++mNumAlignsWithLastId1;
  }

  WriteBuffer( false );

  if ( IndexingIsOn() )
    WriteIndexEntry();

  mpDataOut->flush();
  if ( IndexingIsOn() )
    mpIndexOut->flush();

  // cout << endl;

  if ( in_ptr ) delete in_ptr;
  if ( pb_ptr ) delete pb_ptr;
}
  
  
void
VecAlignmentPlusWriter::TestPackBlock()
{
  cout.setf( ios::hex );
  cout.setf( ios::showbase );
  cout.unsetf( ios::dec );

  unsigned short expected = 0x0000;
  for ( int gap = 0; gap < 16; ++gap, expected += 0x0800 )
  {
    PackBlock( gap, 0 );
    ForceAssertEq( mPackedBlock, expected );
  }

  for ( int gap = -16; gap < 0; ++gap, expected += 0x0800 )
  {
    PackBlock( gap, 0 );
    ForceAssertEq( mPackedBlock, expected );
  }

  cout.setf( ios::dec );
}


void 
OldEfficientRead( istream& in,
                  alignment_plus& ap,
                  int &p,
                  avector<unsigned int>& x,
                  int previous_ap_id1 )
{
  if ( BufSize - p < MaxRecordSize ) {
    for ( int j = 0; j < BufSize - p; j++ )
      x(j) = x( p + j );
    in.read( (char*) &(x( BufSize - p )), p * sizeof(int) );
    p = 0;
  }
  unsigned char control = x(p) >> 24;
  Bool rc2 = x(p) >> 16;
  unsigned short block_length = (unsigned short) x(p++);
  if ( control == 0 ) ap.SetId1( x(p++) );
  else ap.SetId1( previous_ap_id1 );
  ap.SetId2( x(p++) );
  ap.SetRc2(rc2);
  int pos1, pos2, errors;
  if ( control == 0 ) {
    pos1 = x(p++);
    pos2 = x(p++);
    errors = x(p++);
  }
  else {
    pos1 = x(p) % 1024;
    unsigned int y = (x(p) - pos1) / 1024;
    pos2 = y % 1024;
    errors = (y - pos2) / 1024;
    p++;
  }
  
  static avector<int> gaps, lengths;
  if ( gaps.length < block_length )
  {    gaps.Setsize(block_length);
       lengths.Setsize(block_length);    }
  
  // Sante --- Tue Sep 17 17:53:56 EDT 2002
  //  Initialize gaps and lengths in the case block_length=0.
  if ( block_length == 0 ) {
    gaps.Setsize( 0 );
    lengths.Setsize( 0 );
  }
   
  for ( unsigned short j = 0; j < block_length; j++ ) {
    if ( control == 0 ) {
      gaps(j) = x(p++);
      lengths(j) = x(p++);
    }
    else {
      gaps(j) = (x(p) >> 16) - 1024;
      lengths(j) = (unsigned short) x(p);
      p++;
    }
  }
  
  ap.a.Set( pos1, pos2, errors, gaps, lengths, block_length );
  ap.score = *( (float*) &x(p) );
  p += sizeof(float)/sizeof(int);
}

