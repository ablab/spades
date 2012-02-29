///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef INDEXED_ALIGNMENT_PLUS_VECTOR_H
#define INDEXED_ALIGNMENT_PLUS_VECTOR_H

// This collection of classes manages the speedy I/O of indexed
// vectors of alignment_plus objects.


#include "Alignment.h"


// This class encapsulates the header information for both vectors and indices.

class VecAlignmentPlusHeaderIO
{
 public:
  void  WriteDataHeader( ostream *pDataOut );
  void  WriteIndexHeader( ostream *pIndexOut );

  bool  VerifyDataHeader( istream *pDataIn );
  bool  VerifyIndexHeader( istream *pIndexIn );

  void  WriteNumAlignments( const String &strDataFilename, const int numAligns );
  void  WriteNumAlignments( ostream *pDataOut, const int numAligns );
  void  ReadNumAlignments( istream *pDataIn, int &numAligns );

  int   GetDataHeaderSize();
  int   GetIndexHeaderSize();

 private:
  static const String mStrVectorHeader;
  static const String mStrIndexHeader;
};


// This class encapsulates the buffering of input from a file.
  
class ReadBuffer
{
 public:
  ReadBuffer( const int bufferSize );
  ~ReadBuffer();

  void  Attach( istream *pIstream ) { mpIstream = pIstream; }

  template <typename T>
  void  CopyValueFromBuffer( T &value )
  {
    memcpy( &value, mpBufferCursor, sizeof(T) );
    mpBufferCursor += sizeof(T);
  }

  const unsigned char* Cursor()   { return  mpBufferCursor; }
  unsigned char  AtCursor() { return *mpBufferCursor; }

  const char* CursorAsString() { return (const char*) mpBufferCursor; }

  void  CheckBuffer( const int amountToBuffer );

  void  JumpTo( const streampos filePos );

  void  Skip( const unsigned int bytesToSkip ) { mpBufferCursor += bytesToSkip; }

 private:
  istream       *mpIstream;
  int            mBufferSize;
  unsigned char *mpBufferStart;
  unsigned char *mpBufferEnd;
  unsigned char *mpBufferCursor;
  // Track whether the buffer is in sync with the stream (i.e. the
  // curr position of the stream corresponds to the end of the buffer).
  bool           mBufferIsStale;
};


// This class can read in alignments from the given buffer.

class AlignmentPlusReader
{
 public:
  AlignmentPlusReader( ReadBuffer *pReadBuffer );

  void  ReadUnpackedAlignment();
  void  PackAlignment( alignment_plus &anAlign );

  void  ReadAlignmentHeader();
  void  ReadIds();
  void  ReadPos1Pos2Errors();
  void  ReadScore();
  void  ReadNumBlocks();

  void  UnpackBlock( int &gap, int &len );
  void  ReadBlocks();

  int   GetId1() const { return id1; }
  int   GetId2() const { return id2; }

  int   GetLastId1( ) const { return mLastId1; }
  void  SetLastId1( const int lastId1 ) { mLastId1 = lastId1; }
  
 private:
  ReadBuffer *mpReadBuffer;

  int      mLastId1;

  unsigned char mHeaderByte;

  int   id1;
  int   id2;
  Bool  rc2;
  int   pos1, pos2, errors;
  float score;
  int   numBlocks;
  avector<int> gaps, lens;

  unsigned int   mPackedPPE;
  unsigned short mShortNumBlocks;
  unsigned char  mBlockCompressionByteCopy;
  unsigned short mPackedBlock;
};


// This class can read in all or part of a vector of alignment_plus
// objects written to disk.

class VecAlignmentPlusReader
{
 public:

  // A VecAlignmentPlusReader object is instantiated with the filename
  // to be read and the buffer size (the default buffer size seems to
  // work well under most circumstances).  The vector file has the
  // given filename, while the index file (if any) has "_index"
  // appended to the given filename.
  VecAlignmentPlusReader( const String &strVectorFilename,
                          const int bufferSize = 100000 );

  // When instantiated in this fashion, the only legal method to call
  // is ReadVector().  ReadVectorSubset() will fail, as the object
  // will not have access to the index.
  VecAlignmentPlusReader( istream *pVectorIfstream,
                          const int bufferSize = 100000 );

  ~VecAlignmentPlusReader();
  
  int
  GetSize();

  int
  GetNumAlignmentsForId( const int sequenceId );

  void
  ReadAll( vec<alignment_plus> &vecAligns );

  // This is intended for exclusive use by vec_alignment_plus.
  void
  ReadHalf( vec<alignment_plus> &vecAligns,
            vec<Bool> &vecAlignIsFlipped,
            vec<int>  &vecAlignIds );

  // (Replace the contents of/Append to) the given vector with the
  // aligns where id1 is the given sequence id.
  void
  ReadSubset( vec<alignment_plus> &vecAligns,
              const int sequenceId, 
              bool append = false );

  // (Replace the contents of/Append to) the given vector with the
  // aligns where id1 is in the given (not necessarily sorted) vector
  // of sequence ids.
  void
  ReadSubset( vec<alignment_plus> &vecAligns,
              const vec<int> &vecSequenceIds, 
              bool append = false );

  void
  ReadNext( alignment_plus &nextAlign );

 private:
  void  OpenDataFileForReading();
  void  OpenIndexFileForReading();
  
  void  ReadIndexEntry( const int sequenceId );

  String mStrVectorFilename;
  String mStrIndexFilename;

  istream *mpDataIn, *mpIndexIn;

  bool mDataPtrIsShared;

  VecAlignmentPlusHeaderIO mHeaderIO;

  int      mNumAligns;
  longlong mIndexSize;

  longlong mRequestedId1Start;
  int      mNumAlignsWithRequestedId1;

  ReadBuffer          mReadBuffer;
  AlignmentPlusReader mReader;
};


// Forward declaration of class used in method signature of
// VecAlignmentPlusWriter::Write().
class vec_alignment_plus;


// This class can write out vectors of alignment_plus objects, either
// indexed or not, saving them either to new files or appending to
// existing ones.

class VecAlignmentPlusWriter
{
 public:

  // A VecAlignmentPlusWriter object is instantiated with the filename
  // to be written, a bool indicating whether a index should be
  // written, and the buffer size (the default buffer size seems to
  // work well under most circumstances).  The vector file has the
  // given filename, while the index file has "_index" appended to the
  // given filename.
  VecAlignmentPlusWriter( const String &strVectorFilename,
                          const bool createIndex = true,
                          const int bufferSize = 100000 );

  // When instantiated in this fashion, no index file will be created,
  // as we don't know what to call it.  Because of this, there is no
  // requirement that the vectors to be saved are sorted.
  VecAlignmentPlusWriter( ostream *pVectorOstream,
                          const int bufferSize = 100000 );

  ~VecAlignmentPlusWriter();

  void
  Write( const vec<alignment_plus> &vecAligns );

  void
  Write( const vec_alignment_plus &vecAlignPlus );

  void
  Append( const vec<alignment_plus> &vecAligns );

  // Convert an old-style vector to the new format.
  void
  Convert( const String &strOldVectorFilename,
           const int limit = 0 );

  // Code to test the correctness of block-packing.
  void
  TestPackBlock();

 private:
  void  OpenDataFileForWriting();
  void  OpenIndexFileForWriting();
  
  void  OpenDataFileForAppending();
  void  OpenIndexFileForAppending();
  
  void  LoadLastId1FromIndex();
  void  LoadNumAlignsFromData();

  bool  IndexingIsOn();
  void  WriteIndexEntry();

  void  UnpackAlignment( const alignment_plus &anAlign );
  void  WriteUnpackedAlignment();

  void  WriteAlignmentHeader();
  void  WriteIds();
  void  WritePos1Pos2Errors();
  void  WriteScore();
  void  WriteNumBlocks();

  bool  BlockIsPackable( const int gap, const int len );
  void  PackBlock  ( const int gap, const int len );
  void  WriteBlocks();

  template <typename T>
  void  CopyValueToBuffer( T &value )
  {
    memcpy( mpBufferCursor, &value, sizeof(T) );
    mpBufferCursor += sizeof(T);
  }

  void  WriteBuffer( bool cached = true );


  String mStrVectorFilename;
  String mStrIndexFilename;

  ostream *mpDataOut, *mpIndexOut;

  bool mDataPtrIsShared;

  VecAlignmentPlusHeaderIO mHeaderIO;

  unsigned char mHeaderByte;

  int      mNumAligns;
  int      mLastId1;
  longlong mLastId1Start;
  int      mNumAlignsWithLastId1;

  int   id1;
  int   id2;
  Bool  rc2;
  int   pos1, pos2, errors;
  float score;
  int   numBlocks;
  avector<int> gaps, lens;

  unsigned int   mPackedPPE;
  unsigned short mShortNumBlocks;
  int            mGap, mLen;
  int            mReverseMaskedGap;
  unsigned char *mpBlockCompressionByte;
  unsigned short mPackedBlock;

  int            mBufferSize;
  unsigned char *mpBufferStart;
  unsigned char *mpBufferEnd;
  unsigned char *mpBufferCursor;
};

#endif
