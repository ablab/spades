///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file FeudalControlBlock.h
 * \author tsharpe
 * \date Dec 19, 2008
 *
 * \brief The header on Feudal Files.
 */
#ifndef FEUDAL_FEUDALCONTROLBLOCK_H_
#define FEUDAL_FEUDALCONTROLBLOCK_H_

#include "feudal/BinaryStreamTraits.h"
#include "system/Assert.h"
#include "system/file/FileReader.h"
#include <cstddef>

/// Describes the header information for a feudal file.
/// A feudal file pickles a vector of Xs, where X is some vector-like type
/// that contains an indefinite number of elements of type A.
class FeudalControlBlock
{
public:
    /// read one from a file that you won't otherwise be using
    /// final arg is set to the file length, if supplied
    FeudalControlBlock( char const* filename, bool validate=true,
                        size_t* pFileLen=0 )
    { init(FileReader(filename),validate,pFileLen); }


    /// read one from a file that you're going to use subsequently.
    /// final arg is set to the file length, if supplied
    FeudalControlBlock( FileReader const& fr, bool validate=true,
                        size_t* pFileLen=0 )
    { init(fr,validate,pFileLen); }

    /// make a brand new one from scratch.
    FeudalControlBlock( size_t nElements, size_t variableDataLen,
                        unsigned char sizeofFixed, unsigned char sizeofX,
                        unsigned char sizeofA, unsigned int nFiles=1,
                        bool isCompressed=false )
    : mN(nElements), mBitFlags(nFiles&3), mSizeofFixed(sizeofFixed),
      mSizeofX(sizeofX), mSizeofA(sizeofA),
      mVarOffset(variableDataLen+sizeof(*this)),
      mFixedOffset(mVarOffset+(mN+1)*sizeof(size_t))
    { AssertLt(nFiles,4u);
      if ( isCompressed ) mBitFlags |= 4; }

    // compiler-supplied copying and destructor are OK

    /// number of Xs
    size_t getNElements() const
    { return getNFiles()==3 ?
                    mN :
                    ((mFixedOffset - mVarOffset)/sizeof(size_t) - 1UL); }

    /// format of feudal file:  1 file or 3
    unsigned int getNFiles() const
    { return mBitFlags & 3; }

    /// marker for compressed mastervec
    bool isCompressed() const
    { return mBitFlags & 4; }

    /// set marker for compressed mastervec
    void setCompressed( bool isCompressed )
    { if ( isCompressed ) mBitFlags |= 4; else mBitFlags &= ~4; }

    /// not in use, but here to show you where the bit used to be
    bool isBigEndian() const
    { return mBitFlags & 8; }

    /// not in use, but here to show you where the bits used to be
    unsigned int getVersion() const
    { return mBitFlags >> 4; }

    /*
     * first chunk of a feudal file after the header:  the variable-length data
     */

    /// offset of the variable-length data
    size_t getVarDataOffset() const
    { return sizeof(*this); }

    /// total number of bytes of variable-length data
    size_t getVarDataLen() const
    { return mVarOffset - sizeof(*this); }

    /*
     * second chunk of a feudal file after the header:  the table of offsets
     * into the variable length data for each X (plus 1 more offset to delimit
     * the end of the variable-length data).
     */

    /// offset of the table of variable-length data offsets
    size_t getVarTabOffset() const
    { return mVarOffset; }

    /// get the offset of the specified element's variable-length data offset
    size_t getVarTabOffset( size_t idx )
    { return mVarOffset + idx*sizeof(size_t); }

    /// length of the offsets table
    size_t getVarTabLen() const
    { return mFixedOffset - mVarOffset; }

    /*
     * third chunk of a feudal file after the header:  the fixed-length data
     */

    /// offset in file to the fixed-length data
    size_t getFixedOffset() const
    { return mFixedOffset; }

    // we don't actually have any way of calculating the length of the
    // fixed-length data.  if it's not a compressed feudal file, it'll
    // run to the end of the file.  if you know the number of bytes of
    // fixed-len data per X (getSizeofFixed may return 0 or a truncated
    // value -- it's not reliable), you can multiply by the number of X's.

    /*
     * some methods to get the per-element size of stuff, just for
     * sanity-checking that you've hooked up with the right file.
     */
    /// number of bytes of fixed-len data per X.
    /// may be 0, may be truncated to 8 bits!
    unsigned char getSizeofFixed() const
    { return mSizeofFixed; }

    /// sizeof(X) -- may be 0 or truncated to 8 bits!
    unsigned char getSizeofX() const
    { return mSizeofX; }

    /// sizeof(A) -- may be 0 or truncated to 8 bits!
    unsigned char getSizeofA() const
    { return mSizeofA; }

    /// does all the internal consistency checking we can think of
    bool isValid( char const* fileName, size_t fileLen,
                    bool verbose=false ) const;

    /// sets the fixed and variable data offsets to 0, sets number of files to 3
    void to3FileFormat()
    { mFixedOffset = mVarOffset = 0; mBitFlags |= 3; }

    /// opens the file, reads and validates the control block
    /// if verbose is true, reports anomalies to cout
    static bool isGoodFeudalFile( char const* fileName, bool verbose=false );

private:
    void init( FileReader const& fr, bool validate, size_t* pFileLen );

    unsigned int mN; // the number of elements, modulo the size of a uint
    unsigned char mBitFlags; // bit-packed datum: nFiles, compressed flag, big-endian flag and version number
    unsigned char mSizeofFixed; // the per-element size of the fixed length data, modulo the size of a uchar
    unsigned char mSizeofX; // the size of an element, modulo the size of a uchar
    unsigned char mSizeofA; // the size of the inner-vector value_type, modulo the size of a uchar
    size_t mVarOffset; // location of the variable-length data table of offsets
    size_t mFixedOffset; // location of the fixed-length data pool
};

TRIVIALLY_SERIALIZABLE(FeudalControlBlock);

#endif /* FEUDAL_FEUDALCONTROLBLOCK_H_ */
