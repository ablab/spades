///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FEUDAL_FEUDALFILEREADER_H_
#define FEUDAL_FEUDALFILEREADER_H_

#include "feudal/BinaryStream.h"
#include "feudal/FeudalControlBlock.h"
#include "system/Assert.h"
#include "system/SpinLockedData.h"
#include <cstddef>
#include <string>
#include <fstream>

/**
 * \file FeudalFileReader.h
 * \author tsharpe
 * \date Wednesday, December 03, 2008
 *
 * \brief Utility for reading Feudal Files.
 */

/**
 * \class FeudalFileReader
 * \brief Gives you access to the data from a Feudal File.
 */
class FeudalFileReader
{
public:
    FeudalFileReader( char const* filename );
    FeudalFileReader( FeudalFileReader const& that );
    ~FeudalFileReader();

    std::string const& getFilename() const { return mpMapper->getFilename(); }

    FeudalControlBlock const& getFCB() const { return mpMapper->getFCB(); }

    /// Returns the number of inner-vector elements in the file.
    size_t getNElements() const { return mpMapper->getNElements(); }

    size_t getMappedLen() const { return mpMapper->getMappedLen(); }

    /// Get the length of all the variable-length data (in bytes)
    size_t getDataLenTotal() const
    { return mpMapper->getOffset(getNElements())-mpMapper->getOffset(0); }

    /// Get the length of the variable-length data (in bytes)
    /// for the specified element
    size_t getDataLen( size_t ele ) const
    { AssertLt(ele,getNElements());
      return mpMapper->getOffset(ele+1)-mpMapper->getOffset(ele); }

    /// Get a BinaryReader positioned so that it's ready to read the variable-
    /// length data for the specified element
    BinaryReader& getData( size_t ele )
    { AssertLt(ele,getNElements());
      if ( ele != mCurEle++ )
      { mReader.seek(mpMapper->getOffset(ele)); mCurEle = ele+1; }
      return mReader; }

    /// Get a pointer to the fixed-length data for the specified element
    void* getFixedData( size_t ele, size_t bytesPerEle ) const
    { return mpMapper->getFixedData(ele,bytesPerEle); }

    /// Compute a weird number using a method that doesn't belong in this class
    size_t getPreallocation( size_t maxSize, size_t extSize,
                             size_t start, size_t end ) const;

private:
    FeudalFileReader& operator=( FeudalFileReader const& that ); // unimplemented -- no copying

    class Mapper : public SpinLockedData
    {
    public:
        Mapper( char const* filename );
        ~Mapper();

        size_t ref() { return ++mRefCount; }
        size_t deref() { return --mRefCount; }

        std::string const& getFilename() const { return mFilename; }
        FeudalControlBlock const& getFCB() const { return mFCB; }
        size_t getNElements() const { return mFCB.getNElements(); }
        size_t getMappedLen() const { return mMappedLen; }
        void* getFixedData( size_t ele, size_t bytesPerEle ) const
        { AssertLt(ele,getNElements());
          char* result = mFixedData+bytesPerEle*ele;
          Assert(result+bytesPerEle <= mMappedBit+mMappedLen);
          return result; }
        size_t getOffset( size_t ele ) const
        { size_t result;
          memcpy(&result,mOffsetsTable+sizeof(result)*ele,sizeof(result));
          return result; }

    private:
        Mapper( Mapper const& ); // unimplemented -- no copying
        Mapper& operator=( Mapper const& ); // unimplemented -- no copying

        std::string mFilename;
        size_t mRefCount;
        FeudalControlBlock mFCB;
        char* mMappedBit;
        size_t mMappedLen;
        char* mOffsetsTable;
        char* mFixedData;
    };

    Mapper* mpMapper;
    size_t mCurEle;
    BinaryReader mReader;
};

#endif /* FEUDAL_FEUDALFILEREADER_H_ */
