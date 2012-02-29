///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file FeudalFileWriter.h
 * \author tsharpe
 * \date Mar 13, 2009
 *
 * \brief Utility for writing Feudal Files incrementally.
 * If you can iterate over the elements that go into the feudal file, you can
 * use this class to write a feudal file in a single pass without having to have
 * the entire structure in memory, and without writing three files that then
 * have to be merged.
 */
#ifndef FEUDAL_FEUDALFILEWRITER_H_
#define FEUDAL_FEUDALFILEWRITER_H_

#include "feudal/FeudalControlBlock.h"
#include "feudal/BinaryStream.h"
#include <cstddef>
#include <string>
#include <vector>

/**
 * \class FeudalFileWriter
 * \brief Lets you write a Feudal File incrementally.
 *
 * This class holds a file descriptor open, so be careful how many
 * FeudalFileWriters you instantiate. This class streams the variable-length
 * data directly to the file, and keeps an in-memory vector of just the
 * fixed-length data and the starting offset.  In the destructor, it dumps the
 * offsets and fixed-length data, and then goes back to the start and writes the
 * control block.
 */
class FeudalFileWriter
{
public:
    typedef unsigned int size_type;
    typedef unsigned long offset_type;

    FeudalFileWriter( char const* filename,
                      size_type vecSize,
                      size_type eltSize,
                      size_type fixedLenDataLen,
                      unsigned long estimatedNElements = 100000 );

    ~FeudalFileWriter();

    std::string const& getFilename() const
    { return mWriter.getFilename(); }

    /// number of elements written so far
    unsigned long getNElements() const { return mOffsets.size(); }

    /// Get the writer for variable-length data.
    BinaryWriter& getWriter() { return mWriter; }

    /// Tell us that you've added variable-length data for an element.
    void addElement( size_t varLen, void const* fixedLenData );

    /// Flush everything to make a finished feudal file, but leave it open for
    /// further writing.  It's necessary for just a couple of pieces of code
    /// that need to re-read a feudal file as it's being created.  It's slow.
    void checkPoint();

    void close();

private:
    FeudalFileWriter( FeudalFileWriter const& ); // unimplemented -- no copying
    FeudalFileWriter& operator=( FeudalFileWriter const& ); // unimplemented -- no copying

    void finish();

    BinaryWriter mWriter;
    std::vector<offset_type> mOffsets;
    std::vector<char> mFixedLenData;
    offset_type mNextOffset;
    size_type mVecSize;
    size_type mEltSize;
    size_type mFixedLenDataLen;

    static FeudalControlBlock gInvalidFCB;
};

#endif /* FEUDAL_FEUDALFILEWRITER_H_ */
