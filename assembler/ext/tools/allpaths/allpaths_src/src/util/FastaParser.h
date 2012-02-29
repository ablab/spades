/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file FastaParser.h
 * \author tsharpe
 * \date Apr 7, 2009
 *
 * \brief Several utility classes for reading FASTA files.
 */
#ifndef UTIL_FASTAPARSER_H_
#define UTIL_FASTAPARSER_H_

#include "feudal/IncrementalWriter.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "Qualvector.h"
#include "String.h"
#include <fstream>
#include <vector>

/// Base class for the parsers for sequence, quality, and bits which follow.
class FastaParser
{
public:
    FastaParser( char const* fileName )
    : mLocalIS(fileName), mIS(mLocalIS), mLineNumber(0UL)
    { init(); }

    FastaParser( std::istream& is )
    : mIS(is), mLineNumber(0UL)
    { init(); }

    ~FastaParser()
    { if ( mLocalIS.is_open() ) mLocalIS.close(); }

protected:
    ulong getLineNumber()
    { return mLineNumber; }

    char* getBuf()
    { return mBuf; }

    char* readLine();

private:
    FastaParser( FastaParser const& ); // unimplemented:  no copying
    FastaParser& operator=( FastaParser const& ); // unimplemented:  no copying

    void init();

    static int const BUF_SIZE = 8192;

    ifstream mLocalIS;
    istream& mIS;
    char mBuf[BUF_SIZE];
    ulong mLineNumber;
};

/// Reads sequence from a FASTA file.
class SeqParser : public FastaParser
{
public:
    SeqParser( char const* fileName )
    :FastaParser(fileName)
    {}

    SeqParser( std::istream& is )
    :FastaParser(is)
    {}

    /// Ref args are set to the data from the next comma-delimited chunk of the file.
    /// The returned comment has the inital '>' trimmed.
    /// Returns false at eof.
    bool nextChunk( std::string& comment, std::vector<char>& bases );

    /// Ref args are set to the data from the next comma-delimited chunk of the file.
    /// The returned comment has the inital '>' trimmed.
    /// Ambiguous bases are silently replaced by a random base.
    /// Returns false at eof.
    bool nextChunk( std::string& comment, basevector& bases )
    { return nextChunk(comment,bases,0); }

    struct RandomizedBase
    {
        RandomizedBase( uint offset, char originalBase )
        : mOffset(offset), mOriginalBase(originalBase)
        {}

        uint mOffset;
        char mOriginalBase;
    };

    /// Ref args are set to the data from the next comma-delimited chunk of the file.
    /// The returned comment has the inital '>' trimmed.
    /// Ambiguous bases are replaced by a random base, and the offset and original
    /// value are returned.
    /// Returns false at eof.
    bool nextChunk( std::string& comment, basevector& bases, std::vector<RandomizedBase>& randomizedBases )
    { return nextChunk(comment,bases,&randomizedBases); }

    struct NullTrimmer
    {
        void operator()( basevector& /*bases*/ ) {}
    };

    struct LengthTrimmer
    {
        LengthTrimmer( int rightTrim, int truncLen )
        : mRightTrim(rightTrim), mTruncLen(truncLen)
        {}

        void operator()( basevector& bases )
        { bases.resize(max(0,min(mTruncLen,bases.isize()-mRightTrim))); }

        int mRightTrim;
        int mTruncLen;
    };

    struct Trimmer
    {
        Trimmer( int leftTrim, int rightTrim, int truncLen )
        : mLeftTrim(leftTrim), mTotTrim(leftTrim+rightTrim), mTruncLen(truncLen)
        {}

        void operator()( basevector& bases )
        {
            basevector tmp;
            tmp.SetToSubOf(bases,mLeftTrim,max(0,min(mTruncLen,bases.isize()-mTotTrim)));
            tmp.Swap(bases);
        }

        int mLeftTrim;
        int mTotTrim;
        int mTruncLen;
    };

    /// Reads sequence from a FASTA file, and creates a fastb file from it.
    /// Optionally creates a fastamb and/or a names file.
    /// Can also trim reads as they're written.
    template <class T>
    static void fastaToFastb( char const* fastaFile, char const* fastbFile, char const* fastambFile = 0, char const* namesFile = 0, T trimmer = T() );

private:
    bool nextChunk( std::string& comment, basevector& bases, std::vector<RandomizedBase>* pRandomizedBases );
};

/// Reads quality data from a FASTA file.
class QualParser : public FastaParser
{
public:
    QualParser( char const* fileName ) : FastaParser(fileName) {}

    QualParser( std::istream& is ) : FastaParser(is) {}

    /// Ref args are set to the data from the next comma-delimited chunk of the file.
    /// The returned comment has the inital '>' trimmed.
    /// Returns false at eof.
    bool nextChunk( std::string& comment, qualvector& quals );

    static void qualaToQualb( char const* qualaFile, char const* qualbFile );
};

/// Reads boolean data from a FASTA file.
class BitsParser : public FastaParser
{
public:
    BitsParser( char const* fileName ) : FastaParser(fileName) {}

    BitsParser( std::istream& is ) : FastaParser(is) {}

    /// Ref args are set to the data from the next comma-delimited chunk of the file.
    /// The returned comment has the inital '>' trimmed.
    /// Returns false at eof.
    bool nextChunk( std::string& comment, bitvector& bits );

    static void filterToVecbitvec( char const* filterFile, char const* vecbitvecFile );
};

template <class T>
void SeqParser::fastaToFastb( char const* fastaFile, char const* fastbFile, char const* fastambFile, char const* namesFile, T trimmer )
{
    SeqParser reader(fastaFile);
    IncrementalWriter<basevector> seqWriter(fastbFile);

    IncrementalWriter<bitvector>* pMBWriter = 0;
    std::vector<RandomizedBase>* pRBVec = 0;
    if ( fastambFile )
    {
        pMBWriter = new IncrementalWriter<bitvector>(fastambFile);
        pRBVec = new vector<RandomizedBase>();
    }

    IncrementalWriter<String>* pNamesWriter = 0;
    if ( namesFile )
        pNamesWriter = new IncrementalWriter<String>(namesFile);

    std::string comment;
    basevector bases(0,10000000);
    while ( reader.nextChunk(comment,bases,pRBVec) )
    {
        trimmer(bases);

        seqWriter.add(bases);

        if ( pNamesWriter )
            pNamesWriter->add(comment);

        if ( pMBWriter )
        {
            bitvector bits(bases.size());
            std::vector<RandomizedBase>::const_iterator end(pRBVec->end());
            for ( std::vector<RandomizedBase>::const_iterator itr(pRBVec->begin()); itr != end; ++itr )
                bits.Set(itr->mOffset,true);
            pMBWriter->add(bits);
        }
    }

    if ( pNamesWriter )
    {
        pNamesWriter->close();
        delete pNamesWriter;
    }
    if ( pMBWriter )
    {
        pMBWriter->close();
        delete pMBWriter;
        delete pRBVec;
    }
    seqWriter.close();
}

#endif /* UTIL_FASTAPARSER_H_ */
