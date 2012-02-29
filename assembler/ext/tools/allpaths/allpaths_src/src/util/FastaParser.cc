/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file FastaParser.cc
 * \author tsharpe
 * \date Apr 7, 2009
 *
 * \brief Several utility classes for reading FASTA files.
 */
#include "util/FastaParser.h"
#include "dna/Bases.h"
#include "feudal/IncrementalWriter.h"
#include <string.h>
#include <stdlib.h>
#include <sstream>

using std::string;
using std::vector;

char* FastaParser::readLine()
{
    // read until we get a non-empty line or hit eof
    while ( mIS.getline(mBuf,BUF_SIZE) )
    {
        if ( mBuf[0] )
            return mBuf;
    }

    if ( !mIS.eof() )
    { FatalErr("Error while reading fasta input stream."); }

    mBuf[0] = 0;
    return mBuf;
}

void FastaParser::init()
{
    if ( !mIS.good() )
    { FatalErr("Fasta input stream can't be read."); }

    readLine();
}

bool SeqParser::nextChunk( string& comment, vector<char>& bases )
{
    comment.clear();
    bases.clear();

    char* buf = getBuf();
    if ( !buf[0] )
        return false;  // EARLY RETURN!

    if ( buf[0] == '>' )
    {
        comment = string(buf+1);
        buf = readLine();
    }

    while ( buf[0] && buf[0] != '>' )
    {
        char* pBuf = buf;
        char* end = pBuf + strlen(buf);
        while ( pBuf < end )
        {
            char base = *pBuf++;

            if ( !GeneralizedBase::isGeneralizedBase(base) )
            { FatalErr("Fasta input stream contained the bogus character '" << base << "' at position " << pBuf-buf << " of line " << getLineNumber()); }

            bases.push_back(base);
        }
        buf = readLine();
    }

    return true;
}

bool SeqParser::nextChunk( string& comment, basevector& bases, vector<RandomizedBase>* pRandomizedBases )
{
    comment.clear();
    bases.clear();
    if ( pRandomizedBases )
        pRandomizedBases->clear();

    char* buf = getBuf();
    if ( !buf[0] )
        return false;  // EARLY RETURN!

    if ( buf[0] == '>' )
    {
        comment = string(buf+1);
        buf = readLine();
    }

    bool randomizerSeeded = false;
    while ( buf[0] && buf[0] != '>' )
    {
        char* pBuf = buf;
        char* end = pBuf + strlen(buf);

        while ( pBuf < end )
        {
            char chr = *pBuf++;
            if ( Base::isBase(chr) )
                bases.push_back(Base::char2Val(chr));
            else
            {
                if ( !GeneralizedBase::isGeneralizedBase(chr) )
                    FatalErr("Fasta input stream contained the bogus character '" << chr << "' at position " << pBuf-buf << " of line " << getLineNumber());

                // Seed the randomizer, when required, for each sequence by
                // computing a hash of the comment string.  This will
                // produce the same randomized sequence for any unchanged
                // file, but will not suffer the problem of always
                // producing the same "random" bases in the same order for
                // every sequence.
                if ( !randomizerSeeded )
                {
                    uint seed = bases.size();
                    size_t commentLen = comment.size();
                    for ( size_t idx = 0; idx < commentLen; ++idx )
                        seed += (seed * 47) + comment.at(idx);
                    srandom(seed);
                    randomizerSeeded = true;
                }
                if ( pRandomizedBases )
                {
                    pRandomizedBases->push_back(RandomizedBase(bases.size(),chr));
                }
                bases.push_back(GeneralizedBase::char2Val(chr));
            }
        }

        buf = readLine();
    }

    return true;
}

bool QualParser::nextChunk( std::string& comment, qualvector& quals )
{
    comment.clear();
    quals.resize(0U);

    char* buf = getBuf();
    if ( !buf[0] )
        return false;  // EARLY RETURN!

    if ( buf[0] == '>' )
    {
        comment = string(buf+1);
        buf = readLine();
    }

    while ( buf[0] && buf[0] != '>' )
    {
        istringstream ss(buf);
        int val;
        while ( ss >> val )
        {
            if ( val < 0 || val > 255 )
            { FatalErr("Invalid quality value " << val << " on line " << getLineNumber()); }

            quals.push_back(static_cast<uchar>(val));
        }

        if ( !ss.eof() )
        { FatalErr("Can't interpret quality values on line " << getLineNumber()); }

        buf = readLine();
    }

    return true;
}

void QualParser::qualaToQualb( char const* qualaFile, char const* qualbFile )
{
    QualParser reader(qualaFile);
    IncrementalWriter<qualvector> writer(qualbFile);

    std::string comment;
    qualvector quals;
    while ( reader.nextChunk(comment,quals) )
    {
        writer.add(quals);
    }

    writer.close();
}

bool BitsParser::nextChunk( std::string& comment, bitvector& bitsOut )
{
    comment.clear();
    bitvector bits(1000000);
    unsigned int bitCount = 0;

    char* buf = getBuf();
    if ( !buf[0] )
        return false;  // EARLY RETURN!

    if ( buf[0] == '>' )
    {
        comment = string(buf+1);
        buf = readLine();
    }

    while ( buf[0] && buf[0] != '>' )
    {
        istringstream ss(buf);
        bool val;
        while ( ss >> val )
        {
            if( bitCount == bits.size() )
                bits.resize(2*bitCount);
            bits.Set(bitCount++,val);
        }

        if ( !ss.eof() )
        { FatalErr("Can't interpret bit values on line " << getLineNumber()); }

        buf = readLine();
    }

    bits.resize(bitCount);
    bits.swap(bitsOut);
    return true;
}

void BitsParser::filterToVecbitvec( char const* filterFile, char const* vecbitvecFile )
{
    BitsParser reader(filterFile);
    IncrementalWriter<bitvector> writer(vecbitvecFile);

    std::string comment;
    bitvector bits;
    while ( reader.nextChunk(comment,bits) )
    {
        writer.add(bits);
    }

    writer.close();
}
