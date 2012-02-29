/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file LookupTabBuilder.cc
 * \author tsharpe
 * \date Jun 16, 2009
 *
 * \brief
 */
#include "lookup/LookupTabBuilder.h"
#include "util/FastaParser.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/VirtualMasterVec.h"
#include <string.h>
#include <fstream>

LookupTabBuilder::LookupTabBuilder( unsigned int k ) :
    mK(k), mContigNo(0), mpPrevList(0)
{
    if ( k > 15 )
        FatalErr("k can be no larger than 15");

    size_type nListHeads = getNKmers();
    mListHeads = new unsigned int[nListHeads];
    memset(mListHeads, ~0, nListHeads * sizeof(unsigned int));

    // NOTE: This pre-reserves a lot of memory, regardless of the reference size
    // (currently 12 bytes per Entry * 2^28 Entries = 3GB)
    mLocations.reserve(1UL << 28);
}

LookupTabBuilder::~LookupTabBuilder()
{
    delete mpPrevList;
    delete[] mListHeads;
}

/// Add a new contig.
void LookupTabBuilder::add( bvec const& contig, progressReportFunc progressReport )
{
    if ( contig.size() >= mK )
    {
        bvec::const_iterator end(contig.End());
        bvec::const_iterator itr(contig.Begin());
        unsigned int kmer = 0U;
        for ( unsigned int iii = 0U; iii < mK; ++iii )
        {
            kmer = (kmer << 2) | *itr;
            ++itr;
        }
        unsigned int offset = 0U;
        unsigned int mask = getNKmers() - 1U;
        add(kmer, mContigNo, offset++);
        while ( itr != end )
        {
            kmer = (kmer << 2) | *itr;
            ++itr;
            add(kmer & mask, mContigNo, offset++);
        }
    }
    mContigNo += 1;
    progressReport(mContigNo, contig);
}

/// Add all the contigs in a fastb file.
void LookupTabBuilder::addFastb( char const* fileName, progressReportFunc progressReport )
{
    VirtualMasterVec<bvec> fvv(fileName);
    VirtualMasterVec<bvec>::const_iterator end(fvv.end());
    for ( VirtualMasterVec<bvec>::const_iterator itr(fvv.begin()); itr != end; ++itr )
        add(*itr, progressReport);
}

/// Add all the contigs in a fasta file.
void LookupTabBuilder::addFasta( char const* fileName, progressReportFunc progressReport )
{
    SeqParser sp(fileName);
    std::string comment;
    bvec seq;
    while ( sp.nextChunk(comment, seq) )
        add(seq, progressReport);
}

/// Write the lookup table.
void LookupTabBuilder::write( char const* fileName ) const
{
    unsigned int nKmers = getNKmers();
    IncrementalWriter<LocationVec> writer(fileName, nKmers);
    for ( unsigned int kmer = 0; kmer < nKmers; ++kmer )
    {
        std::list<Location> elements(getLocations(kmer));
        LocationVec vec;
        vec.reserve(elements.size());
        std::list<Location>::const_reverse_iterator end(elements.rend());
        for ( std::list<Location>::const_reverse_iterator itr(elements.rbegin()); itr != end; ++itr )
            vec.push_back(*itr);
        writer.add(vec);
    }
    writer.close();
}

void LookupTabBuilder::overflow()
{
    LookupTabBuilder* pNewList = new LookupTabBuilder(mK);

    pNewList->mpPrevList = mpPrevList;
    mpPrevList = pNewList;

    unsigned int* tmp = mListHeads;
    mListHeads = pNewList->mListHeads;
    pNewList->mListHeads = tmp;

    mLocations.swap(pNewList->mLocations);
}

void LookupTabBuilder::fill( unsigned int kmer, std::list<Location>& locations ) const
{
    if ( mpPrevList )
        mpPrevList->fill(kmer, locations);

    unsigned int entryIdx = mListHeads[kmer];
    while ( entryIdx != ~0U )
    {
        Entry const& entry = mLocations[entryIdx];
        locations.push_back(entry);
        entryIdx = entry.getPrevLocation();
    }
}
