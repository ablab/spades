/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file LookupTabBuilder.h
 * \author tsharpe
 * \date Jun 16, 2009
 *
 * \brief Facility for building new-style lookup tables.
 */
#ifndef LOOKUP_LOOKUPTABBUILDER_H_
#define LOOKUP_LOOKUPTABBUILDER_H_
#include "lookup/LookupTab.h"
#include "Basevector.h"
#include <vector>
#include <list>

class LookupTabBuilder
{
    class Entry : public Location
    {
    public:
        Entry( unsigned int contigID, unsigned int offset, unsigned int prevLocation )
        : Location(contigID,offset), mPrevLocation(prevLocation)
        {}

        // compiler-supplied copying and destructor are OK

        unsigned int getPrevLocation() const
        { return mPrevLocation; }

    private:
        unsigned int mPrevLocation;
    };

public:
    /// A function to call after each contig is added to the lookup table.
    typedef void (*progressReportFunc)( unsigned int contigNo, bvec const& bvec );

    LookupTabBuilder( unsigned int k = 12 );
    ~LookupTabBuilder();

    /// Add a new contig.
    void add( bvec const& contig, progressReportFunc progressReport=nullProgressReport );

    /// Add all the contigs in a fastb file.
    void addFastb( char const* fileName, progressReportFunc progressReport=nullProgressReport );

    /// Add all the contigs in a fasta file.
    void addFasta( char const* fileName, progressReportFunc progressReport=nullProgressReport );

    /// Retrieve a list of locations for a kmer.
    std::list<Location> getLocations( unsigned int kmer ) const
    { std::list<Location> result;
      fill(kmer,result);
      return result; }

    /// Write the lookup table.
    void write( char const* fileName ) const;

private:
    LookupTabBuilder( LookupTabBuilder const& ); // unimplemented -- no copying
    LookupTabBuilder& operator=( LookupTabBuilder const& ); // unimplemented -- no copying

    unsigned int getNKmers() const
    { return 1U << 2*mK; }

    // Adds a new location for a kmer.
    void add( unsigned int kmer, unsigned int contigID, unsigned int offset )
    { unsigned int entryIdx = mLocations.size();
      if ( entryIdx == mLocations.capacity() ) overflow();
      mLocations.push_back( Entry(contigID,offset,setListHead(kmer,entryIdx)) ); }

    // This table is full, chain a new one.
    void overflow();

    // Build singly-linked list of locations.
    unsigned int setListHead( unsigned int kmer, unsigned int entryIdx )
    { unsigned int result = mListHeads[kmer];
      mListHeads[kmer] = entryIdx;
      return result; }

    // Fill a list of locations for a kmer.
    void fill( unsigned int kmer, std::list<Location>& locations ) const;

    static void nullProgressReport( unsigned int /*contigNo*/, bvec const& )
    {}

    unsigned int mK; // kmer size
    unsigned int mContigNo; // number of contigs added so far
    LookupTabBuilder* mpPrevList; // each table is 1G Locations, if we need more we chain LookupTabBuilders
    unsigned int* mListHeads; // the head of a singly-linked list of LocationEntries for each kmer
    std::vector<Entry> mLocations;
};

#endif /* LOOKUP_LOOKUPTABBUILDER_H_ */
