///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FragmentFiller.h
 * \author tsharpe
 * \date Oct 28, 2010
 *
 * \brief
 */
#ifndef FRAGMENTFILLER_H_
#define FRAGMENTFILLER_H_

#ifndef NDEBUG
#define LOG(x) std::cout << x << std::endl
#else
#define LOG(x)
#endif

#include "kmers/ReadPather.h"
#include "PairsManager.h"
#include "String.h"
#include <numeric>
#include <ostream>
#include <utility>

/// Little helper class to track observed pair separations.
/// Used to adjust library stats at the end of the run.
class LibSeps
{
public:
    LibSeps( PairsManager const& pm, int libId, int K, double maxStretch );
    LibSeps( LibSeps const& );
    ~LibSeps() { delete [] mBins; }

    LibSeps& operator=( LibSeps const& );
    LibSeps& operator+=( LibSeps const& );

    long getMinSep() const { return mMinSep; }
    long getMaxSep() const { return mMaxSep; }

    /// Check a separation measured in Kmers.
    bool testSep( long separation ) const
    { return mMinSep <= separation && separation <= mMaxSep; }

    /// Record a separation measured in Kmers.
    void postSep( long separation )
    { mBins[separation-mMinSep] += 1; }

    /// Look for greatest density across a 1.0*sigma window.
    /// Returns the number of BASES of separation at the center of this window.
    long getLibrarySeparation() const;

    /// Get revised lib stats assuming filled pairs have a normal distribution
    /// centered on the existing mean as reported by the PairsManager.
    std::pair<int,int> getLibStats() const;

    size_t getNFills() const
    { return std::accumulate(mBins,mBins+getNBins(),0ul); }

private:
    size_t getNBins() const { return mMaxSep - mMinSep + 1; }

    int mSD;
    int mK;
    long mMinSep;
    long mMaxSep;
    unsigned long* mBins;
};

class LibSepsColl
{
public:
    LibSepsColl( PairsManager& pm, int K, double maxStretch );
    LibSepsColl( LibSepsColl const& );
    ~LibSepsColl();

    LibSepsColl& operator+=( LibSepsColl const& );

    PairsManager const& getPM() const { return mPM; }
    LibSeps const& getLibSeps( int libId ) const { return *mColl[libId]; }
    LibSeps& getLibSeps( int libId ) { return *mColl[libId]; }

    void updatePM();

    size_t getNFills() const
    { size_t result = 0;
      LibSeps** ppSep = mColl; LibSeps** ppEnd = ppSep + mPM.nLibraries();
      for ( ; ppSep != ppEnd; ++ppSep )
          if ( *ppSep ) result += ppSep[0]->getNFills();
      return result; }

private:
    LibSepsColl& operator=( LibSepsColl const& ); // unimplemented -- no copying

    PairsManager& mPM;
    LibSeps** mColl;
};

class FragmentFillerInfo
{
public:
    FragmentFillerInfo()
    : mReadId(-1), mIsRc(false) {}

    FragmentFillerInfo( long readId, bool isRc )
    : mReadId(readId), mIsRc(isRc) {}

    // compiler-supplied copying and destructor are OK

    long getReadId() const { return mReadId; }
    bool isRC() const { return mIsRc; }

private:
    long mReadId;
    bool mIsRc;
};

class FragmentFillInfo : public FragmentFillerInfo
{
public:
    FragmentFillInfo( bvec const& fill, size_t pairId )
    : mFill(fill), mPairId(pairId) {}

    FragmentFillInfo( FragmentFillerInfo const& filler,
                        bvec const& fill, size_t pairId )
      : FragmentFillerInfo(filler), mFill(fill), mPairId(pairId) {}

    // compiler-supplied copying and destructor are OK

    bvec const& getFill() const { return mFill; }
    size_t getPairId() const { return mPairId; }

    ostream& write( ostream& os, size_t nFrags ) const
    { os << mPairId << ' ' << nFrags << ' ' << mFill.size() << ' '
        << getReadId() << ' ' << isRC();
      return os; }

    friend int compare( FragmentFillInfo const& ffi1,
                        FragmentFillInfo const& ffi2 )
    { int result = compare(ffi1.mPairId,ffi2.mPairId);
      if ( !result ) result = compare(ffi1.mFill,ffi2.mFill);
      if ( !result ) result = compare(ffi1.getReadId(),ffi2.getReadId());
      if ( !result ) result = compare(ffi1.isRC(),ffi2.isRC());
      return result; }

    friend bool operator<( FragmentFillInfo const& ffi1,
                           FragmentFillInfo const& ffi2 )
    { return compare(ffi1,ffi2) < 0; }

private:
    bvec mFill;
    size_t mPairId;
};

TRIVIALLY_SERIALIZABLE(FragmentFillInfo);

template <unsigned K>
class FragmentFiller
{
public:
    typedef std::vector<FragmentFillInfo> OutVec;

    FragmentFiller( String const& pathInfo, String const& graphInfo,
                    long minOverlap, size_t maxCliq, size_t maxGlue,
                    bool allClosures );

    // compiler-supplied destructor is OK

    /// fill all pairs, write closures to the given file as a vecbvec
    size_t processPairs( LibSepsColl& lsc, String const& filledReadsFilename,
                            size_t nThreads );

    /// fill a specific pair and append closures to vecbvec
    void processPair( LibSepsColl& lsc, size_t pairId, OutVec* pOut );

    typedef EdgeList::const_iterator ELItr;
    typedef EdgeList const ELC;
    struct Matchup
    {
        ELC* pRead;
        ELC* pFill;
        ELItr begFill;
        ELItr endFill;
        bool isRC;

        friend bool operator==( Matchup const& m1, Matchup const& m2 )
        { return m1.begFill == m2.begFill && m1.endFill == m2.endFill &&
                  m1.isRC == m2.isRC; }
        friend bool operator!=( Matchup const& m1, Matchup const& m2 )
        { return !(m1 == m2); }
    };

private:
    FragmentFiller( FragmentFiller const& ); // unimplemented -- no copying
    FragmentFiller& operator=( FragmentFiller const& ); // unimplemented -- no copying

    // attempts to fill the pair <el1,el2> by looking for other read pairs that
    // overlap.  if successful, it writes the closures to pVBV.
    void findGlue( size_t pairId, ELC& el1, size_t descIdx1, ELC& el2,
                    size_t descIdx2, PairsManager const& pm, LibSeps* pLS,
                    OutVec* pOut );

    // checks that a proposed fill meets the stretch criterion
    bool checkSeparation( Matchup const& m1, Matchup const& m2,
                          LibSeps const& ls, EdgeList* pFill ) const;

    // checks that the glue-read's mate pair is consistent with the fill
    bool checkMate( bool isRC, size_t readId, ELC& fill ) const
    { EdgeList rd = mPaths.getEdgeList(readId);
      return isRC? checkLeadingMate(rd,fill): checkTrailingMate(rd.rc(),fill); }
    bool checkLeadingMate( ELC& read, ELC& fill ) const;
    bool checkTrailingMate( ELC& read, ELC& fill ) const;

    long separation( ELItr itr1, ELItr itr2, long ends ) const
    { if ( itr1 < itr2 ) { ends -= mPaths.getEdgeListLen(itr1,itr2); }
      if ( itr1 > itr2 ) { ends += mPaths.getEdgeListLen(itr2,itr1); }
      return ends; }

    PathsWithEvidence<K> mPaths;
    long mMinOverlap; // minimum no. of kmers overlap for valid glue
    size_t mMaxCliq;
    size_t mMaxGlue;
    bool mAllClosures;
};

#endif /* FRAGMENTFILLER_H_ */
