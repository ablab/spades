///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FragmentFillerDefs.h
 * \author tsharpe
 * \date Nov 22, 2010
 *
 * \brief
 */
#ifndef FRAGMENTFILLERDEFS_H_
#define FRAGMENTFILLERDEFS_H_

#include "paths/FragmentFiller.h"
#include "feudal/IncrementalWriter.h"
#include "system/Assert.h"
#include "system/SortInPlace.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include "Basevector.h"
#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <time.h>

namespace
{

template <unsigned K>
class EdgeListKeeper
{
public:
    EdgeListKeeper( PathsWithEvidence<K> const& paths )
    : mPaths(paths) {}

    // compiler-supplied copying and destructor is OK

    ReadID const& getReadID() const { return mReadID; }

    void setReadID( ReadID const& readID )
    { if ( mReadID != readID )
      { mReadID = readID; mFwd.clear(); mRC.clear(); } }

    EdgeList const& getFwd()
    { if ( !mFwd.size() ) mFwd = mPaths.getEdgeList(mReadID); return mFwd; }

    EdgeList const& getRC()
    { if ( !mRC.size() ) { mRC = getFwd(); mRC.rc(); } return mRC; }

    bool isPalindrome( EdgeID const& edgeID ) const
    { return mPaths.getGraph().getEdge(edgeID).isPalindrome(); }

    UnipathEvidenceVec const& getEvidence( EdgeID const& edgeID ) const
    { return mPaths.getEvidence(edgeID); }

    size_t getEdgeListLen( EdgeList::const_iterator itr,
                           EdgeList::const_iterator end )
    { return mPaths.getEdgeListLen(itr,end); }

    size_t getEdgeLen( EdgeID const& edgeID ) const
    { return mPaths.getEdgeLen(edgeID); }

private:
    PathsWithEvidence<K> const& mPaths;
    ReadID mReadID;
    EdgeList mFwd;
    EdgeList mRC;
};

template <unsigned K>
class Matcher
{
public:
    Matcher( EdgeListKeeper<K>& elk, EdgeList const& el, EdgeID const& edgeID,
                long minOverlap )
    : mpELK(&elk), mpEL(&el), mMinOverlap(minOverlap),
      mIsPalindromic(elk.isPalindrome(edgeID))
    { UnipathEvidenceVec const& uev = mpELK->getEvidence(edgeID);
      mItr = uev.begin();
      mEnd = uev.end();
      mLast = mDescIndices;
      for ( size_t idx = 0; idx < el.size(); ++idx )
          if ( el[idx].getEdgeID() == edgeID )
              *mLast++ = idx; }

    Matcher( Matcher const& that )
    { *this = that; }

    Matcher& operator=( Matcher const& that )
    { memcpy(this,&that,sizeof(Matcher));
      mLast = mDescIndices + (that.mLast-that.mDescIndices);
      return *this; }

    // compiler-supplied destructor is OK

    bool isDone() const { return mItr == mEnd; }

    EdgeID const& getEdgeID() const
    { return mpEL[0][mDescIndices[0]].getEdgeID(); }

    ReadID const& getCurReadID() const { return mItr->getReadID(); }
    void advancePastReadID( ReadID const& readID )
    { while ( mItr != mEnd && mItr->getReadID() == readID ) ++mItr; }

    bool leftMatchesUnambiguously();
    bool rightMatchesUnambiguously();

    typename FragmentFiller<K>::Matchup const& getMatchup() const
    { return mMatchup; }

    bool isMatchRC() const { return mMatchup.isRC; }

    size_t getCliqs() const
    { return (mLast-mDescIndices) * (mEnd-mItr) /
                mpELK->getEdgeLen(getEdgeID()); }

private:
    typedef EdgeList::const_iterator ELItr;

    bool leftMatch( EdgeList const& elFill, size_t idxFill, size_t idxRd )
    { bool result = false;
      size_t nTrailingEdgesFill = elFill.size() - idxFill;
      size_t nTrailingEdgesRd = mpEL->size() - idxRd;
      if ( nTrailingEdgesFill > nTrailingEdgesRd ||
           (nTrailingEdgesFill == nTrailingEdgesRd &&
              elFill.getFinalSkip() < mpEL->getFinalSkip()) )
      { size_t nLeadingEdges = std::min(idxFill,idxRd);
        ELItr iFill = elFill.begin() + idxFill;
        ELItr begFill = iFill-nLeadingEdges;
        ELItr begRd = mpEL->begin()+idxRd-nLeadingEdges;
        if ( std::equal(begRd,mpEL->end(),begFill) )
        { ELItr endFill = iFill+nTrailingEdgesRd;
          if ( checkLeftOverlap(begFill,endFill,elFill,idxFill,idxRd) )
          { mMatchup.begFill = begFill;
            mMatchup.endFill = endFill;
            mMatchup.pFill = &elFill;
            mMatchup.pRead = mpEL;
            result = true; } } }
        return result; }

    bool checkLeftOverlap( ELItr const& begFill, ELItr const& endFill,
                            EdgeList const& elFill,size_t idxFill,size_t idxRd )
    { if ( mMinOverlap <= 1 && endFill-begFill > 1 )
        return true;
      long overlap = mpELK->getEdgeListLen(begFill,endFill);
      overlap -= mpEL->getFinalSkip();
      if ( idxFill < idxRd ) overlap -= elFill.getInitialSkip();
      else if ( idxRd < idxFill ) overlap -= mpEL->getInitialSkip();
      else overlap -= std::max(elFill.getInitialSkip(),mpEL->getInitialSkip());
      return overlap >= mMinOverlap; }

    bool rightMatch( EdgeList const& elFill, size_t idxFill, size_t idxRd )
    { bool result = false;
      if ( idxFill > idxRd ||
          (idxFill == idxRd && elFill.getInitialSkip()<mpEL->getInitialSkip()) )
      { size_t nTrailingEdges = std::min(elFill.size()-idxFill,
                                          mpEL->size()-idxRd);
        ELItr iFill = elFill.begin() + idxFill;
        ELItr begFill = iFill-idxRd;
        ELItr endRd = mpEL->begin()+idxRd+nTrailingEdges;
        if ( std::equal(mpEL->begin(),endRd,begFill) )
        { ELItr endFill = iFill+nTrailingEdges;
          if ( checkRightOverlap(begFill,endFill,elFill,idxFill,idxRd) )
          { mMatchup.begFill = begFill;
            mMatchup.endFill = endFill;
            mMatchup.pFill = &elFill;
            mMatchup.pRead = mpEL;
            result = true; } } }
        return result; }


    bool checkRightOverlap( ELItr const& begFill, ELItr const& endFill,
                            EdgeList const& elFill,size_t idxFill,size_t idxRd )
    { if ( mMinOverlap <= 1 && endFill-begFill > 1 )
        return true;
      long overlap = mpELK->getEdgeListLen(begFill,endFill);
      overlap -= mpEL->getInitialSkip();
      idxFill = elFill.size() - idxFill;
      idxRd = mpEL->size() - idxRd;
      if ( idxFill < idxRd ) overlap -= elFill.getFinalSkip();
      else if ( idxRd < idxFill ) overlap -= mpEL->getFinalSkip();
      else overlap -= std::max(elFill.getFinalSkip(),mpEL->getFinalSkip());
      return overlap >= mMinOverlap; }

    typedef UnipathEvidenceVec::const_iterator Itr;
    EdgeListKeeper<K>* mpELK;
    EdgeList const* mpEL;
    long mMinOverlap;
    bool mIsPalindromic;
    Itr mItr;
    Itr mEnd;
    size_t mDescIndices[UnipathEvidence::MAX_NSEGMENTS];
    size_t* mLast;
    typename FragmentFiller<K>::Matchup mMatchup;
};

template <unsigned K>
bool Matcher<K>::leftMatchesUnambiguously()
{
    bool result = false;
    mpELK->setReadID(mItr->getReadID());
    typename FragmentFiller<K>::Matchup firstMatch;
    do
    {
        size_t idxFill = mItr->getSegmentID();
        EdgeList const& elFill = mpELK->getFwd();
        EdgeDesc const& edFill = elFill[idxFill];
        for ( size_t* ppp = mDescIndices; ppp < mLast; ++ppp )
        {
            size_t idxRd = *ppp;
            EdgeDesc const& edRd = mpEL[0][idxRd];
            if ( mIsPalindromic || edRd.getStatus() == edFill.getStatus() )
            {
                if ( leftMatch(elFill,idxFill,idxRd) )
                {
                    if ( result && mMatchup != firstMatch )
                        return false; // EARLY RETURN for ambiguous result
                    result = true;
                    mMatchup.isRC = false;
                    firstMatch = mMatchup;
                }
            }
            if ( mIsPalindromic || edRd.getStatus() != edFill.getStatus() )
            {
                EdgeList const& elFillRC = mpELK->getRC();
                size_t idxFillRC = elFillRC.size() - idxFill - 1;
                if ( leftMatch(elFillRC,idxFillRC,idxRd) )
                {
                    if ( result && mMatchup != firstMatch )
                        return false; // EARLY RETURN for ambiguous result
                    result = true;
                    mMatchup.isRC = true;
                    firstMatch = mMatchup;
                }
            }
        }
    }
    while ( mItr != mEnd && (++mItr)->getReadID() == mpELK->getReadID() );
    return result;
}

template <unsigned K>
bool Matcher<K>::rightMatchesUnambiguously()
{
    bool result = false;
    mpELK->setReadID(mItr->getReadID());
    typename FragmentFiller<K>::Matchup firstMatch;
    do
    {
        size_t idxFill = mItr->getSegmentID();
        EdgeList const& elFill = mpELK->getFwd();
        EdgeDesc const& edFill = elFill[idxFill];
        for ( size_t* ppp = mDescIndices; ppp < mLast; ++ppp )
        {
            size_t idxRd = *ppp;
            EdgeDesc const& edRd = mpEL[0][idxRd];
            if ( mIsPalindromic || edRd.getStatus() == edFill.getStatus() )
            {
                if ( rightMatch(elFill,idxFill,idxRd) )
                {
                    if ( result && mMatchup != firstMatch )
                        return false; // EARLY RETURN for ambiguous result
                    result = true;
                    mMatchup.isRC = false;
                    firstMatch = mMatchup;
                }
            }
            if ( mIsPalindromic || edRd.getStatus() != edFill.getStatus() )
            {
                EdgeList const& elFillRC = mpELK->getRC();
                size_t idxFillRC = elFillRC.size() - idxFill - 1;
                if ( rightMatch(elFillRC,idxFillRC,idxRd) )
                {
                    if ( result && mMatchup != firstMatch )
                        return false; // EARLY RETURN for ambiguous result
                    result = true;
                    mMatchup.isRC = true;
                    firstMatch = mMatchup;
                }
            }
        }
    }
    while ( mItr != mEnd && (++mItr)->getReadID() == mpELK->getReadID() );
    return result;
}

template <unsigned K> class Processor;

template <unsigned K>
class Common : SpinLockedData
{
public:
    typedef typename FragmentFiller<K>::OutVec OutVec;

    Common( LibSepsColl& lsc, String const& filledFragsFilename )
    : mLSC(lsc), mNFrags(0), mNDots(0), mNextBatchToWrite(0),
      mWriter(filledFragsFilename.c_str(),getNPairs()),
      mInfoOS(filledFragsFilename.ReplaceExtension(".fastb",".info").c_str())
    { std::cout << Date() << " There are " << getNPairs() <<
                " read pairs to process." << std::endl; }

    ~Common()
    { ForceAssertEq(mNextBatchToWrite,getNBatches());
      mWriter.close();
      mInfoOS.close();
      size_t nPairs = getNPairs();
      std::cout << std::endl;
      std::cout << '\n' << Date() << " Filled " <<
              100.*mWriter.getNElements()/nPairs << "% of " << nPairs <<
              " pairs." << std::endl; }

    size_t getNPairs() const { return mLSC.getPM().nPairs(); }
    size_t getNBatches() const { return (getNPairs()+BATCH_SIZE-1)/BATCH_SIZE; }
    size_t getNClosures() const { return mWriter.getNElements(); }

    void dump( size_t batchNo, OutVec& output )
    { SpinLocker locker(*this);
      if ( batchNo != mNextBatchToWrite )
      { OutVec& outVec = mWriteMap[batchNo];
        outVec.reserve(output.capacity());
        outVec.swap(output); }
      else
      { doDump(output);
        typedef typename map<size_t,OutVec>::iterator Itr;
        Itr itr(mWriteMap.find(++mNextBatchToWrite));
        while ( itr != mWriteMap.end() )
        { doDump(itr->second);
          mWriteMap.erase(itr);
          itr = mWriteMap.find(++mNextBatchToWrite); } }
      dot(); }

private:
    Common( Common const& ); // unimplemented -- no copying
    Common& operator=( Common const& ); // unimplemented -- no copying

    void doDump( OutVec const& outVec )
    { typedef typename OutVec::const_iterator Itr;
      for ( Itr itr(outVec.begin()), end(outVec.end()); itr != end; ++itr )
      { mWriter.add(itr->getFill());
        itr->write(mInfoOS,mNFrags++) << '\n'; } }

    void dot()
    { std::cout << '.' << std::flush;
      if ( !((++mNDots)%100) ) std::cout << std::endl; }

    LibSepsColl& mLSC;
    size_t mNFrags;
    size_t mNDots;
    size_t mNextBatchToWrite;
    map<size_t,OutVec> mWriteMap;
    IncrementalWriter<bvec> mWriter;
    ofstream mInfoOS;

    static size_t const BATCH_SIZE = 100000;
    friend class Processor<K>;
};

template <unsigned K>
class Processor
{
public:
    Processor( FragmentFiller<K>& ff, Common<K>& common )
    : mFF(ff), mLSC(common.mLSC), mCommon(common)
    {}

    Processor( Processor const& that )
    : mFF(that.mFF), mLSC(that.mLSC), mCommon(that.mCommon)
    { mOutput.reserve(Common<K>::BATCH_SIZE); }

    ~Processor()
    { SpinLocker locker(mCommon); mCommon.mLSC += mLSC; }

    void operator()( size_t batchNo );

private:
    FragmentFiller<K>& mFF;
    LibSepsColl mLSC;
    Common<K>& mCommon;
    typename FragmentFiller<K>::OutVec mOutput;
};

template <unsigned K>
void Processor<K>::operator()( size_t batchNo )
{
    size_t pairId = batchNo * Common<K>::BATCH_SIZE;
    size_t endId = std::min(pairId+Common<K>::BATCH_SIZE,mCommon.getNPairs());
    while ( pairId != endId )
        mFF.processPair(mLSC, pairId++, &mOutput);

    mCommon.dump(batchNo,mOutput);
    mOutput.clear();
}


}

LibSeps::LibSeps( PairsManager const& pm, int libId, int K, double maxStretch )
 : mK(K)
{
    int libSep = pm.getLibrarySep(libId) + K - 1;
    mSD = pm.getLibrarySD(libId);
    int delta = maxStretch*mSD + .5;
    mMinSep = libSep - delta;
    mMaxSep = libSep + delta;
    size_t nBins = getNBins();
    mBins = new unsigned long[nBins];
    memset(mBins, 0, sizeof(unsigned long)*nBins);
}

LibSeps::LibSeps( LibSeps const& that )
 : mSD(that.mSD), mK(that.mK), mMinSep(that.mMinSep), mMaxSep(that.mMaxSep)
{
    size_t nBins = getNBins();
    mBins = new unsigned long[nBins];
    memcpy(mBins, that.mBins, sizeof(unsigned long)*nBins);
}

LibSeps& LibSeps::operator=( LibSeps const& that )
{
    size_t nBins = that.getNBins();
    if ( getNBins() != nBins )
    {
        delete[] mBins;
        mBins = 0;
    }
    mSD = that.mSD;
    mK = that.mK;
    mMinSep = that.mMinSep;
    mMaxSep = that.mMaxSep;
    if ( !mBins )
    {
        mBins = new unsigned long[nBins];
        memcpy(mBins, that.mBins, sizeof(unsigned long)*nBins);
    }
    return *this;
}

LibSeps& LibSeps::operator+=( LibSeps const& that )
{
    size_t nBins = getNBins();
    ForceAssertEq(nBins,that.getNBins());
    unsigned long* itr = mBins;
    unsigned long* end = mBins + nBins;
    unsigned long* src = that.mBins;
    while ( itr != end )
        *itr++ += *src++;
    return *this;
}

long LibSeps::getLibrarySeparation() const
{
    unsigned long result = getNBins() / 2;
    size_t nBins = getNBins();
    if ( nBins > static_cast<size_t>(mSD) )
    {
        typedef unsigned long* Itr;
        Itr beg = mBins;
        Itr end = mBins + nBins;
        Itr stop = beg + mSD;
        unsigned long sum = 0;
        for ( Itr itr(beg); itr != stop; ++itr )
            sum += *itr;
        unsigned long bestSum = sum;
        result = mSD / 2;
        while ( stop != end )
        {
            sum -= *beg++;
            sum += *stop++;
            if ( sum > bestSum )
            {
                bestSum = sum;
                result = (beg - mBins) + mSD / 2;
            }
        }
    }
    return result + mMinSep - mK + 1;
}

std::pair<int,int> LibSeps::getLibStats() const
{
    double sum = 0;
    double sum2 = 0;
    unsigned long nnn = 0;
    unsigned long* ppp = mBins;
    long stop = mMaxSep - mK + 1;
    for ( long idx = mMinSep - mK + 1; idx <= stop; ++idx )
    {
        long count = *ppp++;
        nnn += count;
        sum += idx * count;
        sum2 += idx*idx * count;
    }

    int mean = -1;
    int sd = -1;
    if ( nnn )
    {
        mean = round(sum/nnn);
        sd = round(sqrt((sum2-sum*sum/nnn)/nnn));
    }
    return std::pair<int,int>(mean,sd);
}

LibSepsColl::LibSepsColl( PairsManager& pm, int K, double maxStretch )
: mPM(pm), mColl(0)
{
    size_t nLibs = pm.nLibraries();
    mColl = new LibSeps*[nLibs];
    for ( size_t libId = 0; libId != nLibs; ++libId )
        mColl[libId] = new LibSeps(pm, libId, K, maxStretch);
    pm.getPairsIndex(); // need to pre-load this, because it's not thread-safe
}

LibSepsColl::LibSepsColl( LibSepsColl const& that )
: mPM(that.mPM)
{
    size_t nLibs = mPM.nLibraries();
    mColl = new LibSeps*[nLibs];
    LibSeps** src(that.mColl);
    LibSeps** end(mColl+nLibs);
    for ( LibSeps** itr(mColl); itr != end; ++itr, ++src )
        *itr = new LibSeps(**src);
}

LibSepsColl::~LibSepsColl()
{
    LibSeps** end(mColl+mPM.nLibraries());
    for ( LibSeps** itr(mColl); itr != end; ++itr )
        delete *itr;
    delete[] mColl;
}

LibSepsColl& LibSepsColl::operator+=( LibSepsColl const& that )
{
    ForceAssertEq(&mPM,&that.mPM);
    LibSeps** src(that.mColl);
    LibSeps** end(mColl+mPM.nLibraries());
    for ( LibSeps** itr(mColl); itr != end; ++itr, ++src )
        **itr += **src;
    return *this;
}

void LibSepsColl::updatePM()
{
    size_t nLibs = mPM.nLibraries();
    for ( size_t libId = 0; libId < nLibs; ++libId )
    {
        int sep = getLibSeps(libId).getLibrarySeparation();
        int sd = mPM.getLibrarySD(libId);
        mPM.changeLibrarySepSd(libId, sep, sd);
        //std::pair<int,int> libStats = getLibSeps(libId).getLibStats();
        //mPM.changeLibrarySepSd(libId,libStats.first,libStats.second);
    }
}

template <unsigned K>
FragmentFiller<K>::FragmentFiller( String const& pathInfo,
                                    String const& graphInfo,
                                    long minOverlap, size_t maxCliq,
                                    size_t maxGlue, bool allClosures )
: mPaths(pathInfo,graphInfo), mMinOverlap(minOverlap), mMaxCliq(maxCliq),
  mMaxGlue(maxGlue), mAllClosures(allClosures)
{
    ForceAssertGt(mMinOverlap,0l);
}

template <unsigned K>
size_t FragmentFiller<K>::processPairs( LibSepsColl& lsc,
                                        String const& filledReadsFilename,
                                        size_t nThreads )
{
    Common<K> common(lsc,filledReadsFilename);
    Processor<K> proc(*this,common);
    if ( true )
        WorklistN<Processor<K> > wl(proc,common.getNBatches(),nThreads-1);
    return common.getNClosures();
}

template <unsigned K>
void FragmentFiller<K>::processPair( LibSepsColl& lsc, size_t pairId,
                                        OutVec* pOut )
{
    PairsManager const& pm = lsc.getPM();
    size_t readID1 = pm.ID1(pairId);
    size_t readID2 = pm.ID2(pairId);
    LibSeps& ls = lsc.getLibSeps(pm.libraryID(pairId));
    EdgeList el1 = mPaths.getEdgeList(readID1);
    EdgeList el2 = mPaths.getEdgeList(readID2);
    el2.rc();

    LOG("Working on " << pairId << ": " << el1 << " + " << el2);

    if ( !el1.size() || !el2.size() )
    {
        LOG("Giving up:  One or both of the ends are unpathed.");
        return; // HEY!  Non-structured code!
    }

    // if there's a possibility that an edge not represented in either read
    // could appear between read1 and read2 (i.e., if the number of kmers
    // skipped in the final edge of read1's path plus the number skipped in
    // the first edge of read2's path don't exceed the maximum separation),
    // then just go to look for glue immediately.
    long sep = el1.getFinalSkip() + el2.getInitialSkip();
    if ( sep <= ls.getMaxSep() )
    {
        LOG("Need glue:  pretty close to ends.");
        findGlue(pairId, el1, el1.size()-1, el2, 0, pm, &ls, pOut);
        return; // HEY!  Non-structured code!
    }

    // see if there is an unambiguous way of pairing the reads
    // just from their own evidence.

    ELItr beg1(el1.begin());
    ELItr end1(el1.end());
    ELItr beg2(el2.begin());
    ELItr end2(el2.end());
    ELItr matchPos(end1);
    long matchSep = 0;
    bool isAmbiguous = false;
    ELItr itr(end1 - std::min(el1.size(),el2.size()));
    sep -= mPaths.getEdgeListLen(itr,end1);
    while ( itr != end1 )
    {
        if ( ls.testSep(sep) && std::equal(itr, end1, beg2) )
        {
            if ( matchPos == end1 )
            {
                matchPos = itr;
                matchSep = sep;
            }
            else
            {
                isAmbiguous = true;
                break;
            }
        }
        sep += mPaths.getEdgeLen(*itr);
        ++itr;
    }

    if ( matchPos == end1 ) // if there are no end-matches, we give up
    {
        LOG("Giving up:  No end matches.");
    }
    else if ( !isAmbiguous ) // found a single overlap so we know the result
    {
        EdgeList fill(el1.getInitialSkip(), el2.getFinalSkip());
        fill.reserve((matchPos - beg1) + (end2 - beg2));
        fill.insert(fill.end(), beg1, matchPos);
        fill.insert(fill.end(), beg2, end2);
        if ( mPaths.getEdgeListLen(fill.begin(),fill.end()) <=
                fill.getInitialSkip() + fill.getFinalSkip() )
            LOG("Giving up:  Non-overlapping outie.");
        else
        {
            if ( pOut )
                pOut->push_back(FragmentFillInfo(mPaths.getBases(fill),pairId));
            ls.postSep(matchSep);
            LOG("Have an unambiguous overlap.  We're done.");
        }
    }
    else if ( matchPos != beg1 && end1 - matchPos
            < static_cast<long> (el2.size()) )
    {
        LOG("Need glue:  ambiguous overlap.");
        findGlue(pairId,el1,matchPos-beg1-1,el2,end1-matchPos,pm,&ls,pOut);
    }
    else // one or the other read is totally ambiguous
    {
        LOG("Giving up:  totally ambiguous read.");
    }
}

template <unsigned K>
void FragmentFiller<K>::findGlue( size_t pairId, ELC& el1, size_t descIdx1,
                                  ELC& el2, size_t descIdx2,
                                  PairsManager const& pm, LibSeps* pLS,
                                  OutVec* pOut )
{
    EdgeListKeeper<K> elk(mPaths);
    EdgeID const& edgeID1 = el1[descIdx1].getEdgeID();
    Matcher<K> m1(elk,el1,edgeID1,mMinOverlap);
    size_t cliqs1 = m1.getCliqs();
    if ( descIdx1 && cliqs1 > mMaxCliq/10 && mPaths.getEdgeLen(edgeID1) <= 3 )
    {
        Matcher<K> tmp(elk,el1,el1[descIdx1-1].getEdgeID(),mMinOverlap);
        size_t cliqsTmp = tmp.getCliqs();
        if ( cliqs1 > mMaxCliq || 2*cliqsTmp < cliqs1 )
        {
            m1 = tmp;
            cliqs1 = cliqsTmp;
        }
    }

    EdgeID const& edgeID2 = el2[descIdx2].getEdgeID();
    Matcher<K> m2(elk,el2,edgeID2,mMinOverlap);
    size_t cliqs2 = m2.getCliqs();
    if ( descIdx2 +1 < el2.size() && cliqs2 > mMaxCliq/10 &&
            mPaths.getEdgeLen(edgeID2) <= 3 )
    {
        Matcher<K> tmp(elk,el2,el2[descIdx2+1].getEdgeID(),mMinOverlap);
        size_t cliqsTmp = tmp.getCliqs();
        if ( cliqs2 > mMaxCliq || 2*cliqsTmp < cliqs2 )
        {
            m2 = tmp;
            cliqs2 = cliqsTmp;
        }
    }
    if ( cliqs1 > mMaxCliq || cliqs2 > mMaxCliq )
    {
        LOG("Too much potential glue.");
        return; //EARLY RETURN!!!
    }
    map<EdgeList,FragmentFillerInfo> closures;
    size_t nGlueReads = 0;
    EdgeList fill;

    while ( !m1.isDone() && !m2.isDone() )
    {
        ReadID const& readID1 = m1.getCurReadID();
        ReadID const& readID2 = m2.getCurReadID();
        if ( readID1 == readID2 )
        {
            long readId = readID1.val();
            long pairId = pm.getPairID(readId);
            if ( pairId != -1 &&
                 m1.leftMatchesUnambiguously() &&
                 m2.rightMatchesUnambiguously() &&
                 m1.getMatchup().isRC == m2.getMatchup().isRC &&
                 checkSeparation(m1.getMatchup(),m2.getMatchup(),*pLS,&fill) &&
                 checkMate(m1.isMatchRC(),pm.getPartner(pairId,readId),fill) )
            {
                closures[fill]=FragmentFillerInfo(readId,m1.getMatchup().isRC);
                LOG("Closure by " << readId << ": " << elk.getFwd());
                if ( ++nGlueReads >= mMaxGlue ||
                     (!mAllClosures && closures.size() > 1) )
                    break;
            }
        }
        if ( readID1 <= readID2 )
            m1.advancePastReadID(readID1);
        if ( readID2 <= readID1 )
            m2.advancePastReadID(readID2);
    }

    size_t nClosures = closures.size();
    if ( nClosures == 1 || (mAllClosures && nClosures > 0) )
    {
        typedef map<EdgeList,FragmentFillerInfo>::iterator Itr;
        Itr end(closures.end());
        size_t len = mPaths.getEdgeListLen(el1) + mPaths.getEdgeListLen(el2);
        if ( !pOut )
        {
            for ( Itr itr(closures.begin()); itr != end; ++itr )
                pLS->postSep(mPaths.getEdgeListLen(itr->first) - len);
        }
        else
        {
            size_t startingSize = pOut->size();
            for ( Itr itr(closures.begin()); itr != end; ++itr )
            {
                pOut->push_back(FragmentFillInfo(itr->second,
                                                 mPaths.getBases(itr->first),
                                                 pairId));
                pLS->postSep(mPaths.getEdgeListLen(itr->first) - len);
            }
            if ( nClosures > 1 )
                sortInPlace(pOut->begin()+startingSize,pOut->end());
        }
    }

    LOG("Found " << closures.size() << " unique closures supported by " <<
            nGlueReads << " reads using " << cliqs1 << " cliqs on edge " <<
            m1.getEdgeID() << " and " << cliqs2 << " cliqs on edge " <<
            m2.getEdgeID());
}

template <unsigned K>
bool FragmentFiller<K>::checkSeparation( Matchup const& m1, Matchup const& m2,
                                      LibSeps const& ls, EdgeList* pFill ) const
{
    bool result = false;
    if ( m1.begFill <= m2.endFill )
    {
        ELC& el1 = *m1.pRead;
        ELC& el2 = *m2.pRead;
        long sep = separation(m2.begFill,m1.endFill,
                                el1.getFinalSkip() + el2.getInitialSkip());
        if ( ls.testSep(sep) )
        {
            EdgeList fill(el1.getInitialSkip(), el2.getFinalSkip());
            ELItr beg1 = el1.begin();
            ELItr end1 = el1.end() - (m1.endFill - m1.begFill);
            ELItr begF = m1.begFill;
            ELItr endF = m2.endFill;
            ELItr beg2 = el2.begin() + (m2.endFill - m2.begFill);
            ELItr end2 = el2.end();
            Assert(end1 >= beg1);
            Assert(endF >= begF);
            Assert(end2 >= beg2);
            fill.reserve((end1-beg1) + (endF-begF) + (end2-beg2));
            fill.insert(fill.end(), beg1, end1);
            fill.insert(fill.end(), begF, endF);
            fill.insert(fill.end(), beg2, end2);
            swap(*pFill, fill);

            result = true;
        }
    }
    return result;
}

template <unsigned K>
bool FragmentFiller<K>::checkLeadingMate( ELC& read, ELC& fill ) const
{
    ELItr begR(read.begin());
    ELItr endR(read.end());
    ELItr startR(endR);
    ELItr begF(fill.begin());
    ELItr endF(fill.end());
    ELItr startF(begF);
    ELItr stopF(begF);

    if ( begR == endR )
        return false;

    while ( startR != begR && stopF != endF )
    {
        if ( ++stopF == endF && fill.getFinalSkip() >= read.getFinalSkip() )
            break; // last line-up doesn't extend
        if ( std::equal(--startR, endR, begF) )
        {
            long endSkip = std::max(startR == begR ? read.getInitialSkip() : 0,
                    fill.getInitialSkip());
            endSkip += std::max(read.getFinalSkip(),
                    stopF == endF ? fill.getFinalSkip() : 0);
            if ( -separation(startR, endR, endSkip) >= mMinOverlap )
                return true;
        }
    }

    while ( stopF != endF )
    {
        if ( ++stopF == endF && fill.getFinalSkip() >= read.getFinalSkip() )
            break; // last line-up doesn't extend
        if ( std::equal(begR, endR, ++startF) )
        {
            long endSkip = read.getInitialSkip();
            endSkip += std::max(read.getFinalSkip(),
                    stopF == endF ? fill.getFinalSkip() : 0);
            if ( -separation(begR, endR, endSkip) >= mMinOverlap )
                return true;
        }
    }

    return false;
}

template <unsigned K>
bool FragmentFiller<K>::checkTrailingMate( ELC& read, ELC& fill ) const
{
    ELItr begR(read.begin());
    ELItr endR(read.end());
    ELItr stopR(begR);
    ELItr begF(fill.begin());
    ELItr endF(fill.end());
    ELItr startF(endF);
    ELItr stopF(endF);

    if ( begR == endR )
        return false;

    while ( stopR != endR && startF != begF )
    {
        if ( --startF == begF && fill.getInitialSkip() >= read.getInitialSkip() )
            break; // last line-up doesn't extend
        if ( std::equal(begR, ++stopR, startF) )
        {
            long endSkip = std::max(read.getInitialSkip(),
                    startF == begF ? fill.getInitialSkip() : 0);
            endSkip += std::max(stopR == endR ? read.getFinalSkip() : 0,
                    fill.getFinalSkip());
            if ( -separation(begR, stopR, endSkip) >= mMinOverlap )
                return true;
        }
    }

    while ( startF != begF )
    {
        if ( --startF == begF && fill.getInitialSkip() >= read.getInitialSkip() )
            break; // last line-up doesn't extend
        --stopF;
        if ( std::equal(begR, endR, startF) )
        {
            long endSkip = std::max(read.getInitialSkip(),
                    startF == begF ? fill.getInitialSkip() : 0);
            endSkip += read.getFinalSkip();
            if ( -separation(begR, endR, endSkip) >= mMinOverlap )
                return true;
        }
    }

    return false;
}

#endif // FRAGMENTFILLERDEFS_H_
