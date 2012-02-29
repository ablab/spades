///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ReadPatherDefs.h
 * \author tsharpe
 * \date May 4, 2010
 *
 * \brief
 */
#ifndef KMERS_READ_PATHER_DEFS_H_
#define KMERS_READ_PATHER_DEFS_H_

#ifndef NDEBUG
#define LOG(x) std::cout << x << std::endl
#else
#define LOG(x) ((void)0)
#endif

#include "feudal/ChunkDumper.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/VirtualMasterVec.h"
#include "kmers/ReadPather.h"
#include "system/SpinLockedData.h"
#include "system/SysConf.h"
#include "system/System.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <time.h>
#include <unistd.h>

namespace
{
// worklist processor class to kmerize reads
template <unsigned K>
class KmerizationProcessor
{
public:
    KmerizationProcessor( VirtualMasterVec<bvec> const& vmv, KmerDict<K>& dict,
                            Dotter& dotter )
    : mVMV(vmv), mDict(dict), mDotter(dotter)
    {}

    KmerizationProcessor( KmerizationProcessor const& that )
    : mVMV(that.mVMV), mDict(that.mDict), mDotter(that.mDotter)
    { mReads.reserve(BATCHSIZE); }

    // copying prohibited by references, compiler-supplied destructor is OK

    void operator()( size_t );

    static size_t getBatchSize() { return BATCHSIZE; }
    static size_t getNBatches( size_t nReads )
    { return (nReads + BATCHSIZE - 1)/BATCHSIZE; }

private:
    static size_t const BATCHSIZE = 10000;

    VirtualMasterVec<bvec> mVMV;
    KmerDict<K>& mDict;
    vecbvec mReads;
    Dotter& mDotter;
};

// output iterator that just increments counts for canonicalized kmers
template <unsigned K>
class DictionaryOutputIterator
{
public:
    DictionaryOutputIterator( KmerDict<K>* pDict )
    : mpDict(pDict) {}

    // compiler-supplied copying and destructor are OK

    DictionaryOutputIterator& operator++() { return *this; }
    DictionaryOutputIterator& operator++( int ) { return *this; }
    DictionaryOutputIterator& operator*() { return *this; }
    KMer<K> const& operator=( KMer<K> const& kmer )
    { mpDict->refCanonical(kmer).incrCount(); return kmer; }

private:
    KmerDict<K>* mpDict;
};

template <unsigned K>
void KmerizationProcessor<K>::operator()( size_t batchNo )
{
    bvec bv;
    size_t readId = batchNo*BATCHSIZE;
    size_t endId = std::min(readId+BATCHSIZE,mVMV.size());
    DictionaryOutputIterator<K> oItr(&mDict);
    while ( readId != endId )
    {
        bv = mVMV[readId++];
        KMer<K>::kmerize(bv.begin(),bv.end(),oItr);
    }
    mDotter.batchDone();
}

// creates the files that represent a unipath graph
template <unsigned K>
class GraphBuilder
{
public:
    GraphBuilder( String const& graphInfoFilename, KmerDict<K>& dict )
    : mGraphInfoFilename(graphInfoFilename), mDict(dict), mStepper(dict),
      mSeq(UnipathGraph<K>::getSeqFilename(graphInfoFilename).c_str()),
      mNComponents(0)
    { unlink(graphInfoFilename.c_str()); }

    // compiler-supplied destructor is OK

    void process();

private:
    GraphBuilder( GraphBuilder const& ); // unimplemented -- no copying
    GraphBuilder& operator=( GraphBuilder const& ); // unimplemented--no copying

    typedef typename KmerDict<K>::Entry DictEntry;
    struct EdgeStart
    {
        EdgeStart( DictEntry const& entry, bool isRC, bool isExtension )
        : mpEntry(&entry), mRC(isRC), mExtension(isExtension) {}

        DictEntry const* mpEntry;
        bool mRC;
        bool mExtension;
    };
    typedef std::vector<EdgeStart> EdgeStartVec;

    void buildEdges( EdgeStartVec& toBuild );
    void joinGraph( HugeBVec const& bv );

    EdgeID fromRef( UnipathEdge const& unipath ) const
    { return EdgeID(&unipath - &mEdges[0]); }

    UnipathEdge& toRef( EdgeID const& edgeID )
    { AssertNot(edgeID.isNull()); return mEdges[edgeID.val()]; }

    static HugeBVec::const_iterator getBases( HugeBVec const& bv,
                                              KmerID const& kmerID )
    { return bv.begin(kmerID.val()); }

    void joinSuccessor( KMer<K> const& kmer, KDef* pDef, unsigned predCode,
                            EdgeID const& predID, HugeBVec const& bv )
    { UnipathEdge& pred = toRef(predID);
      EdgeID const& succID = pDef->getEdgeID();
      UnipathEdge& succ = toRef(succID);
      KmerID succKmerID = succ.getKmerID(pDef->getEdgeOffset());
      unsigned succCode = kmer.back();
      using std::equal;
      if ( succKmerID == succ.getInitialKmerID() &&
           equal(kmer.begin(),kmer.end(),getBases(bv,succ.getInitialKmerID())) )
      { pred.setSuccessor(succCode,succID,false);
        succ.setPredecessor(predCode,predID,false); }
      else
      { AssertEq( succKmerID, succ.getFinalKmerID() );
        Assert( equal(kmer.rcbegin(),kmer.rcend(),
                      getBases(bv,succ.getFinalKmerID())) );
        pred.setSuccessor(succCode,succID,true);
        succ.setSuccessor(GetComplementaryBase(predCode),predID,true); } }

    void joinPredecessor( unsigned succCode, EdgeID const& succID,
                            KMer<K> const& kmer, KDef* pDef,
                            HugeBVec const& bv )
    { UnipathEdge& succ = toRef(succID);
      EdgeID const& predID = pDef->getEdgeID();
      UnipathEdge& pred = toRef(predID);
      KmerID predKmerID = pred.getKmerID(pDef->getEdgeOffset());
      unsigned predCode = kmer.front();
      using std::equal;
      if ( predKmerID == pred.getFinalKmerID() &&
           equal(kmer.begin(),kmer.end(),getBases(bv,pred.getFinalKmerID())) )
      { pred.setSuccessor(succCode,succID,false);
        succ.setPredecessor(predCode,predID,false); }
      else
      { AssertEq( predKmerID, pred.getInitialKmerID() );
        Assert( equal(kmer.rcbegin(),kmer.rcend(),
                      getBases(bv,pred.getInitialKmerID())) );
        pred.setPredecessor(GetComplementaryBase(succCode),succID,true);
        succ.setPredecessor(predCode,predID,true); } }

    // median edge length in bases
    size_t medianEdgeLen() const
    { std::vector<size_t> lengths;
      size_t nnn = mEdges.size();
      lengths.reserve(nnn);
      typedef UnipathEdgeVec::const_iterator Itr;
      for ( Itr itr(mEdges.begin()), end(mEdges.end()); itr != end; ++itr )
        lengths.push_back(itr->getLength());
      std::vector<size_t>::iterator med(lengths.begin()+nnn/2);
      std::nth_element(lengths.begin(),med,lengths.end());
      return *med + K - 1; }

    String mGraphInfoFilename;
    KmerDict<K>& mDict;
    KmerStepper<K> mStepper;
    HugeBVec::Builder mSeq;
    UnipathEdgeVec mEdges;
    size_t mNComponents;
};

template <unsigned K>
void GraphBuilder<K>::process()
{
    std::cout << Date() << " Building unipaths." << std::endl;

    mEdges.reserve(mDict.size()/100);

    typedef typename KmerDict<K>::OCItr OCItr;
    typedef typename KmerDict<K>::ICItr ICItr;
    OCItr oEnd(mDict.cend());
    EdgeStartVec toBuild;

    DictEntry const* entries[4];

    for ( OCItr oItr(mDict.cbegin()); oItr != oEnd; ++oItr )
    {
        for ( ICItr itr(oItr->begin()), end(oItr->end()); itr != end; ++itr )
        {
            DictEntry const& entry = *itr;
            if ( entry.getKDef().isNull() )
            {
                if ( mStepper.getPredecessors(entry,entries) != 1 )
                {
                    LOG("IF: " << static_cast<KMer<K> const&>(entry));
                    toBuild.push_back(EdgeStart(entry,false,false));
                    buildEdges(toBuild);
                    mNComponents += 1;
                }
                else if ( mStepper.getSuccessors(entry,entries) != 1 )
                {
                    LOG("IR: " << static_cast<KMer<K> const&>(entry));
                    toBuild.push_back(EdgeStart(entry,true,false));
                    buildEdges(toBuild);
                    mNComponents += 1;
                }
            }
        }
    }

    // leaves only smooth rings, but we'd better scan for that.  you never know.
    for ( OCItr oItr(mDict.cbegin()); oItr != oEnd; ++oItr )
    {
        for ( ICItr itr(oItr->begin()), end(oItr->end()); itr != end; ++itr )
        {
            DictEntry const& entry = *itr;
            if ( entry.getKDef().isNull() )
            {
                LOG("I1: " << static_cast<KMer<K> const&>(entry));
                toBuild.push_back(EdgeStart(entry,false,false));
                buildEdges(toBuild);
                mNComponents += 1;
            }
        }
    }

    mSeq.close();
    HugeBVec seq(UnipathGraph<K>::getSeqFilename(mGraphInfoFilename).c_str());
    joinGraph(seq);

    String dictFile = UnipathGraph<K>::getDictFilename(mGraphInfoFilename);
    BinaryWriter::writeFile(dictFile.c_str(),mDict);
    String edgesFile = UnipathGraph<K>::getEdgesFilename(mGraphInfoFilename);
    BinaryWriter::writeFile(edgesFile.c_str(),mEdges);
    GraphInfo gi;
    gi.K = K;
    gi.kSeqSize = mSeq.size();
    gi.nComponents = mNComponents;
    gi.nKmers = mDict.size();
    gi.nEdges = mEdges.size();
    BinaryWriter writer(mGraphInfoFilename.c_str());
    writer.write(gi);
    writer.close();
    std::cout << "There are " << mEdges.size() << " unipaths in "
              << mNComponents << " connected components.\n"
              "The median unipath length is " << medianEdgeLen() << " bases."
              << std::endl;
}

template <unsigned K>
void GraphBuilder<K>::buildEdges( EdgeStartVec& toBuild )
{
    KMer<K> kmer;

    while ( !toBuild.empty() )
    {
        EdgeStart const& es = toBuild.back();
        KDef* pDef = const_cast<KDef*>(&es.mpEntry->getKDef());
        if ( !pDef->isNull() )
        {
            toBuild.pop_back();
            continue;
        }

        kmer = *es.mpEntry;
        if ( es.mRC ) kmer.rc();
        if ( es.mExtension )
            mSeq.push_back(kmer.back());
        else
            mSeq.append(kmer.begin(),kmer.end());

        EdgeID edgeID(mEdges.size());
        bool isPalindrome = kmer.isPalindrome();
        mEdges.push_back(UnipathEdge(KmerID(mSeq.size()-kmer.size()),
                                     ComponentID(mNComponents),
                                     isPalindrome));
        pDef->set( edgeID, EdgeOffset(0) );
        toBuild.pop_back();

        LOG("BF: " << kmer << ' ' << edgeID);

        DictEntry const* entries[4];

        // queue predecessors
        if ( mStepper.getPredecessors(kmer,entries) )
        {
            for ( unsigned predCode = 0; predCode < 4u; ++predCode )
            {
                DictEntry const* pEntry = entries[predCode];
                if ( pEntry && pEntry->getKDef().isNull() )
                {
                    KMer<K>& kmerPred = mStepper.getSteppedKmer();
                    kmerPred.setFront(predCode);
                    LOG("QR: " << KMer<K>(kmerPred).rc() << ' ' << edgeID);
                    bool isRC = kmerPred.getCanonicalForm() != REV;
                    toBuild.push_back(EdgeStart(*pEntry,isRC,false));
                }
            }
        }

        // while we're on the straight and narrow, with just a single successor,
        // try to extend the edge
        unsigned nSuccessors;
        while ( (nSuccessors = mStepper.getSuccessors(kmer,entries)) == 1 )
        {
            // if we're a palindrome we must queue successors
            if ( isPalindrome )
                break;
            nSuccessors = 0; // we'll handle our single successor now

            DictEntry const* pEntry;
            for ( unsigned succCode = 0; succCode < 4u; ++succCode )
            {
                if ( (pEntry = entries[succCode]) )
                {
                    kmer = mStepper.getSteppedKmer();
                    kmer.setBack(succCode);
                    break;
                }
            }

            // if single successor is already built, we're done extending
            if ( !pEntry->getKDef().isNull() )
                break;

            // if successor has multiple predecessors, or is a palindrome
            // start a new edge
            CanonicalForm form = kmer.getCanonicalForm();
            if ( form==PALINDROME ||
                    mStepper.getPredecessors(kmer,entries) > 1 )
            {
                LOG("QF: " << kmer << ' ' << edgeID);
                toBuild.push_back(EdgeStart(*pEntry,form==REV,true));
                break;
            }

            // no branching -- extend the edge we're building
            LOG("XF: " << kmer << ' ' << edgeID);
            mSeq.push_back(kmer.back());
            pDef = const_cast<KDef*>(&pEntry->getKDef());
            pDef->set( edgeID, EdgeOffset(toRef(edgeID).extend()) );
        }

        // if there are unhandled successors (multiple successors, or the single
        // successor of a palindrome), queue them up
        if ( nSuccessors )
        {
            for ( unsigned succCode = 0; succCode < 4u; ++succCode )
            {
                DictEntry const* pEntry = entries[succCode];
                if ( pEntry && pEntry->getKDef().isNull() )
                {
                    KMer<K>& kmerSucc = mStepper.getSteppedKmer();
                    kmerSucc.setBack(succCode);
                    LOG("QF: " << kmerSucc << ' ' << edgeID);
                    bool isRC = kmerSucc.getCanonicalForm() == REV;
                    toBuild.push_back(EdgeStart(*pEntry,isRC,!--nSuccessors));
                }
            }
        }
    }
}

template <unsigned K>
void GraphBuilder<K>::joinGraph( HugeBVec const& bv )
{
    typedef UnipathEdgeVec::iterator Itr;
    DictEntry const* entries[4];
    KMer<K> kmer;
    for ( size_t idx = 0; idx < mEdges.size(); ++idx )
    {
        UnipathEdge& edge = mEdges[idx];
        kmer.assign(getBases(bv,edge.getInitialKmerID()));
        unsigned succCode = kmer.back();
        unsigned predCode;
        if ( mStepper.getPredecessors(kmer,entries) )
        {
            EdgeID succID = fromRef(edge);
            for ( predCode = 0; predCode < 4u; ++predCode )
            {
                DictEntry const* pEntry = entries[predCode];
                if ( pEntry && edge.getPredecessor(predCode).isNull() )
                {
                    KDef* pDef = const_cast<KDef*>(&pEntry->getKDef());
                    KMer<K>& kmerPred = mStepper.getSteppedKmer();
                    kmerPred.setFront(predCode);
                    LOG("JP: " << kmerPred << ' ' << succID <<
                            " onto " << pDef->getEdgeID());
                    joinPredecessor(succCode,succID,kmerPred,pDef,bv);
                }
            }
        }
        kmer.assign(getBases(bv,edge.getFinalKmerID()));
        predCode = kmer.front();
        if ( mStepper.getSuccessors(kmer,entries) )
        {
            EdgeID predID = fromRef(edge);
            for ( succCode = 0; succCode < 4u; ++succCode )
            {
                DictEntry const* pEntry = entries[succCode];
                if ( pEntry && edge.getSuccessor(succCode).isNull() )
                {
                    KDef* pDef = const_cast<KDef*>(&pEntry->getKDef());
                    KMer<K>& kmerSucc = mStepper.getSteppedKmer();
                    kmerSucc.setBack(succCode);
                    LOG("JS: " << kmerSucc << ' ' << predID <<
                            " onto " << pDef->getEdgeID());
                    joinSuccessor(kmerSucc,pDef,predCode,predID,bv);
                }
            }
        }
    }
}

template <unsigned K> class PathingProcessor;

// Paths sequence onto UnipathGraph.
template <unsigned K>
class PathBuilder
{
public:
    PathBuilder( String const& infoFile, UnipathGraph<K> const& graph );

    void processReads( String const& fastb, KmerDict<K> const& dict,
                        unsigned nThreads );

    void writeEvidence( unsigned nThreads );

    // compiler-supplied destructor is OK

private:
    PathBuilder( PathBuilder const& ); // unimplemented -- no copying
    PathBuilder& operator=( PathBuilder const& ); // unimplemented -- no copying

    // the rest of these methods are for use by the friendly PathingProcessor
    friend class PathingProcessor<K>;

    bool isRC( KMer<K> const& kmer, KmerID const& kmerID ) const
    { using std::equal;
      HugeBVec::const_iterator seqItr(mGraph.getBases(kmerID));
      bool result = !equal(kmer.begin(),kmer.end(),seqItr);
      Assert(!result || equal(kmer.rcbegin(),kmer.rcend(),seqItr));
      return result; }

    UnipathEdge const& getEdge( EdgeID const& edgeID ) const
    { return mGraph.getEdge(edgeID); }

    bool nextEdge( unsigned baseCode, EdgeID* pEdgeID, bool isRC ) const
    { return mGraph.nextEdge(baseCode,pEdgeID,isRC); }

    PathID addPathSeq( bvec const& bases )
    { SpinLocker locker(mPathSeqLock);
      PathID pathID(mPathSeq.size());
      mPathSeq.append(bases.begin(),bases.end());
      return pathID; }

    String mInfoFile;
    UnipathGraph<K> const& mGraph;
    HugeBVec::Builder mPathSeq; // indexed by PathID
    SpinLockedData mPathSeqLock;
};

// worklist processor class to kmerize reads
template <unsigned K>
class PathingProcessor
{
public:
    PathingProcessor( VirtualMasterVec<bvec> const& vmv,
                      KmerDict<K> const& dict,
                      PathBuilder<K>& pathBuilder,
                      ChunkDumper<UnipathingVec>& pathingsDumper,
                      Dotter& dotter )
    : mVMV(vmv), mDict(dict), mPathBuilder(pathBuilder),
      mPathingsDumper(pathingsDumper), mDotter(dotter), mTotalReadLen(0)
    {}

    PathingProcessor( PathingProcessor const& that )
    : mVMV(that.mVMV), mDict(that.mDict), mPathBuilder(that.mPathBuilder),
      mPathingsDumper(that.mPathingsDumper), mDotter(that.mDotter),
      mTotalReadLen(0)
    { mReads.reserve(BATCHSIZE); }

    ~PathingProcessor()
    { SpinLocker locker(gLock); gTotalReadLen += mTotalReadLen; }

    // copy assignment prohibited by references, compiler-supplied destructor OK

    static size_t getBatchSize() { return BATCHSIZE; }
    static size_t getNBatches( size_t nReads )
    { return (nReads+BATCHSIZE-1)/BATCHSIZE; }

    void operator()( size_t batchNo )
    { size_t readId = batchNo*BATCHSIZE;
      size_t endId = std::min(readId+BATCHSIZE,mVMV.size());
      size_t nnn = endId - readId;
      mPathings.reserve(nnn);
      for ( ;readId != endId; ++readId )
          path(readId,mVMV[readId]);
      mPathingsDumper.dumpChunk(batchNo,mPathings);
      mDotter.batchDone(); }

    static void clearTotalReadLen() { gTotalReadLen = 0; }
    static size_t getTotalReadLen() { return gTotalReadLen; }

private:
    static size_t const BATCHSIZE = 10000;

    void path( size_t readId, bvec const& read );

    VirtualMasterVec<bvec> mVMV;
    KmerDict<K> const& mDict;
    PathBuilder<K>& mPathBuilder;
    ChunkDumper<UnipathingVec>& mPathingsDumper;
    Dotter& mDotter;
    vecbvec mReads;
    bvec mPathSeq;
    UnipathingVec mPathings;
    size_t mTotalReadLen;

    static size_t gTotalReadLen; // sum of read lengths in kmers
    static SpinLockedData gLock;
};

template <unsigned K> SpinLockedData PathingProcessor<K>::gLock;
template <unsigned K> size_t PathingProcessor<K>::gTotalReadLen;

template <unsigned K>
void PathingProcessor<K>::path( size_t readId, bvec const& read )
{
    size_t readLen = read.size();
    if ( readLen < K )
    {
        mPathings.push_back(Unipathing());
        return; // LOOKOUT:  early return
    }

    mTotalReadLen += readLen - K + 1;
    KMer<K> kmer(read.begin());
    KDef const* pDef = mDict.lookup(kmer);
    if ( !pDef )
        FatalErr("Can't find kmer " << kmer << " in dictionary.");

    UnipathEdge const* pEdge = &mPathBuilder.getEdge(pDef->getEdgeID());
    KmerID kmerID = pEdge->getKmerID(pDef->getEdgeOffset());
    bool isRC = mPathBuilder.isRC(kmer,kmerID);
    size_t skipped = isRC ? pEdge->getFinalKmerID().val() - kmerID.val() :
                            kmerID.val() - pEdge->getInitialKmerID().val();
    Unipathing up(pDef->getEdgeID(),isRC,skipped);
    size_t readIdx = K + pEdge->getLength() - skipped - 1;

    if ( readIdx < readLen )
    {
        mPathSeq.clear();
        EdgeID edgeID = pDef->getEdgeID();
        do
        {
            unsigned char baseCode = read[readIdx];
            if ( isRC ) baseCode = GetComplementaryBase(baseCode);
            mPathSeq.push_back(baseCode);
            isRC = mPathBuilder.nextEdge(baseCode,&edgeID,isRC);
            if ( edgeID.isNull() )
                FatalErr("Can't path read " << readId << " " << read.ToString()
                         << " at " << readIdx );
            pEdge = &mPathBuilder.getEdge(edgeID);
            readIdx += pEdge->getLength();
        }
        while ( readIdx < readLen );

        unsigned nBases = mPathSeq.size();
        up.setSegments( nBases+1,
                        nBases <= 8 ?
                           PathID(8*mPathSeq.extractKmer(0,nBases)) :
                           mPathBuilder.addPathSeq(mPathSeq) );
    }

    up.setFinalSkip(readIdx-readLen);
    mPathings.push_back(up);
}

template <unsigned K>
PathBuilder<K>::PathBuilder( String const& infoFile,
                                UnipathGraph<K> const& graph )
: mInfoFile(infoFile), mGraph(graph),
  mPathSeq(PathCollection<K>::getPathseqFilename(infoFile).c_str())
{
    unlink(infoFile.c_str());

    // seed traversal sequence with all 8-base sequences
    bvec bv;
    for ( unsigned idx = 0; idx < 65536; ++idx )
    {
        bv.assignBaseBits(8,&idx);
        mPathSeq.append(bv.begin(),bv.end());
    }
    PathingProcessor<K>::clearTotalReadLen();
}

template <unsigned K>
void PathBuilder<K>::processReads( String const& fastb, KmerDict<K> const& dict,
                                    unsigned nThreads )
{
    std::cout << Date() << " Pathing reads from " << fastb << std::endl;
    VirtualMasterVec<bvec> vmv(fastb.c_str());
    size_t nReads = vmv.size();
    size_t nBatches = PathingProcessor<K>::getNBatches(nReads);
    Dotter dotter(nBatches);
    String pathsFilename = PathCollection<K>::getPathsFilename(mInfoFile);
    ChunkDumper<UnipathingVec> pathingsDumper(pathsFilename.c_str(),
                                              nReads,nBatches);
    PathingProcessor<K> proc(vmv,dict,*this,pathingsDumper,dotter);
    std::cout << "Pathing " << nReads << " reads in "
                << nBatches << " batches of "
                << PathingProcessor<K>::getBatchSize() << '.' << std::endl;

    if ( !nThreads )
        nThreads = processorsOnline();
    WorklistN<PathingProcessor<K> > wl(proc,nBatches,nThreads-1);
}

template <unsigned K>
class EvidenceProcessor
{
public:
    struct Common : private Dotter
    {
        Common( UnipathGraph<K> const& graph,
                      HugeBVec const& pathSeq,
                      UnipathingVec const& pathings,
                      ChunkDumper<VecUnipathEvidenceVec>& dumper,
                      size_t nBatches, size_t batchSize,
                      std::vector<unsigned> const& evCounts )
        : Dotter(nBatches), mGraph(graph), mPathSeq(pathSeq), mPathings(pathings),
          mDumper(dumper), mBatchSize(batchSize), mEvCounts(evCounts)
        {}

        void dumpChunk( size_t chunk, VecUnipathEvidenceVec& vuev )
        { mDumper.dumpChunk(chunk,vuev); batchDone(); }

        UnipathGraph<K> const& mGraph;
        HugeBVec const& mPathSeq;
        UnipathingVec const& mPathings;
        ChunkDumper<VecUnipathEvidenceVec>& mDumper;
        size_t mBatchSize;
        std::vector<unsigned> const& mEvCounts;
    };

    EvidenceProcessor( Common& common )
    : mCommon(common)
    {}

    EvidenceProcessor( EvidenceProcessor const& that )
    : mCommon(that.mCommon)
    {}

    // copy assignment prohibited by references, compiler-supplied destructor OK

    void operator()( size_t chunk );

private:
    Common& mCommon;
    VecUnipathEvidenceVec mVUEV;
};

template <unsigned K>
void EvidenceProcessor<K>::operator()( size_t chunk )
{
    typedef VecUnipathEvidenceVec::iterator Itr;
    size_t minEdgeId = chunk*mCommon.mBatchSize;
    size_t maxEdgeId = std::min(minEdgeId+mCommon.mBatchSize,
                                mCommon.mGraph.getNEdges());
    mVUEV.clear().resize(maxEdgeId-minEdgeId);
    size_t idx = minEdgeId;
    for ( Itr itr(mVUEV.begin()), end(mVUEV.end()); itr != end; ++itr )
        itr->reserve(mCommon.mEvCounts[idx++]);

    size_t nReads = mCommon.mPathings.size();
    for ( size_t readId = 0; readId < nReads; ++readId )
    {
        ReadID rid(readId);
        Unipathing const& pathing = mCommon.mPathings[readId];
        EdgeID edgeID = pathing.getInitialEdgeID();
        bool isRC = pathing.isInitialEdgeRC();
        size_t edgeId = edgeID.val();
        if ( minEdgeId <= edgeId && edgeId < maxEdgeId )
        {
            UnipathEvidenceVec& ev = mVUEV[edgeId-minEdgeId];
            ev.push_back(UnipathEvidence(rid,0u,isRC));
        }

        unsigned nSegs = pathing.getNSegments();
        if ( nSegs > 1 )
        {
            size_t off = pathing.getPathID().val();
            HugeBVec::const_iterator itr = mCommon.mPathSeq.begin(off-1);
            for ( unsigned segNo = 1; segNo < nSegs; ++segNo )
            {
                isRC = mCommon.mGraph.nextEdge(*++itr,&edgeID,isRC);
                if ( edgeID.isNull() )
                    FatalErr("Invalid successor base in path sequence for read "
                             << readId );
                edgeId = edgeID.val();
                if ( minEdgeId <= edgeId && edgeId < maxEdgeId )
                {
                    UnipathEvidenceVec& ev = mVUEV[edgeId-minEdgeId];
                    ev.push_back(UnipathEvidence(rid,segNo,isRC));
                }
            }
        }
    }
    mCommon.dumpChunk(chunk,mVUEV);
}

template <unsigned K>
void PathBuilder<K>::writeEvidence( unsigned nThreads )
{
    std::cout << Date() << " Creating edge location dictionary." << std::endl;

    mGraph.ungetDict();

    UnipathingVec pathings;
    String pathsFilename = PathCollection<K>::getPathsFilename(mInfoFile);
    BinaryReader::readFile(pathsFilename.c_str(),&pathings);

    mPathSeq.close();
    String pathSeqFilename = PathCollection<K>::getPathseqFilename(mInfoFile);
    HugeBVec pathSeq(pathSeqFilename.c_str());

    size_t nEdges = mGraph.getNEdges();
    std::vector<unsigned> evCounts;
    evCounts.resize(nEdges);

    size_t nReads = pathings.size();
    size_t evTot = 0;
    for ( size_t readId = 0; readId < nReads; ++readId )
    {
        Unipathing const& pathing = pathings[readId];
        if ( pathing.isNull() )
            continue; // Non-structured code!

        EdgeID edgeID = pathing.getInitialEdgeID();
        evCounts[edgeID.val()] += 1;
        evTot += 1;
        unsigned nSegs = pathing.getNSegments();
        if ( nSegs > 1 )
        {
            bool isRC = pathing.isInitialEdgeRC();
            size_t off = pathing.getPathID().val();
            HugeBVec::const_iterator itr = pathSeq.begin(off-1);
            for ( unsigned segNo = 1; segNo < nSegs; ++segNo )
            {
                isRC = mGraph.nextEdge(*++itr,&edgeID,isRC);
                evCounts[edgeID.val()] += 1;
                evTot += 1;
            }
        }
    }

    if ( true )
    {
        size_t bytesPer = (evTot*sizeof(UnipathEvidence)+nEdges-1)/nEdges +
                                    sizeof(UnipathEvidenceVec);
        WorklistParameterizer wp(nEdges,bytesPer,100,nThreads);
        size_t nBatches = wp.getNBatches();
        size_t batchSize = wp.getBatchSize();
        String evFilename = PathCollection<K>::getEvidenceFilename(mInfoFile);
        ChunkDumper<VecUnipathEvidenceVec> dumper(evFilename.c_str(),
                                                  nEdges,nBatches);
        typename EvidenceProcessor<K>::Common common(mGraph, pathSeq, pathings,
                                                     dumper, nBatches,
                                                     batchSize, evCounts);
        EvidenceProcessor<K> proc(common);
        std::cout << "Recording evidence for " << nEdges << " unipaths in "
                    << nBatches << " batches of "
                    << batchSize << '.' << std::endl;
        WorklistN<EvidenceProcessor<K> > wl(proc,nBatches,wp.getNThreads()-1);
    }

    PathInfo pi;
    pi.nReads = pathings.size();
    pi.pathSeqSize = pathSeq.size();
    BinaryWriter::writeFile(mInfoFile.c_str(), pi);
}

} // end of anonymous namespace


template <unsigned K>
void KmerDict<K>::process( String const& fastb, unsigned nThreads )
{
    std::cout << Date() << " Kmerizing reads from " << fastb << std::endl;
    VirtualMasterVec<bvec> vmv(fastb.c_str());
    size_t nReads = vmv.size();
    size_t nBatches = KmerizationProcessor<K>::getNBatches(nReads);
    Dotter dotter(nBatches);
    KmerizationProcessor<K> proc(vmv,*this,dotter);
    std::cout << "Kmerizing " << nReads << " reads in "
                << nBatches << " batches of "
                << KmerizationProcessor<K>::getBatchSize() << '.' << std::endl;
    if ( !nThreads )
        nThreads = processorsOnline();
    WorklistN<KmerizationProcessor<K> > wl(proc,nBatches,nThreads-1);
}

template <unsigned K>
void UnipathGraph<K>::create( String const& graphInfoFilename,
                                String const& fastb,
                                unsigned nThreads,
                                size_t nKmersEstimate )
{
    String countsFile = UnipathGraph::getCountsFilename(graphInfoFilename);
    bool countsFileExists = IsRegularFile(countsFile);
    KmerDict<K> dict(countsFileExists?0:nKmersEstimate);

    if ( countsFileExists )
        BinaryReader::readFile(countsFile.c_str(),&dict);
    else
    {
        dict.process(fastb,nThreads);
        BinaryWriter::writeFile(countsFile.c_str(),dict);
    }

    create(graphInfoFilename,dict);
}

template <unsigned K>
void UnipathGraph<K>::create( String const& graphInfoFilename,
                                KmerDict<K>& dict )
{
    GraphBuilder<K> builder(graphInfoFilename,dict);
    builder.process();
}

template <unsigned K>
void PathCollection<K>::create( String const& fastb, bool validate,
                             unsigned nThreads, size_t nKmersEst,
                             bool writeKmerCounts )
{
    String graphInfoFilename = UnipathGraph<K>::getInfoFilename(fastb);
    String countsFile = UnipathGraph<K>::getCountsFilename(graphInfoFilename);
    KmerDict<K>* pDict;
    if ( IsRegularFile(countsFile) )
    {
        std::cout << Date() << " Reading kmer counts." << std::endl;
        pDict = new KmerDict<K>(0);
        BinaryReader::readFile(countsFile.c_str(),pDict);
    }
    else
    {
        if ( !nKmersEst )
            // assume 50x coverage
            nKmersEst = 4*MasterVec<bvec>::MastervecFileRawCount(fastb)/50;
        std::cout << Date() << " Allocating kmer dictionary of size "
                    << nKmersEst << '.' << std::endl;
        pDict = new KmerDict<K>(nKmersEst);
        pDict->process(fastb,nThreads);
        std::cout << Date() << " There were " << pDict->size() << " kmers."
                    << std::endl;
        if ( writeKmerCounts )
            BinaryWriter::writeFile(countsFile.c_str(),*pDict);
    }

    UnipathGraph<K>::create(graphInfoFilename,*pDict);
    UnipathGraph<K> graph(graphInfoFilename);
    if ( validate )
        graph.validateUnipaths(graph.getAllBases(),graph.getAllEdges(),*pDict);

    PathBuilder<K> pb(getInfoFilename(fastb),graph);
    pb.processReads(fastb,*pDict,nThreads);

    delete pDict;

    pb.writeEvidence(nThreads);
}

template <unsigned K>
void PathCollection<K>::create( String const& fastbFilename,
                               String const& graphInfoFilename,
                               bool validate, unsigned nThreads )
{
    UnipathGraph<K> graph(graphInfoFilename);
    if ( validate )
        graph.validate();
    String pathInfoFile = getInfoFilename(fastbFilename);
    PathBuilder<K> pb(pathInfoFile,graph);
    pb.processReads(fastbFilename,graph.getDict(),nThreads);
    pb.writeEvidence(nThreads);
}

#endif // KMERS_READ_PATHER_DEFS_H_
