///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file HKPMerger.cc
 * \author tsharpe
 * \date Jan 29, 2010
 *
 * \brief
 */
#include "paths/HKPMerger.h"
#include "paths/InternalMerge.h"
#include "pairwise_aligners/PerfectAlignerLG.h"
#include "system/Assert.h"
#include "system/ThreadsafeIO.h"
#define cout ThreadsafeOStreamFetcher::cout()

HKPMerger::HKPMerger( const vec<HyperKmerPath> & HKPs,
        const vecbasevector & HKP_bases, const vec<int> & base_to_HKP_ID,
        const String sub_dir, int n_threads, const int max_kmer_freq,
        const int min_align_length, const int max_Q_size,
        NegativeGapValidator const* ngv, const vec<tagged_rpint>& uniqdb,
        long checkpointInterval, String const& checkpointFile )
: mMaxQSize(max_Q_size), mNGV(ngv), mUniqDB(uniqdb),
  mCheckpointInterval(checkpointInterval), mCheckpointFile(checkpointFile),
  mWorklist(mProcessor,n_threads)
{
    // HKP_bases contains the basevectors that correspond to the KmerPaths in
    // all the HyperKmerPaths in HKPs.
    // Find all K-mers (with K = min_align_length) that appear multiple times in
    // HKP_bases.  Then use these K-mers to create a graph of adjacencies
    // between the HyperKmerPaths.  Two HyperKmerPaths are marked as adjacent
    // (indicating they should be merged) if they share a non-repetitive K-mer.
    PerfectAlignerLG aligner(min_align_length, PerfectAlignerLG::findImproper);
    vec<alignment_plus> aligns;
    aligner.Align( HKP_bases, aligns, n_threads, -1, max_kmer_freq );

    size_t nNodes = HKPs.size();
    mNodes.reserve(nNodes);
    for ( size_t iii = 0; iii < nNodes; ++iii )
        mNodes.push_back(new Node(HKPs[iii],iii));

    // Fill the graph with adjacencies.  Note that adjacencies are
    // symmetric: iff the edge a->b exists, b->a also exists.
    // The graph represented by the nodes' adjacencies is an undirected graph.
    typedef vec<alignment_plus>::iterator AlignItr;
    for ( AlignItr itr(aligns.begin()), end(aligns.end()); itr != end; ++itr )
    {
        Node* pNode1 = mNodes[base_to_HKP_ID[itr->Id1()]];
        Node* pNode2 = mNodes[base_to_HKP_ID[itr->Id2()]];
        if ( pNode1 != pNode2 )
        {
            pNode1->addAdjacency(pNode2);
            pNode2->addAdjacency(pNode1);
        }
    }
}

HyperKmerPath HKPMerger::localMerge( int max_group_size,
                                      int min_overlap, int min_proper_overlap )
{
    ForceAssertGt(max_group_size,1);

    // if there's nothing to merge, just return an empty HyperKmerPath
    if ( !mNodes.size() )
        return HyperKmerPath();

    Comparator comparator;
    Checkpointer chkPtr;
    FinalMerger fMerger(chkPtr,mNodes.size());

    // while there are nodes that have adjacencies
    while ( mNodes.size() )
    {
        using std::sort;
        sort(mNodes.begin(),mNodes.end(),comparator);
        Node* pNode = mNodes.back();
        mNodes.resize(mNodes.size()-1);

        if ( !pNode->getAdjacencyCount() )
        {
            fMerger.addNode(pNode);
            continue;
        }

        LocalMerger* pLMerger = new LocalMerger(chkPtr,max_group_size,*this,
                                                min_overlap,min_proper_overlap);
        pLMerger->addNode(pNode);
        for ( int iii = 1; iii < max_group_size; ++iii )
        {
            std::set<Node*> const& xadjs = pLMerger->getAdjacencies();
            if ( xadjs.empty() )
                break;

            std::vector<Node*> adjs( xadjs.begin(), xadjs.end() );
            sort(adjs.begin(),adjs.end(),comparator);
            pNode = adjs.back();
            pLMerger->addNode(pNode);
            using std::find;
            mNodes.erase( find(mNodes.begin(),mNodes.end(),pNode) );
        }
        mNodes.push_back(pLMerger);
        pLMerger->schedule();
    }

    cout << Date() << ": Merge order is " << fMerger.getName() << std::endl;

    while ( !fMerger.wait(mCheckpointInterval) )
    {
        HyperKmerPath hkp = fMerger.checkpoint();
        HKPCleanup(hkp);
        BinaryOverwrite(mCheckpointFile,hkp);
    }

    return fMerger.getPath();
}

HyperKmerPath HKPMerger::Checkpointer::checkpoint()
{
    cout << Date() << ": Checkpointing " << mHKPs.size() << " HyperKmerPaths." << std::endl;
    typedef std::set<HyperKmerPath const*>::iterator Itr;
    vec<HyperKmerPath> hkps;
    hkps.reserve(mHKPs.size());
    for ( Itr itr(mHKPs.begin()), end(mHKPs.end()); itr != end; ++itr )
        hkps.push_back(**itr);
    return HyperKmerPath(hkps[0].K(),hkps);
}

HKPMerger::OutputAccumulator::~OutputAccumulator()
{
    delete &getPath();
    typedef vec<HyperKmerPath>::iterator Itr;
    for ( Itr itr(mToMerge.begin()), end(mToMerge.end()); itr != end; ++itr )
        mChkPtr.removeHKP(*itr);
}

void HKPMerger::OutputAccumulator::addNode( Node* pNode )
{
    Locker locker(mLock);

    String& name = getName();
    if ( !name.size() )
        name = '[' + pNode->getName() + ']';
    else
    {
        name.resize(name.size()-1);
        name += ',' + pNode->getName() + ']';
    }

    if ( pNode->isPathAvailable() )
        addMerge(pNode);
    else
    {
        mPending.insert(pNode);
        pNode->setListener(this);
    }
}

void HKPMerger::OutputAccumulator::internalPathAvailable( Node* pNode )
{
    mPending.erase(pNode);
    addMerge(pNode);
}

void HKPMerger::OutputAccumulator::addMerge( Node* pNode )
{
    mToMerge.push_back(pNode->getPath());
    Locker chkptLocker(mChkPtr);
    delete pNode;
    mChkPtr.addHKP(mToMerge.back());
}

void HKPMerger::LocalMerger::addNode( Node* pNode )
{
    removeAdjacency(pNode);
    pNode->removeAdjacency(this);

    typedef std::set<Node*>::const_iterator Itr;
    std::set<Node*> const& adjs = pNode->getAdjacencies();
    for ( Itr itr(adjs.begin()), end(adjs.end()); itr != end; ++itr )
    {
        Node* pNode2 = *itr;
        pNode2->removeAdjacency(pNode);
        pNode2->addAdjacency(this);
        addAdjacency(pNode2);
    }
    pNode->clearAdjacencies();
    OutputAccumulator::addNode(pNode);
}

void HKPMerger::LocalMerger::schedule()
{
    if ( !getNPending() )
    {
        LocalMerger* me = this;
        mHKPM.mWorklist.add(me);
    }
}

void HKPMerger::LocalMerger::doMerge()
{
    cout << Date() << ": Merging " << getName() << std::endl;
    vec<HyperKmerPath> const& merges = getMerges();
    HyperKmerPath* pHKP = new HyperKmerPath(merges[0].K(),merges);
    mHKPM.doInternalMerge(*pHKP,mMinOverlap,mMinProperOverlap);
    setPath(pHKP);
}

void HKPMerger::LocalMerger::pathAvailable( Node* pNode )
{
    Locker locker(getLock());
    internalPathAvailable(pNode);
    schedule();
}

HKPMerger::FinalMerger::FinalMerger( Checkpointer& chkPtr, size_t maxMerges )
: OutputAccumulator(chkPtr,maxMerges), mCondVar(getLock())
{
}

HKPMerger::FinalMerger::~FinalMerger()
{
}

bool HKPMerger::FinalMerger::wait( long nSecs )
{
    Locker locker(getLock());
    if ( !isPathAvailable() )
    {
        if ( getNPending() )
        {
            locker.timedWait(mCondVar,nSecs);
        }
        else
        {
            vec<HyperKmerPath> const& merges = getMerges();
            setPath(new HyperKmerPath(merges[0].K(),merges));
        }
    }

    return isPathAvailable();
}

void HKPMerger::FinalMerger::pathAvailable( Node* pNode )
{
    Locker locker(getLock());
    internalPathAvailable(pNode);
    if ( !getNPending() )
        mCondVar.signal();
}
