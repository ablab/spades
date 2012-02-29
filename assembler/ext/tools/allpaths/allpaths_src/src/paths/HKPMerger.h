///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file HKPMerger.h
 * \author tsharpe
 * \date Jan 29, 2010
 *
 * \brief
 */
#ifndef HKPMERGER_H_
#define HKPMERGER_H_

#include "Basevector.h"
#include "String.h"
#include "Vec.h"
#include "paths/HyperKmerPath.h"
#include "paths/InternalMergeImpl.h"
#include "paths/KmerPathInterval.h"
#include "paths/NegativeGapValidator.h"
#include "system/Assert.h"
#include "system/LockedData.h"
#include "system/Worklist.h"
#include <set>
#include <vector>

class HKPMerger
{
public:
    HKPMerger( const vec<HyperKmerPath> & HKPs,
            const vecbasevector & HKP_bases, const vec<int> & base_to_HKP_ID,
            const String sub_dir, int n_threads, const int max_kmer_freq,
            const int min_align_length, const int max_Q_size,
            NegativeGapValidator const* ngv, const vec<tagged_rpint>& uniqdb,
            long checkpointInterval, String const& checkpointFile );

    HyperKmerPath localMerge( int max_group_size,
                                int min_overlap, int min_proper_overlap );

    void finalMerge( int min_proper_overlap_final, HyperKmerPath& hkp )
    { doInternalMerge(hkp,min_proper_overlap_final,min_proper_overlap_final); }

private:
    HKPMerger( HKPMerger const& ); // unimplemented -- no copying
    HKPMerger& operator=( HKPMerger const& ); // unimplemented -- no copying

    void doInternalMerge( HyperKmerPath& hkp, int min_overlap, int min_proper )
    { InternalMergeImpl(hkp,mNGV,min_overlap,min_proper,
                        mMaxQSize,True,True,mUniqDB); }

    class Node
    {
    public:
        Node() : mpHKP(0), mpListener(0) {}
        Node( HyperKmerPath const& hkp, unsigned id )
        : mpHKP(&hkp), mpListener(0), mName(ToString(id)) {}
        virtual ~Node() {}

        bool isPathAvailable() const { return mpHKP; }
        HyperKmerPath const& getPath() const { return *mpHKP; }

        void addAdjacency( Node* pNode ) { mAdjacencies.insert(pNode); }
        void removeAdjacency( Node* pNode ) { mAdjacencies.erase(pNode); }
        void clearAdjacencies() { mAdjacencies.clear(); }
        std::set<Node*> const& getAdjacencies() const { return mAdjacencies; }
        size_t getAdjacencyCount() const { return mAdjacencies.size(); }

        void setListener( Node* pNode )
        { Assert(!mpListener); mpListener = pNode; }

        String const& getName() const { return mName; }
        String& getName() { return mName; }

    protected:
        void setPath( HyperKmerPath const* pHKP )
        { Assert(!mpHKP); mpHKP = pHKP;
          if ( mpListener ) mpListener->pathAvailable(this); }
        virtual void pathAvailable( Node* pNode ) {}

    private:
        Node( Node const& ); // unimplemented:  no copying
        Node& operator=( Node const& ); // unimplemented:  no copying

        HyperKmerPath const* mpHKP;
        std::set<Node*> mAdjacencies;
        Node* mpListener;
        String mName;
    };

    struct Comparator
    {
      bool operator()( Node* p1, Node* p2 )
      { return p1->getAdjacencyCount() > p2->getAdjacencyCount() ||
           p1->getAdjacencyCount() == p2->getAdjacencyCount() && p1 > p2; }
    };

    class Checkpointer : public LockedData
    {
    public:
        void addHKP( HyperKmerPath const& hkp ) { mHKPs.insert(&hkp); }
        void removeHKP( HyperKmerPath const& hkp ) { mHKPs.erase(&hkp); }
        HyperKmerPath checkpoint();

    private:
        std::set<HyperKmerPath const*> mHKPs;
    };

    class OutputAccumulator : public Node
    {
    public:
        OutputAccumulator( Checkpointer& chkPtr, size_t maxMerges )
        : mChkPtr(chkPtr) { mToMerge.reserve(maxMerges); }
        virtual ~OutputAccumulator();
        void addNode( Node* pNode );
        HyperKmerPath checkpoint()
        { Locker locker(mChkPtr); return mChkPtr.checkpoint(); }

    protected:
        void internalPathAvailable( Node* pNode );
        size_t getNPending() { return mPending.size(); }
        vec<HyperKmerPath> const& getMerges() { return mToMerge; }
        LockedData& getLock() { return mLock; }

    private:
        void addMerge( Node* pNode );

        Checkpointer& mChkPtr;
        vec<HyperKmerPath> mToMerge;
        std::set<Node*> mPending;
        LockedData mLock;
    };

    class LocalMerger : public OutputAccumulator
    {
    public:
        LocalMerger( Checkpointer& chkPtr, size_t maxMerges, HKPMerger& hkpm,
                     int minOverlap, int minProperOverlap )
        : OutputAccumulator(chkPtr,maxMerges), mHKPM(hkpm),
          mMinOverlap(minOverlap), mMinProperOverlap(minProperOverlap)
        {}

        void addNode( Node* pNode );
        void schedule();
        void doMerge();

    protected:
        virtual void pathAvailable( Node* pNode );

    private:
        HKPMerger& mHKPM;
        int mMinOverlap;
        int mMinProperOverlap;
    };

    class FinalMerger : public OutputAccumulator
    {
    public:
        FinalMerger( Checkpointer& chkPtr, size_t maxMerges );
        ~FinalMerger();
        bool wait( long checkpointIntervalNSecs );

    protected:
        virtual void pathAvailable( Node* pNode );

    private:
        Condition mCondVar; // signals that everyone has reported back
    };

    struct Processor
    {
        void operator()( LocalMerger* pLocalMerger )
        { pLocalMerger->doMerge(); }
    };

    int mMaxQSize;
    NegativeGapValidator const* mNGV;
    vec<tagged_rpint> const& mUniqDB;
    std::vector<Node*> mNodes;
    long mCheckpointInterval;
    String mCheckpointFile;
    Processor mProcessor;
    Worklist<LocalMerger*,Processor> mWorklist;
    friend class LocalMerger;
};

#endif /* HKPMERGER_H_ */
