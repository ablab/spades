///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ReadPather.h
 * \author tsharpe
 * \date May 3, 2010
 *
 * \brief Kmerization and Unipathing of a mess of sequences.
 *
 */
#ifndef KMERS_READ_PATHER_H_
#define KMERS_READ_PATHER_H_

#include "Basevector.h"
#include "String.h"
#include "feudal/BinaryStream.h"
#include "feudal/HashSet.h"
#include "feudal/HugeBVec.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include "kmers/KMer.h"
#include "system/Assert.h"
#include "system/ID.h"
#include "system/System.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <ostream>
#include <unistd.h>
#include <vector>

/// An ID for a KMer.
/// Kmers are assigned IDs in such a way that each edge of a unipath graph can
/// be described by a consecutive sequence of KmerIDs.
class KmerID : public ID<5>
{
public:
    KmerID() {}
    explicit KmerID( size_t id ) : ID<5>(id) {}

    // compiler-supplied copying and destructor are ok
};

TRIVIALLY_SERIALIZABLE(KmerID);

/// An ID for an edge of a unipath graph.
class EdgeID
{
public:
    EdgeID() : mID(NULLVAL) {}
    explicit EdgeID( size_t id ) { setVal(id); }

    // compiler-supplied copying and destructor are ok

    size_t val() const { return mID; }
    void setVal( size_t id ) { ForceAssertLe(id,NULLVAL); mID = id; }
    bool isNull() const { return mID == NULLVAL; }

    friend bool operator==( EdgeID const& eid1, EdgeID const& eid2 )
    { return eid1.mID == eid2.mID; }

    friend bool operator!=( EdgeID const& eid1, EdgeID const& eid2 )
    { return eid1.mID != eid2.mID; }

    friend bool operator<( EdgeID const& eid1, EdgeID const& eid2 )
    { return eid1.mID < eid2.mID; }

    friend int compare( EdgeID const& eid1, EdgeID const& eid2 )
    { return compare(eid1.mID,eid2.mID); }

    friend ostream& operator<<( ostream& os, EdgeID const& eid )
    { return os << eid.mID; }

private:
    static unsigned const NULLVAL = ~0u;
    unsigned mID;
};

TRIVIALLY_SERIALIZABLE(EdgeID);

/// An ID for a connected set of edges in a unipath graph.
class ComponentID : public ID<2>
{
public:
    ComponentID() {}
    explicit ComponentID( size_t id ) : ID<2>(id) {}

    // compiler-supplied copying and destructor are ok
};

TRIVIALLY_SERIALIZABLE(ComponentID);

/// A length along an edge, given as a number of kmers.
class EdgeOffset
{
public:
    EdgeOffset() : mOff(0) {}
    explicit EdgeOffset( size_t off ) { setVal(off); }

    // compiler-supplied copying and destructor are ok

    size_t val() const { return mOff; }
    void setVal( size_t off ) { ForceAssertLe(off,NULLVAL); mOff = off; }
    void threadsafeIncrement()
    { unsigned oldVal;
      do oldVal = mOff;
      while ( oldVal != NULLVAL &&
              !__sync_bool_compare_and_swap(&mOff,oldVal,oldVal+1) ); }

    bool isNull() const { return mOff == NULLVAL; }

    friend bool operator==( EdgeOffset const& eo1, EdgeOffset const& eo2 )
    { return eo1.mOff == eo2.mOff; }

    friend bool operator!=( EdgeOffset const& eo1, EdgeOffset const& eo2 )
    { return eo1.mOff != eo2.mOff; }

    friend bool operator<( EdgeOffset const& eo1, EdgeOffset const& eo2 )
    { return eo1.mOff < eo2.mOff; }

    friend int compare( EdgeOffset const& eo1, EdgeOffset const& eo2 )
    { return compare(eo1.mOff,eo2.mOff); }

    friend ostream& operator<<( ostream& os, EdgeOffset const& eo )
    { return os << eo.mOff; }

private:
    static unsigned const NULLVAL = ~0u;
    unsigned mOff;
};

TRIVIALLY_SERIALIZABLE(EdgeOffset);

/// What you get if you look up a KMer that exists in the dictionary.
/// Namely, you get an EdgeID in which it occurs, and the offset of the kmer
/// in question.
class KDef
{
public:
    KDef() {}

    // compiler-supplied copying and destructor are ok

    EdgeID const& getEdgeID() const { return mEdgeID; }
    EdgeOffset const& getEdgeOffset() const { return mEdgeOffset; }

    bool isNull() const { return mEdgeID.isNull(); }
    void set( EdgeID const& edgeID, EdgeOffset const& edgeOffset )
    { Assert(isNull()); mEdgeID = edgeID; mEdgeOffset = edgeOffset; }

    // the UnipathGraph stores a kmer count in mEdgeOffset during construction
    size_t getCount() const { return mEdgeOffset.val(); }

    void incrCount() { mEdgeOffset.threadsafeIncrement(); }

    friend int compare( KDef const& rl1, KDef const& rl2 )
    { int result = compare(rl1.mEdgeID,rl2.mEdgeID);
      if ( !result ) result = compare(rl1.mEdgeOffset,rl2.mEdgeOffset);
      return result; }

    friend bool operator<( KDef const& rl1, KDef const& rl2 )
    { return compare(rl1,rl2) < 0; }

    friend bool operator<=( KDef const& rl1, KDef const& rl2 )
    { return compare(rl1,rl2) <= 0; }

    friend bool operator>( KDef const& rl1, KDef const& rl2 )
    { return compare(rl1,rl2) > 0; }

    friend bool operator>=( KDef const& rl1, KDef const& rl2 )
    { return compare(rl1,rl2) >= 0; }

    friend bool operator==( KDef const& kd1, KDef const& kd2 )
    { return !memcmp(&kd1,&kd2,sizeof(KDef)); }

    friend bool operator!=( KDef const& kd1, KDef const& kd2 )
    { return memcmp(&kd1,&kd2,sizeof(KDef)); }

    friend ostream& operator<<( ostream& os, KDef const& kd )
    { os << '[' << kd.mEdgeID << '+' << kd.mEdgeOffset << ']'; return os; }

private:
    EdgeID mEdgeID;
    EdgeOffset mEdgeOffset;
    static unsigned short const MAX_COUNT = ~0;
};

TRIVIALLY_SERIALIZABLE(KDef);

template <unsigned K>
class KmerDictEntry : public KMer<K>
{
public:
    KmerDictEntry() {}
    KmerDictEntry( KMer<K> const& kmer ) : KMer<K>(kmer) {}

    // compiler-supplied copying and destructor are OK

    KDef& getKDef() { return mKDef; }
    KDef const& getKDef() const { return mKDef; }

private:
    KDef mKDef;
};

template <unsigned K>
struct Serializability< KmerDictEntry<K> > : public TriviallySerializable {};

/// A set of KmerDictEntry's, used as a map from KMer onto KDef.
template <unsigned K>
class KmerDict
{
public:
    typedef KmerDictEntry<K> Entry;
    typedef typename KMer<K>::Hasher Hasher;
    typedef std::equal_to<KMer<K> > Comparator;
    typedef HashSet<Entry,Hasher,Comparator> Set;
    typedef typename Set::const_iterator OCItr;
    typedef typename Set::iterator OItr;
    typedef typename Set::ICItr ICItr;

    KmerDict( size_t dictSize )
    : mKSet(dictSize) {}

    // compiler-supplied destructor is OK

    Entry const* findEntryCanonical( KMer<K> const& kmer ) const
    { return mKSet.lookup(kmer); }

    /// Returns null pointer if kmer isn't in dictionary.
    Entry const* findEntry( KMer<K> const& kmer ) const
    { return mKSet.lookup( kmer.getCanonicalForm() == REV ?
                            KMer<K>(kmer).rc() :
                            kmer ); }

    /// Returns null pointer if kmer isn't in dictionary.
    KDef* lookup( KMer<K> const& kmer )
    { Entry const* pEnt = findEntry(kmer);
      return pEnt ? const_cast<KDef*>(&pEnt->getKDef()) : 0; }

    KDef const* lookup( KMer<K> const& kmer ) const
    { Entry const* pEnt = findEntry(kmer);
      return pEnt ? &pEnt->getKDef() : 0; }

    /// Canonicalizes and inserts kmer, if necessary.
    KDef& operator[]( KMer<K> const& kmer )
    { Entry const* pEnt;
      if ( kmer.getCanonicalForm() != REV ) pEnt = &mKSet[kmer];
      else pEnt = &mKSet[KMer<K>(kmer).rc()];
      return const_cast<KDef&>(pEnt->getKDef()); }

    /// Inserts kmer, if necessary.  Kmer is assumed to be canonical.
    KDef& refCanonical( KMer<K> const& kmer )
    { return const_cast<KDef&>(mKSet[kmer].getKDef()); }

    /// Inserts a kmer known to be in canonical form.
    void insertCanonical( KMer<K> const& kmer )
    { mKSet.add(kmer); }

    OCItr begin() const { return mKSet.begin(); }
    OCItr end() const { return mKSet.end(); }
    OCItr cbegin() { return mKSet.cbegin(); }
    OCItr cend() { return mKSet.cend(); }
    OItr begin() { return mKSet.begin(); }
    OItr end() { return mKSet.end(); }

    size_t size() const { return mKSet.size(); }

    void process( String const& fastb, unsigned nThreads );
    void clear() { mKSet.clear(); }

    friend void swap( KmerDict& kd1, KmerDict& kd2 )
    { swap(kd1.mKSet,kd2.mKSet); }

    friend bool operator==( KmerDict const& kd1, KmerDict const& kd2 )
    { return kd1.mKSet == kd2.mKSet; }

    friend bool operator!=( KmerDict const& kd1, KmerDict const& kd2 )
    { return !(kd1==kd2); }

    size_t writeBinary( BinaryWriter& writer ) const
    { return writer.write(mKSet); }

    void readBinary( BinaryReader& reader )
    { reader.read(&mKSet); }

    static size_t externalSizeof() { return 0; }

private:
    KmerDict( KmerDict const& ); // unimplemented -- no copying
    KmerDict& operator=( KmerDict const& ); // unimplemented -- no copying

    Set mKSet;
};

template <unsigned K>
struct Serializability< KmerDict<K> > : public SelfSerializable {};

/// An unbranched sequence of kmers.
/// Described by a close-ended range of kmer IDs.
/// Also includes info on predecessor and successor edges.
class UnipathEdge
{
public:
    UnipathEdge() : mRCFlags(0), mIsPalindrome(false) {}

    UnipathEdge( KmerID const& kmerID, ComponentID componentID,
                 bool isPalindrome )
    : mRCFlags(0),
      mFirstKmerID(kmerID), mLastKmerID(kmerID),
      mComponentID(componentID), mIsPalindrome(isPalindrome)
    {}

    // compiler-supplied copying and destructor are ok

    KmerID const& getInitialKmerID() const { return mFirstKmerID; }
    KmerID const& getFinalKmerID() const { return mLastKmerID; }
    KmerID getKmerID( EdgeOffset const& offset ) const
    { return KmerID(mFirstKmerID.val()+offset.val()); }

    size_t getLength() const { return mLastKmerID.val()-mFirstKmerID.val()+1; }
    ComponentID const& getComponentID() const { return mComponentID; }

    size_t extend()
    { size_t newVal = mLastKmerID.val()+1;
      mLastKmerID.setVal(newVal);
      return newVal - mFirstKmerID.val(); }

    EdgeID const& getPredecessor( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return mConnections[baseCode]; }

    bool isPredecessorRC( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return (mRCFlags >> baseCode) & 1; }

    void setPredecessor( unsigned baseCode, EdgeID id, bool rc )
    { AssertLt(baseCode,4u);
      if ( !mConnections[baseCode].isNull() )
      { AssertEq(mConnections[baseCode],id); }
      else
      { mConnections[baseCode] = id;
        if ( rc ) mRCFlags |= 1 << baseCode; } }

    EdgeID const& getSuccessor( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return mConnections[baseCode+4]; }

    bool isSuccessorRC( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return (mRCFlags >> (baseCode+4)) & 1; }

    void setSuccessor( unsigned baseCode, EdgeID id, bool rc )
    { AssertLt(baseCode,4u);
      baseCode += 4u;
      if ( !mConnections[baseCode].isNull() )
      { AssertEq(mConnections[baseCode],id); }
      else
      { mConnections[baseCode] = id;
        if ( rc ) mRCFlags |= 1 << baseCode; } }

    bool isPalindrome() const
    { return mIsPalindrome; }

    void setPalindrome( bool isPalindrome ) { mIsPalindrome = isPalindrome; }
    friend bool operator==( UnipathEdge const& ue1, UnipathEdge const& ue2 )
    { return !memcmp(&ue1,&ue2,sizeof(UnipathEdge)); }

    friend bool operator!=( UnipathEdge const& ue1, UnipathEdge const& ue2 )
    { return !(ue1 == ue2); }

    friend ostream& operator<<( ostream& os, UnipathEdge const& edge )
    { os << '[' << edge.getInitialKmerID() << '-' << edge.getFinalKmerID() <<
        '@' << edge.getComponentID();
      char sep = '>';
      for ( size_t idx = 0; idx < 8; ++idx )
      { os << sep; sep = ',';
        if ( edge.mRCFlags & (1 << idx) ) os << '~';
        os << edge.mConnections[idx]; }
      return os << ']'; }

private:
    unsigned char mRCFlags;
    KmerID mFirstKmerID;
    KmerID mLastKmerID;
    ComponentID mComponentID;
    EdgeID mConnections[8];
    bool mIsPalindrome;
};

TRIVIALLY_SERIALIZABLE(UnipathEdge);
typedef std::vector<UnipathEdge> UnipathEdgeVec;

/// A bi-directional graph of UnipathEdges.
template <unsigned K>
class UnipathGraph
{
public:
    typedef KmerDict<K> KDict;

    /// Create a new graph by kmerizing a mess of reads.
    static void create( String const& graphInfoFilename,
                        String const& fastb,
                        unsigned nThreads,
                        size_t nKmersEstimate );

    /// Create a new graph from a kmer dictionary.
    static void create( String const& graphInfoFilename, KDict& dict );

    /// Instantiate a graph from its binary files.
    /// If you pass a KDict, you're relinquishing control over it to the
    /// graph -- don't delete it yourself!
    UnipathGraph( String const& graphInfoFile, KDict* pDict = 0 );

    ~UnipathGraph() { ungetDict(); }

    size_t getNKmers() const { return mNKmers; }
    size_t getNEdges() const { return mEdges.size(); }

    UnipathEdge const& getEdge( EdgeID const& edgeID ) const
    { Assert(edgeID.val()<mEdges.size()); return mEdges[edgeID.val()]; }

    UnipathEdgeVec const& getAllEdges() const { return mEdges; }

    HugeBVec::const_iterator getBases( KmerID const& kmerID ) const
    { return mSeq.begin(kmerID.val()); }
    HugeBVec::const_iterator getBases( EdgeID const& edgeID ) const
    { return mSeq.begin(getEdge(edgeID).getInitialKmerID().val()); }
    HugeBVec::const_iterator getBases( KDef const& kDef ) const
    { return mSeq.begin(getEdge(kDef.getEdgeID()).getInitialKmerID().val()
                        +kDef.getEdgeOffset().val()); }
    HugeBVec::const_iterator getBasesEnd( KmerID const& kmerID ) const
    { return mSeq.begin(kmerID.val()+mK); }
    HugeBVec::const_iterator getBasesEnd( EdgeID const& edgeID ) const
    { return mSeq.begin(getEdge(edgeID).getFinalKmerID().val()+mK); }
    HugeBVec::const_iterator getBasesEnd( KDef const& kDef ) const
    { return mSeq.begin(getEdge(kDef.getEdgeID()).getInitialKmerID().val()
                        +kDef.getEdgeOffset().val()+mK); }

    HugeBVec::const_rc_iterator getRCBases( KmerID const& kmerID ) const
    { return mSeq.rcbegin(kmerID.val()+mK); }
    HugeBVec::const_rc_iterator getRCBases( EdgeID const& edgeID ) const
    { return mSeq.rcbegin(getEdge(edgeID).getFinalKmerID().val()+mK); }
    HugeBVec::const_rc_iterator getRCBases( KDef const& kDef ) const
    { return mSeq.rcbegin(getEdge(kDef.getEdgeID()).getInitialKmerID().val()
                        +kDef.getEdgeOffset().val()+mK); }
    HugeBVec::const_rc_iterator getRCBasesEnd( KmerID const& kmerID ) const
    { return mSeq.rcbegin(kmerID.val()+mK,mK); }
    HugeBVec::const_rc_iterator getRCBasesEnd( EdgeID const& edgeID ) const
    { UnipathEdge const& edge = getEdge(edgeID);
      size_t len = edge.getLength()+mK-1;
      return mSeq.rcbegin(edge.getFinalKmerID().val()+mK,len); }
    HugeBVec::const_rc_iterator getRCBasesEnd( KDef const& kDef ) const
    { return mSeq.rcbegin(getEdge(kDef.getEdgeID()).getInitialKmerID().val()
                        +kDef.getEdgeOffset().val()+mK,mK); }

    HugeBVec const& getAllBases() const { return mSeq; }

    bool nextEdge( unsigned baseCode, EdgeID* pEdgeID, bool isRC ) const
    { UnipathEdge const& edge = getEdge(*pEdgeID);
      if ( isRC )
      { *pEdgeID = edge.getPredecessor(baseCode);
        if ( edge.isPredecessorRC(baseCode) ) isRC = false; }
      else
      { *pEdgeID = edge.getSuccessor(baseCode);
        if ( edge.isSuccessorRC(baseCode) ) isRC = true; }
      AssertNot(pEdgeID->isNull());
      return isRC; }

    // loads dictionary (which is huge) lazily
    KDict const& getDict() const
    { if ( !mpDict ) loadDict(); return *mpDict; }

    // unloads dictionary to recover lots of memory
    void ungetDict() const
    { delete mpDict; mpDict = 0; }

    void validate() { validateUnipaths(mSeq,mEdges,getDict()); }

    friend bool operator==( UnipathGraph const& ug1, UnipathGraph const& ug2 )
    { return ug1.mSeq == ug2.mSeq &&
             ug1.mEdges == ug2.mEdges; }

    friend bool operator!=( UnipathGraph ug1, UnipathGraph ug2 )
    { return !(ug1 == ug2); }

    static String getInfoFilename( String const& fastb )
    { return fastb.ReplaceExtension(".fastb",".k"+ToString(K)+".ug.info"); }

    static String getSeqFilename( String const& graphInfoFilename )
    { return graphInfoFilename.ReplaceExtension(".info",".seq"); }
    static String getEdgesFilename( String const& graphInfoFilename )
    { return graphInfoFilename.ReplaceExtension(".info",".edges"); }
    static String getDictFilename( String const& graphInfoFilename )
    { return graphInfoFilename.ReplaceExtension(".info",".dict"); }
    static String getCountsFilename( String const& graphInfoFilename )
    { return graphInfoFilename.ReplaceExtension(".info",".counts"); }
    static void removeFiles( String const& graphInfoFilename )
    { unlink(graphInfoFilename.c_str());
      unlink(getSeqFilename(graphInfoFilename).c_str());
      unlink(getEdgesFilename(graphInfoFilename).c_str());
      unlink(getDictFilename(graphInfoFilename).c_str());
      unlink(getCountsFilename(graphInfoFilename).c_str()); }

    static void validateUnipaths( HugeBVec const& seq,
                                  UnipathEdgeVec const& edges,
                                  KDict const& dict );

private:
    UnipathGraph( UnipathGraph const& ); // unimplemented -- no copying
    UnipathGraph& operator=( UnipathGraph const& ); // unimplemented -- no copying

    void loadDict() const;

    String mGraphInfoFile;
    HugeBVec mSeq;
    UnipathEdgeVec mEdges;
    mutable KDict* mpDict;
    size_t mNKmers;
    unsigned mK;
};

// Data about the graph files stashed in the .info file.
struct GraphInfo
{
    size_t nComponents;
    size_t nKmers;
    size_t nEdges;
    size_t kSeqSize;
    unsigned K;
};

TRIVIALLY_SERIALIZABLE(GraphInfo);

template <unsigned K>
UnipathGraph<K>::UnipathGraph( String const& graphInfoFile, KmerDict<K>* pDict )
: mGraphInfoFile(graphInfoFile), mSeq(getSeqFilename(graphInfoFile).c_str()),
  mpDict(pDict)
{
    GraphInfo gi;
    BinaryReader::readFile(graphInfoFile.c_str(),&gi);

    mNKmers = gi.nKmers;
    mK = gi.K;

    if ( mSeq.size() != gi.kSeqSize )
        FatalErr("Length of kmer bases sequence (" << mSeq.size() <<
                 ") doesn't jibe with ug.info's value (" << gi.kSeqSize <<
                 ").");

    String edgesFile = getEdgesFilename(graphInfoFile);
    BinaryReader::readFile(edgesFile.c_str(),&mEdges);
    if ( mEdges.size() != gi.nEdges )
        FatalErr("Number of graph edges (" << mEdges.size() <<
                 ") doesn't jibe with ug.info's value (" << gi.nEdges << ").");
}

/// Figures out successors or predecessors to a given kmer.
template <unsigned K>
class KmerStepper
{
public:
    typedef typename KmerDict<K>::Entry DictEntry;

    KmerStepper( KmerDict<K> const& dict )
    : mDict(dict) {}

    // compiler-supplied copy ctor and dtor are OK

    KMer<K>& getSteppedKmer() { return mKMer; }

    unsigned getSuccessors( KMer<K> const& kmer, DictEntry const** pEntries )
    { unsigned result = 0;
      mKMer = kmer;
      mKMer.toSuccessor(0);
      for ( unsigned succCode = 0; succCode < 4u; ++succCode )
      { mKMer.setBack(succCode);
        KMer<K> const& kkk = mKMer.getCanonicalForm() == REV ?
                KMer<K>(mKMer).rc() : mKMer;
        if ( (*pEntries++ = mDict.findEntryCanonical(kkk)) ) result += 1; }
      return result; }

    unsigned getPredecessors( KMer<K> const& kmer, DictEntry const** pEntries )
    { unsigned result = 0;
      mKMer = kmer;
      mKMer.toPredecessor(0);
      for ( unsigned predCode = 0; predCode < 4u; ++predCode )
      { mKMer.setFront(predCode);
        KMer<K> const& kkk = mKMer.getCanonicalForm() == REV ?
              KMer<K>(mKMer).rc() : mKMer;
        if ( (*pEntries++ = mDict.findEntryCanonical(kkk)) ) result += 1; }
      return result; }

private:
    KmerDict<K> const& mDict;
    KMer<K> mKMer;
};

template <unsigned K>
void UnipathGraph<K>::validateUnipaths( HugeBVec const& bv,
                                     UnipathEdgeVec const& edges,
                                     KmerDict<K> const& dict )
{
    std::cout << Date() << " Validating unipaths." << std::endl;
    size_t nErrors = 0;
    typedef UnipathEdgeVec::const_iterator Itr;
    typename KmerDict<K>::Entry const* entries[4];
    KmerStepper<K> stepper(dict);
    KMer<K> kmer;
    using std::equal;
    for ( Itr itr(edges.begin()), end(edges.end()); itr != end; ++itr )
    {
        UnipathEdge const& edge = *itr;
        kmer.assign(bv.begin(edge.getFinalKmerID().val()));
        stepper.getSuccessors(kmer,entries);
        for ( unsigned succCode = 0; succCode < 4u; ++succCode )
        {
            typename KmerDict<K>::Entry const* pEntry = entries[succCode];
            EdgeID edgeID = edge.getSuccessor(succCode);
            if ( !pEntry )
            {
                if ( !edgeID.isNull() )
                {
                    std::cout << "Edge " << (itr-edges.begin()) << " has a " <<
                        Base::val2Char(succCode) << " successor, edge " <<
                        edgeID << ", that wasn't found by getSuccessors." <<
                        std::endl;
                    nErrors += 1;
                }
            }
            else
            {
                KDef const& def = pEntry->getKDef();
                if ( edgeID.isNull() )
                {
                    std::cout << "Edge " << (itr-edges.begin()) <<
                            " should have a " << Base::val2Char(succCode) <<
                            " successor, edge " << def.getEdgeID() <<
                            ", but it's not marked." << std::endl;
                    nErrors += 1;
                }
                else
                {
                    if ( edgeID != def.getEdgeID() )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) <<
                                " should have a " << Base::val2Char(succCode) <<
                                " successor, edge " << def.getEdgeID() <<
                                ", but it's marked as edge " << edgeID << '.'
                                << std::endl;
                        nErrors += 1;
                        edgeID = def.getEdgeID();
                    }
                    UnipathEdge const& succ = edges[edgeID.val()];
                    if ( edge.getComponentID() != succ.getComponentID() )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) << " has a "
                                << Base::val2Char(succCode) << " successor, " <<
                                edgeID << ", in a different component: " <<
                                edge.getComponentID() << " vs. " <<
                                succ.getComponentID() << std::endl;
                        nErrors += 1;
                    }
                    bool isRC = edge.isSuccessorRC(succCode);
                    KmerID kid = isRC ? succ.getFinalKmerID() :
                                        succ.getInitialKmerID();
                    KmerID kid2 = succ.getKmerID(def.getEdgeOffset());
                    if ( kid != kid2 )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) << " has a "
                                << Base::val2Char(succCode) << " successor, " <<
                                edgeID << ", but its adjacent kmerID " << kid <<
                                " looks up as " << kid2 << '.' << std::endl;
                        nErrors += 1;
                    }
                    KMer<K>& kmerSucc = stepper.getSteppedKmer();
                    kmerSucc.setBack(succCode);
                    if ( isRC )
                    {
                        if ( !equal(kmerSucc.rcbegin(),kmerSucc.rcend(),
                                    bv.begin(kid.val())) )
                        {
                            std::cout << "Edge " << (itr-edges.begin()) <<
                                    " has a " << Base::val2Char(succCode) <<
                                    " successor, " << edgeID <<
                           ", but the adjacent kmer has an incorrect sequence."
                                    << std::endl;
                            nErrors += 1;
                        }
                    }
                    else
                    {
                        if ( !equal(kmerSucc.begin(),kmerSucc.end(),
                                    bv.begin(kid.val())) )
                        {
                            std::cout << "Edge " << (itr-edges.begin()) <<
                                    " has a " << Base::val2Char(succCode) <<
                                    " successor, " << edgeID <<
                           ", but the adjacent kmer has an incorrect sequence."
                                    << std::endl;
                            nErrors += 1;
                        }
                    }
                }
            }
        }

        kmer.assign(bv.begin(edge.getInitialKmerID().val()));
        stepper.getPredecessors(kmer,entries);
        for ( unsigned predCode = 0; predCode < 4u; ++predCode )
        {
            typename KmerDict<K>::Entry const* pEntry = entries[predCode];
            EdgeID edgeID = edge.getPredecessor(predCode);
            if ( !pEntry )
            {
                if ( !edgeID.isNull() )
                {
                    std::cout << "Edge " << (itr-edges.begin()) << " has a " <<
                        Base::val2Char(predCode) << " predecessor, edge " <<
                        edgeID << ", that wasn't found by getPredecessors." <<
                        std::endl;
                    nErrors += 1;
                }
            }
            else
            {
                KDef const& def = pEntry->getKDef();
                if ( edgeID.isNull() )
                {
                    std::cout << "Edge " << (itr-edges.begin()) <<
                            " should have a " << Base::val2Char(predCode) <<
                            " predecessor, edge " << def.getEdgeID() <<
                            ", but it's not marked." << std::endl;
                    nErrors += 1;
                }
                else
                {
                    if ( edgeID != def.getEdgeID() )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) <<
                                " should have a " << Base::val2Char(predCode) <<
                                " predecessor, edge " << def.getEdgeID() <<
                                ", but it's marked as edge " << edgeID << '.'
                            << std::endl;
                        nErrors += 1;
                        edgeID = def.getEdgeID();
                    }
                    UnipathEdge const& pred = edges[edgeID.val()];
                    if ( edge.getComponentID() != pred.getComponentID() )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) << " has a "
                                << Base::val2Char(predCode) << " predecessor, "
                                << edgeID << ", in a different component: " <<
                                edge.getComponentID() << " vs. " <<
                                pred.getComponentID() << std::endl;
                        nErrors += 1;
                    }
                    bool isRC = edge.isPredecessorRC(predCode);
                    KmerID kid = isRC ? pred.getInitialKmerID():
                                        pred.getFinalKmerID();
                    KmerID kid2 = pred.getKmerID(def.getEdgeOffset());
                    if ( kid != kid2 )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) << " has a "
                                << Base::val2Char(predCode) << " predecessor, "
                                << edgeID << ", but its adjacent kmerID " <<
                                kid << " looks up as " << kid2 << '.' <<
                                std::endl;
                        nErrors += 1;
                    }
                    KMer<K>& kmerPred = stepper.getSteppedKmer();
                    kmerPred.setFront(predCode);
                    if ( isRC )
                    {
                        if ( !equal(kmerPred.rcbegin(),kmerPred.rcend(),
                                    bv.begin(kid.val())) )
                        {
                            std::cout << "Edge " << (itr-edges.begin()) <<
                                    " has a " << Base::val2Char(predCode) <<
                                    " predecessor, " << edgeID <<
                           ", but the adjacent kmer has an incorrect sequence."
                                    << std::endl;
                            nErrors += 1;
                        }
                    }
                    else
                    {
                        if ( !equal(kmerPred.begin(),kmerPred.end(),
                                    bv.begin(kid.val())) )
                        {
                            std::cout << "Edge " << (itr-edges.begin()) <<
                                    " has a " << Base::val2Char(predCode) <<
                                    " predecessor, " << edgeID <<
                           ", but the adjacent kmer has an incorrect sequence."
                                    << std::endl;
                            nErrors += 1;
                        }
                    }
                }
            }
        }
    }
    ForceAssertEq(nErrors,0ul);
}

template <unsigned K>
void UnipathGraph<K>::loadDict() const
{
    std::cout << Date() << " Loading kmer dictionary." << std::endl;
    GraphInfo gi;
    BinaryReader::readFile(mGraphInfoFile.c_str(),&gi);
    mpDict = new KmerDict<K>(0);
    BinaryReader::readFile(getDictFilename(mGraphInfoFile).c_str(),mpDict);
    ForceAssertEq(gi.nKmers,mpDict->size());
}

/// An ID for a Read.
class ReadID : public ID<5>
{
public:
    ReadID() {}
    explicit ReadID( size_t id ) : ID<5>(id) {}

    // compiler-supplied copying and destructor are ok
};

TRIVIALLY_SERIALIZABLE(ReadID);

/// Describes a read that traverses some unipath graph edge.
/// (Which edge is not specified by this class: that's an external association.)
/// Also tells you which segment of the read's pathing has the edge, and whether
/// it's reverse-complement or not.
class UnipathEvidence
{
public:
    UnipathEvidence() : mSegID(0), mIsRC(false) {}
    UnipathEvidence( ReadID const& readID, unsigned segID, bool isRC )
    : mReadID(readID), mSegID(segID), mIsRC(isRC)
    { ForceAssertLt(segID,1u<<SEGBITS); }

    // compiler-supplied copying and destructor are ok

    ReadID const& getReadID() const { return mReadID; }
    unsigned getSegmentID() const { return mSegID; }
    bool isRC() const { return mIsRC; }

    friend int compare( UnipathEvidence const& ev1, UnipathEvidence const& ev2 )
    { int result = compare(ev1.mReadID,ev2.mReadID);
      if ( !result ) result = compare(ev1.mSegID,ev2.mSegID);
      return result; }

    static unsigned const SEGBITS = 15;
    static unsigned const MAX_NSEGMENTS = 1u<<SEGBITS;

private:
    ReadID mReadID;
    unsigned short mSegID : SEGBITS;
    unsigned short mIsRC : 1;
};

typedef SerfVec<UnipathEvidence> UnipathEvidenceVec;
typedef MasterVec<UnipathEvidenceVec> VecUnipathEvidenceVec;

TRIVIALLY_SERIALIZABLE(UnipathEvidence);

/// An index into a sequence of bases.  Each base shows which turn to take at
/// each node of a UnipathGraph traversed by some sequence of bases (usually a
/// read) placed onto the graph.  This is a very terse way of representing a
/// path.
class PathID : public ID<5>
{
public:
    PathID() {}
    explicit PathID( size_t id ) : ID<5>(id) {}

    // compiler-supplied copying and destructor are ok
};

TRIVIALLY_SERIALIZABLE(PathID);

/// Describes a traversal of a UnipathGraph.
/// Used, for example, to indicate where a read is placed on the graph.
/// This is not convenient to compute with, but is a good storage format.
class Unipathing
{
public:
    Unipathing() : mIsRC(false), mNSegments(0) {}
    Unipathing( EdgeID const& edgeID, bool isRC, size_t skip )
    : mEdgeID(edgeID), mInitialSkip(skip), mIsRC(isRC), mNSegments(1)
    {}

    // compiler-supplied copying and destructor are ok

    bool isNull() const { return !mNSegments; }
    unsigned getNSegments() const { return mNSegments; }
    EdgeID const& getInitialEdgeID() const { return mEdgeID; }
    bool isInitialEdgeRC() const { return mIsRC; }
    PathID const& getPathID() const { return mPathID; }
    size_t getInitialSkip() const { return mInitialSkip.val(); }
    size_t getFinalSkip() const { return mFinalSkip.val(); }

    void setSegments( unsigned nSegments, PathID const& pathID )
    { ForceAssertLt(nSegments,32768u); mNSegments = nSegments;
      mPathID = pathID; }

    void setFinalSkip( size_t skip )
    { mFinalSkip = EdgeOffset(skip); }

    friend std::ostream& operator<<( std::ostream& os,
                                     Unipathing const& pathing )
    { os << "Edge=" << pathing.mEdgeID <<
            " IsRC=" << pathing.mIsRC <<
            " IniSkip=" << pathing.mInitialSkip <<
            " FinSkip=" << pathing.mFinalSkip <<
            " nSegs=" << static_cast<unsigned>(pathing.mNSegments) <<
            " PathID=" << pathing.mPathID;
      return os; }

private:
    EdgeID mEdgeID; // initial edge ID
    EdgeOffset mInitialSkip;
    EdgeOffset mFinalSkip;
    unsigned short mIsRC : 1; // initial edge is RC
    unsigned short mNSegments : 15;
    PathID mPathID;
};

TRIVIALLY_SERIALIZABLE(Unipathing);
typedef std::vector<Unipathing> UnipathingVec;

class EdgeDesc
{
public:
    EdgeDesc( EdgeID const& edgeID, CanonicalForm status )
    : mEdgeID(edgeID), mStatus(status) {}

    // compiler-supplied copying and destructor are OK

    EdgeID const& getEdgeID() const { return mEdgeID; }
    CanonicalForm getStatus() const { return mStatus; }

    EdgeDesc operator~() const
    { return EdgeDesc(mEdgeID,gCompStatus[mStatus]); }

    friend bool operator<( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) < 0; }
    friend bool operator<=( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) <= 0; }
    friend bool operator>( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) > 0; }
    friend bool operator>=( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) >= 0; }
    friend bool operator==( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) == 0; }
    friend bool operator!=( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) != 0; }
    friend int compare( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { int result = compare(ed1.mEdgeID,ed2.mEdgeID);
      if ( !result ) result = compare(ed1.mStatus,ed2.mStatus);
      return result; }

    friend ostream& operator<<( ostream& os, EdgeDesc const& ed )
    { if ( ed.mStatus == REV ) os << '~';
      return os << ed.mEdgeID; }

private:
    EdgeID mEdgeID;
    CanonicalForm mStatus;

    static CanonicalForm gCompStatus[3];
};

/// Describes a unipath graph traversal as a sequence of EdgeDescs.
/// Also has the number of kmers not covered in the initial and final segments.
class EdgeList : public std::vector<EdgeDesc>
{
    typedef std::vector<EdgeDesc> Base;
public:
    EdgeList()
    : mInitialSkip(0), mFinalSkip(0)
    {}

    EdgeList( size_t initialSkip, size_t finalSkip )
    : mInitialSkip(initialSkip), mFinalSkip(finalSkip)
    {}

    // compiler-supplied copying and destructor are OK

    // Reverse-complement in place.
    EdgeList& rc()
    { using std::swap; swap(mInitialSkip,mFinalSkip);
      iterator head = begin(); iterator tail = end();
      while (head != tail)
      { EdgeDesc tmp = ~*--tail;
        if ( head == tail ) { *head = tmp; break; }
        *tail = ~*head; *head = tmp; ++head; }
      return *this; }

    size_t getInitialSkip() const { return mInitialSkip; }
    size_t getFinalSkip() const { return mFinalSkip; }

    friend bool operator<( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) < 0; }
    friend bool operator<=( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) <= 0; }
    friend bool operator>( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) > 0; }
    friend bool operator>=( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) >= 0; }
    friend bool operator==( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) == 0; }
    friend bool operator!=( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) != 0; }
    friend int compare( EdgeList const& el1, EdgeList const& el2 )
    { int result = compare(static_cast<Base const&>(el1),
                            static_cast<Base const&>(el2));
      if ( !result ) result = compare(el1.mInitialSkip,el2.mInitialSkip);
      if ( !result ) result = compare(el1.mFinalSkip,el2.mFinalSkip);
      return result; }

    friend ostream& operator<<( ostream& os, EdgeList const& el )
    { os << '[' << el.getInitialSkip();
      char sep = ']';
      EdgeList::const_iterator end(el.end());
      for ( EdgeList::const_iterator itr(el.begin()); itr != end; ++itr )
      { os << sep << *itr; sep = ':'; }
      os << '[' << el.getFinalSkip() << ']';
      return os; }

    friend void swap( EdgeList& el1, EdgeList& el2 )
    { el1.swap(el2);
      using std::swap;
      swap(el1.mInitialSkip,el2.mInitialSkip);
      swap(el1.mFinalSkip,el2.mFinalSkip); }

private:
    size_t mInitialSkip;
    size_t mFinalSkip;
};

/// Describes how reads align to a UnipathGraph.
template <unsigned K>
class PathCollection
{
public:
    typedef UnipathGraph<K> UGraph;

    /// Create a new graph, and path onto it.
    static void create( String const& fastbFilename,
                        bool validate,
                        unsigned nThreads,
                        size_t nKmersEst = 0,
                        bool writeKmerCounts = false );

    /// Path onto an existing graph.
    static void create( String const& fastbFilename,
                        String const& graphInfoFilename,
                        bool validate, unsigned nThreads );

    PathCollection( String const& pathInfoFilename,
                    String const& graphInfoFilename );

    // copying prohibited.  compiler-supplied destructor is OK

    UGraph const& getGraph() const { return mGraph; }

    size_t getNReads() const { return mPathings.size(); }

    EdgeList getEdgeList( size_t readId ) const;
    EdgeList getEdgeList( ReadID const& readID ) const
    { return getEdgeList(readID.val()); }

    // number of segments in the path
    size_t getEdgeListSize( size_t readId ) const
    { return mPathings[readId].getNSegments(); }
    size_t getEdgeListSize( ReadID const& readID ) const
    { return mPathings[readID.val()].getNSegments(); }

    // number of kmers in the path
    size_t getEdgeLen( EdgeDesc const& ed ) const
    { return getEdgeLen(ed.getEdgeID()); }
    size_t getEdgeLen( EdgeID const& edgeID ) const
    { return mGraph.getEdge(edgeID).getLength(); }
    size_t getEdgeListLen( EdgeList const& el ) const
    { size_t len = getEdgeListLen(el.begin(),el.end());
      return len - el.getInitialSkip() - el.getFinalSkip(); }
    size_t getEdgeListLen( EdgeList::const_iterator itr,
                           EdgeList::const_iterator end ) const
    { size_t len = 0;
      while ( itr != end )
      { len += getEdgeLen(itr->getEdgeID()); ++itr; }
      return len; }

    EdgeDesc getNextEdgeDesc( EdgeDesc const& ed, unsigned base ) const
    { UnipathEdge const& edge = mGraph.getEdge(ed.getEdgeID());
      EdgeID id;
      CanonicalForm status = FWD;
      if ( ed.getStatus() == REV )
      { id = edge.getPredecessor(base);
        if ( !edge.isPredecessorRC(base) ) status = REV; }
      else
      { id = edge.getSuccessor(base);
        if ( edge.isSuccessorRC(base) ) status = REV; }
      if ( !id.isNull() && mGraph.getEdge(id).isPalindrome() )
        status = PALINDROME;
      return EdgeDesc(id,status); }

    bvec getBases( EdgeList const& el ) const
    { bvec result; getBases(el,&result); return result; }

    bvec& getBases( EdgeList const& el, bvec* pBV ) const;

    void validate( String const& fastbFilename );

    static String getInfoFilename( String const& fastb )
    { return fastb.ReplaceExtension(".fastb",".k"+ToString(K)+".pc.info"); }

    static String getPathsFilename( String const& infoFile )
    { return infoFile.ReplaceExtension(".info",".paths"); }

    static String getPathseqFilename( String const& infoFile )
    { return infoFile.ReplaceExtension(".info",".seq"); }

    static String getEvidenceFilename( String const& infoFile )
    { return infoFile.ReplaceExtension(".info",".ev"); }

private:
    PathCollection( PathCollection const& ); // unimplemented -- no copying
    PathCollection& operator=( PathCollection const& ); // unimplemented -- no copying

    UGraph mGraph;
    HugeBVec mPathSeq; // indexed by PathID
    UnipathingVec mPathings; // indexed by ReadID
};

// Data about the pathing files stuffed into the .info file.
struct PathInfo
{
    size_t nReads;
    size_t pathSeqSize;
};

TRIVIALLY_SERIALIZABLE(PathInfo);

template <unsigned K>
PathCollection<K>::PathCollection( String const& pathInfoFilename,
                                  String const& graphInfoFilename )
: mGraph(graphInfoFilename),
  mPathSeq(getPathseqFilename(pathInfoFilename).c_str())
{
    String pathsFile = getPathsFilename(pathInfoFilename);
    BinaryReader::readFile(pathsFile.c_str(),&mPathings);
    PathInfo pi;
    BinaryReader::readFile(pathInfoFilename.c_str(),&pi);
    if ( mPathings.size() != pi.nReads )
        FatalErr("Number of traversals (" << mPathings.size() <<
                 ") doesn't jibe with pathinfo's value (" << pi.nReads << ").");
    if ( mPathSeq.size() != pi.pathSeqSize )
        FatalErr("Length of traversal bases sequence (" <<
                 mPathSeq.size() <<
                 ") doesn't jibe with pathinfo's value (" << pi.pathSeqSize <<
                 ").");
}

template <unsigned K>
EdgeList PathCollection<K>::getEdgeList( size_t readId ) const
{
    Unipathing const& pathing = mPathings[readId];
    EdgeList result(pathing.getInitialSkip(),pathing.getFinalSkip());
    unsigned nSegs = pathing.getNSegments();
    if ( nSegs )
    {
        result.reserve(nSegs);
        EdgeID edgeID = pathing.getInitialEdgeID();
        bool isRC = pathing.isInitialEdgeRC();
        bool isPalindrome = mGraph.getEdge(edgeID).isPalindrome();
        CanonicalForm status = isPalindrome ? PALINDROME : isRC ? REV : FWD;
        result.push_back(EdgeDesc(edgeID,status));
        if ( nSegs > 1 )
        {
            size_t off = pathing.getPathID().val()-1;
            HugeBVec::const_iterator itr(mPathSeq.begin(off));
            while ( --nSegs )
            {
                isRC = mGraph.nextEdge(*++itr,&edgeID,isRC);
                if ( edgeID.isNull() )
                    FatalErr("Invalid successor base in path sequence for read "
                             << readId );
                bool isPalindrome = mGraph.getEdge(edgeID).isPalindrome();
                CanonicalForm status = isPalindrome ? PALINDROME :
                                        isRC ? REV : FWD;
                result.push_back(EdgeDesc(edgeID,status));
            }
        }
    }
    return result;
}

template <unsigned K>
bvec& PathCollection<K>::getBases( EdgeList const& el, bvec* pBV ) const
{
    pBV->clear().reserve(getEdgeListLen(el));

    typedef HugeBVec::const_iterator Itr;
    for ( unsigned idx = 0; idx < el.size(); ++idx )
    {
        EdgeDesc const& ed = el[idx];
        UnipathEdge const& edge = mGraph.getEdge(ed.getEdgeID());
        Itr bItr = mGraph.getBases(edge.getInitialKmerID());
        Itr bEnd = mGraph.getBases(edge.getFinalKmerID()) + K;
        if ( ed.getStatus() != REV )
        {
            bItr += idx ? (K-1) : el.getInitialSkip();
            if ( idx+1 == el.size() ) bEnd -= el.getFinalSkip();
            pBV->append(bItr,bEnd);
        }
        else
        {
            bEnd -= idx ? (K-1) : el.getInitialSkip();
            if ( idx+1 == el.size() ) bItr += el.getFinalSkip();
            while ( bEnd != bItr )
                pBV->push_back(GetComplementaryBase(*--bEnd));
        }
    }
    return *pBV;
}

template <unsigned K>
void PathCollection<K>::validate( String const& fastbFilename )
{
    std::cout << Date() << " Validating pathings." << std::endl;
    VirtualMasterVec<bvec> vmv(fastbFilename.c_str());
    size_t nnn = vmv.size();
    bvec scratch;
    for ( size_t idx = 0; idx < nnn; ++idx )
    {
        EdgeList el(getEdgeList(idx));
        if ( el.size() )
            ForceAssertEq(vmv[idx],getBases(el,&scratch));
        else
            ForceAssertLt(vmv[idx].size(),K);
    }
}

template <unsigned K>
class PathsWithEvidence : public PathCollection<K>
{
public:
    PathsWithEvidence( String const& pathInfoFilename,
                       String const& graphInfoFilename );

    // copying prohibited by base class.  compiler-supplied destructor is OK

    size_t getNEdges() const { return mEvidence.size(); }

    UnipathEvidenceVec const& getEvidence( EdgeID const& edgeID ) const
    { return mEvidence[edgeID.val()]; }

private:
    VecUnipathEvidenceVec mEvidence; // indexed by EdgeID
};

template <unsigned K>
PathsWithEvidence<K>::PathsWithEvidence( String const& pathInfoFilename,
                                        String const& graphInfoFilename )
: PathCollection<K>(pathInfoFilename,graphInfoFilename)
{
    String evidenceFile =
            PathCollection<K>::getEvidenceFilename(pathInfoFilename);
    BinaryReader::readFile(evidenceFile.c_str(),&mEvidence);
    size_t nEdges = PathCollection<K>::getGraph().getNEdges();
    if ( mEvidence.size() != nEdges )
        FatalErr("Length of evidence vector (" << mEvidence.size() <<
                 ") doesn't jibe with the number of edges in the graph (" <<
                 nEdges << ").");
}

#endif /* KMERS_READ_PATHER_H_ */
