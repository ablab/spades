/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/AddSuperReads.h"

#include "TaskTimer.h"

#include "paths/ExtendUnipathSeqs.h"
#include "paths/KmerPathDatabase.h"
#include "paths/OrientedKmerPathId.h"
#include "paths/Unipath.h"
#include "paths/UnipathSeq.h"
#include "paths/UnipathSeqBuilder.h"
#include "paths/UnipathSeqDatabase.h"

#include <queue>

// class SuperSeqBuilder

class SuperSeqBuilder {
 public:
  // Given the unipath sequences in unipathSeqsFw and unipathSeqsRc, 

  void Build( const vecUnipathSeq& unipathSeqsFw,
              const vecUnipathSeq& unipathSeqsRc,
              vecUnipathSeq& superSeqs,
              const vec<Mux>& muxesFw,
              const vec<Mux>& muxesRc,
              MuxGraph& muxGraph,
              SubsumptionList& subList ) const;
};


void FindSuperSeqMuxes( const vecKmerPath& unipaths, 
                        const vecUnipathSeq& superSeqs,
                        MuxGraph& muxGraph,
                        MuxGraph& inverseMuxGraph );

longlong FindSuperSeqSubs( const vecKmerPath& unipaths, 
                           const vecUnipathSeq& superSeqs,
                           SubsumptionList& subList );


void CalculateMaxOffsets( const vecKmerPath& pathsFw,
                          const vecKmerPath& pathsRc,
                          const MuxGraph& muxGraph,
                          const vec<read_pairing>& pairs,
                          const Float sdMult,
                          const int maxSep,
                          const int K,
                          vec<int>& maxOffsets );


void PrintSuperDotFile( const String& dotfile,
                        const MuxGraph& muxGraph,
                        const int numReads,
                        const vecUnipathSeq& unipathSeqsFw,
                        const vecUnipathSeq& unipathSeqsRc );


void AddSuperReads( const vecKmerPath& pathsFw, 
                    const vecKmerPath& pathsRc, 
                    const vec<read_pairing>& pairs,
                    const int K,
                    const int maxSep,
                    const Float sdMult,
                    vecKmerPath& allPathsFw, 
                    vecKmerPath& allPathsRc,
                    MuxGraph& allMuxes, 
                    SubsumptionList& allSubs,
                    OffsetTracker* pTracker )
{
  vec<tagged_rpint> pathsRpints;
  cout << "Building read path database... " << flush;
  for ( size_t i = 0; i < pathsFw.size(); ++i )
    pathsFw[i].AppendToDatabase( pathsRpints, i );
  for ( size_t i = 0; i < pathsRc.size(); ++i )
    pathsRc[i].AppendToDatabase( pathsRpints, -i-1 );
  Prepare( pathsRpints );
  cout << "done." << endl;

  int numReads = pathsFw.size();

  cout << "Building unipaths... " << flush;
  vecKmerPath unipaths;
  vec<tagged_rpint> unipathRpints;
  Unipath( pathsFw, pathsRc, pathsRpints, unipaths, unipathRpints );
  KmerPathDatabase unipathDB( &unipathRpints );
  cout << "done." << endl;

  vecUnipathSeq unipathSeqsFw, unipathSeqsRc;
  vec<Mux> muxesFw, muxesRc;

  // Convert the kmer paths of the reads into unipath sequences,
  // tracking the left extension of pathsFw[i] in muxesFw[i] and
  // pathsRc[j] in muxesRc[j].
  UnipathSeqBuilder seqBuilder( &unipaths, &unipathDB );
  seqBuilder.Build( pathsFw, unipathSeqsFw, muxesFw );
  seqBuilder.Build( pathsRc, unipathSeqsRc, muxesRc );

  // Create "new" reads that are the unique unipath sequences from the
  // previous step (stored in unextendedSuperSeqs).  Track the left
  // extensions of unipathSeqsFw and unipathSeqsRc in origMuxGraph.
  vecUnipathSeq unextendedSuperSeqs;
  MuxGraph origMuxGraph( numReads );
  SubsumptionList origSubList( numReads );
  SuperSeqBuilder().Build( unipathSeqsFw, unipathSeqsRc, 
                           unextendedSuperSeqs,
                           muxesFw, muxesRc,
                           origMuxGraph, origSubList );

  //   for ( int i = 0; i < unipathSeqsFw.size(); ++i )
  //     cout << i << "fw" << "\t" << unipathSeqsFw[i] << endl;
  //   cout << endl;
  //   for ( int i = 0; i < unipathSeqsFw.size(); ++i )
  //     cout << i << "rc" << "\t" << unipathSeqsRc[i] << endl;
  //   cout << endl;
  //   for ( int i = 0; i < unextendedSuperSeqs.size(); ++i )
  //     cout << i << "unext" << "\t" << unextendedSuperSeqs[i] << endl;
  //   cout << endl;

  // Extend the "new" reads where such extension is unambiguous.
  vecUnipathSeq extendedSuperSeqs;
  vec<Mux> superMuxes;
  ExtendUnipathSeqs( unipaths, unextendedSuperSeqs, extendedSuperSeqs, superMuxes );

  //   for ( int i = 0; i < extendedSuperSeqs.size(); ++i )
  //     cout << i << "ext" << "\t" << extendedSuperSeqs[i] << endl;
  //   cout << endl;

  // Some of these extended reads may now be redundant ( both A.B and
  // B.C might extend to A.B.C), so we condense down again by finding
  // the unique super sequences, which are stored in superSeqs.
  vecUnipathSeq superSeqs, bogusSeqs;
  vec<Mux> bogusMuxes;
  MuxGraph unextendedSuperMuxGraph( unextendedSuperSeqs.size() );
  SubsumptionList unextendedSuperSubList( unextendedSuperSeqs.size() );
  SuperSeqBuilder().Build( extendedSuperSeqs, bogusSeqs,
                           superSeqs,
                           superMuxes, bogusMuxes,
                           unextendedSuperMuxGraph, unextendedSuperSubList );

  int numSuperSeqs = superSeqs.size();
  cout << "Found " << numSuperSeqs << " unique unipath seqs." << endl;

  vecUnipathSeq allSeqsFw,allSeqsRc;
  allSeqsFw.reserve( numReads + numSuperSeqs );
  allSeqsRc.reserve( numReads + numSuperSeqs );

  vecKmerPath superPaths;
  ConvertUnipathSeqsToKmerPaths( superSeqs, unipaths, superPaths );

  allPathsFw.clear();
  allPathsRc.clear();

  allPathsFw.reserve( numReads + numSuperSeqs );
  allPathsRc.reserve( numReads + numSuperSeqs );

  allMuxes.resize( numReads + numSuperSeqs );
  allSubs.resize( numReads + numSuperSeqs );

  // cout << "Finding subsumptions among superseqs... " << flush;
  SubsumptionList superSubList( numSuperSeqs );
  longlong numSubs = FindSuperSeqSubs( unipaths, superSeqs, superSubList );
  // cout << "done." << endl;
  
  // cout << "Found " << numSubs << " subsumptions." << endl;

  // We now want to add in all the muxes and subsumptions that map
  // from the original reads to these unique extended super reads.
  for ( int i = 0; i < numReads; ++i ) {
    allSeqsFw.push_back( unipathSeqsFw[i] );
    allSeqsRc.push_back( unipathSeqsRc[i] );

    allPathsFw.push_back( pathsFw[i] );
    allPathsRc.push_back( pathsRc[i] );

    for ( int rc = 0; rc < 2; ++rc ) {
      OrientedKmerPathId okpid( i, (rc==1) );

      if ( okpid.GetPathPtr( allPathsFw, allPathsRc )->IsEmpty() )
        continue;

      int numKmers = 0;
      int segment = 0;

      // First, we look up the mux of each read to the
      // unextended unipath sequence it is subsumed by.
      vec<Mux> muxes;
      origMuxGraph.GetMuxesOf( okpid, muxes );
      Mux origMux = muxes.front();
      numKmers += origMux.GetNumKmers();
      segment += origMux.GetSegment();
      
      // We then look up the mux of that unextended unipath sequence
      // to the unique extended unipath sequence it is subsumed by.
      unextendedSuperMuxGraph.GetMuxesOf( origMux.GetPathId(), muxes );
      Mux superMux = muxes.front();
      numKmers += superMux.GetNumKmers();
      segment += superMux.GetSegment();

      int superId = numReads + superMux.GetPathId().GetId();

      allMuxes.SetMuxOf( okpid, Mux( OrientedKmerPathId( superId, false ), segment, numKmers ) );

      int leftOverhang = 0;
      
      vec<SubsumptionRecord> subRecs;

      origSubList.GetFullRecordsFor( okpid, subRecs );
      SubsumptionRecord origSubRec = subRecs.front();
      leftOverhang += origSubRec.GetLeftOverhang();

      unextendedSuperSubList.GetFullRecordsFor( origSubRec.GetSuperPathId(), subRecs );
      SubsumptionRecord unextendedSuperSubRec = subRecs.front();
      leftOverhang += unextendedSuperSubRec.GetLeftOverhang();
      
      ForceAssertEq( unextendedSuperSubRec.GetSuperPathId().GetId() + numReads, superId );

      vec<BriefSubsumptionRecord> newSubRecs;
      
      newSubRecs.push_back( BriefSubsumptionRecord( OrientedKmerPathId( superId, false ),
                                                    leftOverhang ) );

      superSubList.GetFullRecordsFor( unextendedSuperSubRec.GetSuperPathId(), subRecs );
      for ( unsigned int r = 0; r < subRecs.size(); ++r ) {
        OrientedKmerPathId superPathId( numReads + subRecs[r].GetSuperPathId().GetId(), false );
        int transitiveLeftOverhang = leftOverhang + subRecs[r].GetLeftOverhang();
        newSubRecs.push_back( BriefSubsumptionRecord( superPathId,
                                                      transitiveLeftOverhang ) );
      }
      allSubs.SetBriefRecordsFor( okpid, newSubRecs );
    }
  }

  // cout << "Find muxes among superseqs... " << flush;
  MuxGraph superMuxGraph( numSuperSeqs );
  MuxGraph superInverseMuxGraph( numSuperSeqs );
  FindSuperSeqMuxes( unipaths, superSeqs, superMuxGraph, superInverseMuxGraph );
  // cout << "done." << endl;

  //PrintSuperDotFile( "tmp/super.dot", superMuxGraph, numReads, unipathSeqsFw, unipathSeqsRc );

  MuxGraph allInverseMuxes( numReads + numSuperSeqs );
    
  for ( int i = 0; i < numSuperSeqs; ++i ) {
    allSeqsFw.push_back( superSeqs[i] );
    allSeqsRc.push_back( UnipathSeq() );
    
    allPathsFw.push_back( superPaths[i] );
    allPathsRc.push_back( KmerPath() );
    
    vec<Mux> muxes;
    superMuxGraph.GetMuxesOf( OrientedKmerPathId( i, false ), muxes );
    
    for ( unsigned int m = 0; m < muxes.size(); ++m )
      muxes[m].SetPathId( OrientedKmerPathId( numReads + muxes[m].GetPathId().GetId(), false ) );
    
    allMuxes.SetMuxesOf( OrientedKmerPathId( numReads + i, false ), muxes );
    
      
    superInverseMuxGraph.GetMuxesOf( OrientedKmerPathId( i, false ), muxes );
    
    for ( unsigned int m = 0; m < muxes.size(); ++m )
      muxes[m].SetPathId( OrientedKmerPathId( numReads + muxes[m].GetPathId().GetId(), false ) );
    
    allInverseMuxes.SetMuxesOf( OrientedKmerPathId( numReads + i, false ), muxes );
  }
  
  if ( pTracker ) {
    vec<int> maxOffsets;
    CalculateMaxOffsets( pathsFw, pathsRc,
                         allMuxes, pairs,
                         sdMult, maxSep, K,
                         maxOffsets );

    TaskTimer offsetsTimer;
    offsetsTimer.Start();
    pTracker->ConvertFrom( MutableOffsetTracker( allSeqsFw, 
                                                 allInverseMuxes, 
                                                 numReads, 
                                                 maxOffsets ) );
    offsetsTimer.Stop();
    
    PRINT( offsetsTimer );
  }
}



// SuperSeqBuilder methods.

struct UnipathSeqPtr_Okpid {
  UnipathSeqPtr_Okpid( const UnipathSeq* pS, const OrientedKmerPathId& o )
    : pSeq( pS ), okpid( o ) {}

  const UnipathSeq* pSeq;
  OrientedKmerPathId okpid;
};

inline
bool operator< ( const UnipathSeqPtr_Okpid& lhs,
                 const UnipathSeqPtr_Okpid& rhs ) {
  return ( *(lhs.pSeq) < *(rhs.pSeq) );
}


// Presume that each of the UnipathSeqs in unipathSeqs{Fw,Rc} was
// created from KmerPaths offset by the Muxes stored in muxes{Fw,Rc}.
// Find the unique UnipathSeqs in unipathSeqs{Fw,Rc} and put them in
// superSeqs and fill out muxGraph and subList mapping the relation
// between the KmerPaths from which unipathSeqs{Fw,Rc} were made and
// the unique UnipathSeqs.
void SuperSeqBuilder::Build( const vecUnipathSeq& unipathSeqsFw,
                             const vecUnipathSeq& unipathSeqsRc,
                             vecUnipathSeq& superSeqs,
                             const vec<Mux>& muxesFw,
                             const vec<Mux>& muxesRc,
                             MuxGraph& muxGraph,
                             SubsumptionList& subList ) const {
  vec<UnipathSeqPtr_Okpid> seqPtrs;

  for ( vecUnipathSeq::size_type i = 0; i < unipathSeqsFw.size(); ++i )
    if ( unipathSeqsFw[i].size() )
      seqPtrs.push_back( UnipathSeqPtr_Okpid( &(unipathSeqsFw[i]), 
                                              OrientedKmerPathId( i, false ) ) );

  for ( vecUnipathSeq::size_type i = 0; i < unipathSeqsRc.size(); ++i )
    if ( unipathSeqsRc[i].size() )
      seqPtrs.push_back( UnipathSeqPtr_Okpid( &(unipathSeqsRc[i]), 
                                              OrientedKmerPathId( i, true ) ) );

  // cout << "Sorting unipath seqs... " << flush;
  sort( seqPtrs.begin(), seqPtrs.end() );
  // cout << "done." << endl;
 
  vec<UnipathSeqPtr_Okpid>::iterator rangeBegin, rangeEnd;

  // cout << "Finding unique unipath seqs... " << flush;
  for ( rangeBegin = seqPtrs.begin(); 
        rangeBegin != seqPtrs.end();
        rangeBegin = rangeEnd ) {
    rangeEnd = upper_bound( rangeBegin, seqPtrs.end(), *rangeBegin );
    
    int superSeqId = superSeqs.size();
    superSeqs.push_back_reserve( *(rangeBegin->pSeq), 0, 2.0 );

    OrientedKmerPathId superOkpid( superSeqId, false );

    for ( ; rangeBegin != rangeEnd; ++rangeBegin ) {
      OrientedKmerPathId readOkpid( rangeBegin->okpid );
      int id = readOkpid.GetId();

      const Mux& theMux = ( readOkpid.IsFw() ? muxesFw[id] : muxesRc[id] );

      vec<BriefSubsumptionRecord> finalSubRecs;
      finalSubRecs.push_back( BriefSubsumptionRecord( superOkpid, theMux.GetNumKmers() ) );
      subList.SetBriefRecordsFor( readOkpid, finalSubRecs );
    
      muxGraph.SetMuxOf( readOkpid, Mux( superOkpid, theMux.GetSegment(), theMux.GetNumKmers() ) );
    }
  }
  // cout << "done." << endl;
}



// Is otherRec an extension of thisRec?
bool FirstExtendsSecond( const vecUnipathSeq& superSeqs,
                         const UnipathSeqDatabase::Record& first,
                         const UnipathSeqDatabase::Record& second, 
                         bool verbose ) {
  if ( verbose ) {
    cout << "FirstExtendsSecond( ";
    first.Print( cout, superSeqs );
    cout << ", ";
    second.Print( cout, superSeqs );
    cout << " ): " << flush;
  }

  int thisId = second.seqId;
  int thisIndex = second.index;
  
  int otherId = first.seqId;
  int otherIndex = first.index;

  if ( thisIndex > otherIndex ) {
    if ( verbose ) 
      cout << "no" << endl;
    return false;
  }

  int maxOverlapLeft = min( thisIndex, otherIndex );
  int overlapLeft = 0;
  for ( ; overlapLeft < maxOverlapLeft; ++overlapLeft )
    if ( superSeqs[thisId][thisIndex-overlapLeft-1] != 
         superSeqs[otherId][otherIndex-overlapLeft-1] ) {
      if ( verbose ) 
        cout << "no" << endl;
      return false;
    }

  int maxOverlapRight = min( superSeqs[thisId].size() - thisIndex,
                             superSeqs[otherId].size() - otherIndex );
  int overlapRight = 1;
  for ( ; overlapRight < maxOverlapRight; ++overlapRight )
    if ( superSeqs[thisId][thisIndex+overlapRight] != 
         superSeqs[otherId][otherIndex+overlapRight] ) {
      if ( verbose ) 
        cout << "no" << endl;
      return false;
    }

  // Deal with additional conditions on coterminal extensions.
  if ( thisIndex == otherIndex ) {

    // A superSeq cannot be a coterminal extension of itself.
    if ( otherId == thisId ) {
      if ( verbose ) 
        cout << "no" << endl;
      return false;
    }

    // A superSeq cannot have a coterminal extension that is longer than itself.
    if ( superSeqs[thisId].size() < superSeqs[otherId].size() ) {
      if ( verbose ) 
        cout << "no" << endl;
      return false;
    }
        
    // A coterminal extension should never be the same size as the
    // original, since they would then be identical.  SuperSeqs
    // are supposed to be unique.
    ForceAssert( superSeqs[thisId].size() != superSeqs[otherId].size() );
  }
  
  if ( verbose ) 
    cout << "yes" << endl;
  return true;
}


// FindSuperSeqMuxes

void FindSuperSeqMuxes( const vecKmerPath& unipaths,
                        const vecUnipathSeq& superSeqs,
                        MuxGraph& muxGraph,
                        MuxGraph& inverseMuxGraph ) {
  UnipathSeqDatabase superSeqDB( superSeqs );

  muxGraph.resize( superSeqs.size() );

  bool verbose = false;

  for ( vecUnipathSeq::size_type thisId = 0; thisId < superSeqs.size(); ++thisId ) {
    int thisIndex = 0;

    vec<UnipathSeqDatabase::Record> records, extensions;
    superSeqDB.Find( superSeqs[thisId][thisIndex], records );

    UnipathSeqDatabase::Record thisRecord( superSeqs[thisId][thisIndex], thisId, thisIndex );

    // Find all extensions.
    if ( verbose ) {
      cout << "Finding all extensions of ";
      thisRecord.Print( cout, superSeqs );
      cout << endl;
    }

    for ( unsigned int r = 0; r < records.size(); ++r )
      if ( FirstExtendsSecond( superSeqs, records[r], thisRecord, verbose ) )
        extensions.push_back( records[r] );

    if ( verbose ) {
      cout << "Eliminating non-minimal extensions of ";
      thisRecord.Print( cout, superSeqs );
      cout << endl;
    }

    // Find minimal extensions.
    vec<Bool> isMinimal( extensions.size(), True );

    for ( unsigned int i = 0; i < extensions.size(); ++i ) 
      for ( unsigned int j = i+1; j < extensions.size(); ++j ) {
        if ( FirstExtendsSecond( superSeqs, extensions[i], extensions[j], verbose ) )
          isMinimal[i] = False;
        else if ( FirstExtendsSecond( superSeqs, extensions[j], extensions[i], verbose ) ) 
          isMinimal[j] = False;
      }

    OrientedKmerPathId thisOkpid( thisId, false );

    vec<Mux> muxes, inverseMuxes;
    for ( unsigned int i = 0; i < extensions.size(); ++i ) {
      if ( ! isMinimal[i] )
        continue;

      int otherId = extensions[i].seqId;
      int index = extensions[i].index;

      int numKmers = 0;
      int segment = 0;
      for ( int i = 0; i < index; ++i ) {
        numKmers += unipaths[ superSeqs[otherId][i] ].KmerCount();
        segment += unipaths[ superSeqs[otherId][i] ].NSegments();
        if ( i > 0 && 
             unipaths[ superSeqs[otherId][i-1] ].End().GetKmer() + 1 ==
             unipaths[ superSeqs[otherId][i] ].Begin().GetKmer() )
          --segment;
      }
          
      OrientedKmerPathId otherOkpid( otherId, false );

      muxes.push_back( Mux( otherOkpid, segment, numKmers ) );
      
      inverseMuxGraph.GetMuxesOf( otherOkpid, inverseMuxes );
      inverseMuxes.push_back( Mux( thisOkpid, segment, numKmers ) );
      inverseMuxGraph.SetMuxesOf( otherOkpid, inverseMuxes );
    }
    muxGraph.SetMuxesOf( thisOkpid, muxes );
  }
}


longlong FindSuperSeqSubs( const vecKmerPath& unipaths, 
                           const vecUnipathSeq& superSeqs,
                           SubsumptionList& subList ) 
{
  UnipathSeqDatabase superSeqDB( superSeqs );

  subList.resize( superSeqs.size() );

  int subsFound = 0;

  for ( vecUnipathSeq::size_type thisId = 0; thisId < superSeqs.size(); ++thisId ) {
    int thisIndex = 0;

    vec<UnipathSeqDatabase::Record> records;
    vec<BriefSubsumptionRecord> subsumptions;
    superSeqDB.Find( superSeqs[thisId][thisIndex], records );

    for ( unsigned int r = 0; r < records.size(); ++r ) {
      vecUnipathSeq::size_type otherId = records[r].seqId;
      UnipathSeq::size_type otherIndex = records[r].index;

      if ( otherId == thisId ) 
        continue;
      
      if ( superSeqs[otherId].size() - otherIndex < superSeqs[thisId].size() )
        continue;

      UnipathSeq::size_type agreementLength = 1;
      for ( ; agreementLength < superSeqs[thisId].size(); ++agreementLength )
        if ( superSeqs[thisId][agreementLength] != superSeqs[otherId][otherIndex+agreementLength] )
          break;

      if ( agreementLength != superSeqs[thisId].size() )
        continue;

      int numKmers = 0;
      for ( UnipathSeq::size_type i = 0; i < otherIndex; ++i )
        numKmers += unipaths[ superSeqs[otherId][i] ].KmerCount();

      subsumptions.push_back( BriefSubsumptionRecord( OrientedKmerPathId( otherId, false ), 
                                                      numKmers ) );
      ++subsFound;
    }
    
    subList.SetBriefRecordsFor( OrientedKmerPathId( thisId, false ), subsumptions );
  }

  return subsFound;
}


void CalculateMaxOffsets( const vecKmerPath& pathsFw,
                          const vecKmerPath& pathsRc,
                          const MuxGraph& muxGraph,
                          const vec<read_pairing>& pairs,
                          const Float sdMult,
                          const int maxSep,
                          const int K,
                          vec<int>& maxOffsets ) {
  maxOffsets.clear();
  maxOffsets.resize( muxGraph.size(), 0 );
  
  for ( unsigned int p = 0; p < pairs.size(); ++p ) {
    int id1 = pairs[p].id1;
    int id2 = pairs[p].id2;
    int sep = pairs[p].sep;
    int dev = pairs[p].sd;

    if ( maxSep > 0 && sep > maxSep )
      continue;
    
    if ( pathsFw[id1].IsEmpty() || pathsRc[id2].IsEmpty() )
      continue;

    for ( int rc = 0; rc < 2; ++rc ) {
      OrientedKmerPathId fwOkpid( (rc?id2:id1), false );
      OrientedKmerPathId rcOkpid( (rc?id1:id2), true );

      int offset = fwOkpid.GetPathPtr( pathsFw, pathsRc )->KmerCount() + K-1 + sep;

      vec<Mux> muxes;
      muxGraph.GetMuxesOf( fwOkpid, muxes );
      Mux fwMux = muxes.front();
      offset += fwMux.GetNumKmers();
      int from = fwMux.GetPathId().GetId();

      muxGraph.GetMuxesOf( rcOkpid, muxes );
      Mux rcMux = muxes.front();
      offset -= rcMux.GetNumKmers();

      offset += (int) ceil( sdMult * Float( dev ) );

      if ( maxOffsets[from] < offset )
        maxOffsets[from] = offset;
    }
  }
}


void PrintSuperDotFile( const String& dotfile,
                        const MuxGraph& muxGraph,
                        const int numReads,
                        const vecUnipathSeq& unipathSeqsFw,
                        const vecUnipathSeq& unipathSeqsRc )
{
  ofstream dot( dotfile.c_str() );
  
  dot << "digraph G {" << endl;
  dot << "  node [shape=box];" << endl;
  dot << "  rankdir=LR;" << endl;
  dot << endl;

  vec<Mux> muxes;
  for ( int thisId = 0; thisId < muxGraph.size(); ++thisId ) {
    for ( int rc = 0; rc < 2; ++rc ) {
      OrientedKmerPathId thisOkpid( thisId, (rc==1) );
      muxGraph.GetMuxesOf( thisOkpid, muxes );
    
      const UnipathSeq& thisUnipathSeq = ( rc==1 ? unipathSeqsRc[thisId] : unipathSeqsFw[thisId] );
      for ( vec<Mux>::iterator iMux = muxes.begin(); iMux != muxes.end(); ++iMux ) {
        int otherId = iMux->GetPathId().GetId();
        const UnipathSeq& otherUnipathSeq = ( iMux->GetPathId().IsFw() ?
                                              unipathSeqsFw[otherId] :
                                              unipathSeqsRc[otherId] );
        if ( otherId < numReads ) 
          dot << '"' << iMux->GetPathId() << '"';
        else 
          dot << '"' << otherUnipathSeq << '"';
        dot << " -> ";
        if ( thisId < numReads )
          dot << '"' << thisOkpid << '"';
        else
          dot << '"' << thisUnipathSeq << '"';
        dot << " [label=\"" << iMux->GetNumKmers() << "\"];" << endl;
      }
    }
  }
                                         
  dot << "}" << endl;
}

