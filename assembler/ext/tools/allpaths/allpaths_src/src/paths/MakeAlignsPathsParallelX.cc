/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// For documentation, see MakeAlignsPathsParallelX.h

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "kmers/KmerRecord.h"
#include "pairwise_aligners/BalancedMutmerGraph.h"
#include "paths/KmerPath.h"
#include "system/Worklist.h"
#include "kmers/KmerParcels.h"

#include "paths/MakeAlignsPathsParallelX.h"

static inline 
String Tag(String S = "MAPPX") { return Date() + " (" + S + "): "; } 


class BitVecVecLocked : private LockedData
{
private:
  BitVecVec & _data;

public:
  BitVecVecLocked(BitVecVec & data) : _data(data) {}

  // copy constructor and = not defined because of LockedData 
  BitVecVecLocked(const BitVecVecLocked & a);
  BitVecVecLocked & operator=(const BitVecVecLocked & a);

  void SetLocked(const KmerLoc & kmer_loc, 
                 const bool val) 
  {
    Locker lock(*this);
    _data[kmer_loc.GetReadID()].Set(kmer_loc.GetUnsignedPos(), val);
  }
};






// PathingProcessor class.
// This class is designed for use with the Worklist parallelization paradigm,
// in the MakeAlignsPathsParallelX function.
template<int I>
class PathingProcessor
{
private:

  // -------- MEMBER VARIABLES
  // The pointers all refer to externally owned data structures.
  
  // Inputs
  const BaseVecVec       & _bases;
  const KmerParcelsStore & _parcels;
  // Outputs
  BitVecVecLocked        * _kmer_chosen;
  BMG<I>                 * _bmg;
    
  
public:
  
  PathingProcessor(const BaseVecVec       & bases, 
                   const KmerParcelsStore & parcels,
                   BitVecVecLocked        * kmer_chosen, 
                   BMG<I>                 * bmg)
    : _bases(bases),
      _parcels(parcels),
      _kmer_chosen(kmer_chosen),
      _bmg(bmg)
  {};
  
  
  // operator(): This is the function called in parallel by the Worklist.
  // When called, it will process the kmers in parcel parcel_ID and modify
  // kmer_chosen_ and bmg_.  Modifications are regulated via LOCK_.
  void operator() (const size_t parcel_ID) 
  {
    const String tag = "MAPPX_P[" + ToString(parcel_ID) + "]";

    // Create the KmerParcelReader with this parcel_ID.
    KmerParcelReader parcel_reader(_parcels, parcel_ID);

    vec<KmerLoc> kmer_locs;
    while (parcel_reader.GetNextKmerBatch()) {
      const KmerBatch & batch = parcel_reader.CurrentKmerBatch();
      const size_t n_kmers = batch.size();

      if (n_kmers > 0) {
        kmer_locs = batch;
        
        // ribeiro:  kmer_locs for palindromes need to be added in RC form too 
        if (batch.IsPalindrome()) {
          for (size_t i = 0; i < n_kmers; i++)
            kmer_locs.push_back(batch[i].GetRC());
        }
        
        // ribeiro:  we need to sort the list of kmer locs with 
        // ribeiro:  { id1 <=> id2 and abs(pos1) <=> abs(pos2) } 
        // ribeiro:  otherwise the asserts in ReadsToPathsCoreX fail... for some reason 
        sort(kmer_locs.begin(), kmer_locs.end(), KmerLoc::ReadIDUnsignedPosLT);
        
        // ribeiro:  pack all nodes to merge in a vec
        const int64_t read_id1 = kmer_locs[0].GetReadID();
        int               pos1 = kmer_locs[0].GetSignedPos();
        
        if (pos1 >= 0) pos1++; // increment because MergeKmer expects pos != 0
        
        
        for (size_t i = 1; i < kmer_locs.size(); i++) {
          const int64_t read_id2 = kmer_locs[i].GetReadID();
          int               pos2 = kmer_locs[i].GetSignedPos(); 
          
          if (pos2 >= 0) pos2++; // increment because MergeKmer expects pos != 0
          
          _bmg->MergeKmer(pos1, pos2, _parcels.GetK(),
                          read_id1, read_id2,
                          _bases[read_id1], _bases[read_id2],
                          True);
        }

        // add canonical kmer
        (*_kmer_chosen).SetLocked(kmer_locs[0], true);
      }
    }

    //cout << Tag(tag) << "end... sizeof(bmg) = " + ToString((*_bmg).SizeOf()/1024) + " KB" << endl;
  }
};




// Output: a BMG (Balanced Mutmer Graph) and a 'kmer_chosen'
// BitVecVec indicating which kmers are canonical.
template<size_t I>
void MakeAlignsPathsParallelX(const size_t       K,
                              const BaseVecVec & bases, 
                              BitVecVec        * kmer_chosen, 
                              BMG<I>           * bmg, 
                              const String     & PARCEL_HEAD,
                              const size_t       NUM_THREADS,
                              const String       CHECKPOINT_HEAD)
{
  // Check that the reads are not too big for the BalancedMutmerGraph.
  // This should succeed because ReadsToPathsCoreY has broken up the reads.
  for ( size_t id = 0; id < bases.size(); id++ )
    ForceAssertLe( bases[id].size(), static_cast<unsigned>(USHRT_MAX) ); // USHRT_MAX = 65535


  String bmg_fn = CHECKPOINT_HEAD + ".bmg";
  String kmer_chosen_fn = CHECKPOINT_HEAD + ".chosen";

  if (CHECKPOINT_HEAD != "" && 
      IsRegularFile(bmg_fn) &&
      IsRegularFile(kmer_chosen_fn)) {

    cout << Tag() << "CHECKPOINT: retrieving Balanced Mutmer Graph 'bmg'." << endl;
    bmg->FromFile(bmg_fn);
    cout << Tag() << "CHECKPOINT: retrieved " << bmg->size() << " BMG nodes." << endl;

    cout << Tag() << "CHECKPOINT: retrieving 'kmer_chosen'." << endl;
    kmer_chosen->ReadAll(kmer_chosen_fn);
    cout << Tag() << "CHECKPOINT: retrieved " << kmer_chosen->size() << " kmer_chosen entries." << endl;

  }
  else {
    size_t n_parcels = 0;
    
    // ribeiro 2010-05-10: This below is really ugly.  
    // The parcels should be passed as parameters to this function instead 
    // of declared here as pointers.  
    // But today I'm too lazy to do anything about it.
    
    // Create the KmerParcelsStore.  If a PARCEL_HEAD is given, we'll be storing
    // (or reusing) the parcels on disk; otherwise they'll remain in memory.
    KmerParcelsStore * parcels;
    bool parcels_on_disk = ( PARCEL_HEAD != "" );
    
    if ( parcels_on_disk ) {
      parcels = new KmerParcelsDiskStore( K, PARCEL_HEAD );
      
      if (CHECKPOINT_HEAD != "" && 
	  ( ( KmerParcelsDiskStore * ) parcels )->ExistOnDisk() ) {
	// Reuse KmerParcels from disk.
	cout << Tag() << "CHECKPOINT: Reusing previously computed kmer parcels." << endl;
      }
      else {
	// Create KmerParcels on disk.
	if (CHECKPOINT_HEAD == "")
	  ( ( KmerParcelsDiskStore * ) parcels )->SetKeepOnDisk(false);
	KmerParcelsBuilder builder(bases, *parcels, NUM_THREADS);
	builder.Build();
      }
    }
    else {
      // Create KmerParcels in memory.
      parcels = new KmerParcelsMemStore( K );
      KmerParcelsBuilder builder( bases, *parcels, NUM_THREADS );
      builder.Build();
    }
    
    n_parcels = parcels->GetNumParcels();
    
    
    // "Balanced" Mutmer Graph: This object takes up a LOT of memory,
    // especially as kmers are loaded into it!
    bmg->resize(bases.size());
    
    // PARALLEL PROCESSING
    // Place each pass into the worklist.
    // As soon as it gets the first item, the worklist will immediately begin
    // processing the function PathingProcessor::operator( ), which processes
    // each KmerParcel and the ensuing work.  The worklist
    // runs its functions in multi-threaded parallel.
    // Each parallelized pass affects the variables kmer_chosen and bmg (through
    // their pointers in the PathingProcessor classes.
    
    //cout << Tag() << "begin " << n_parcels << " parcels ("
    //     << NUM_THREADS << "-way parallelized)" << endl;
    
    //cout << Date() << ": MemUsage before parallel runs: " << MemUsage( )/1000.0 << "Mb" << endl;
    {
      BitVecVecLocked kmer_chosen_locked(*kmer_chosen);
      
      PathingProcessor<I> processor(bases, *parcels, &kmer_chosen_locked, bmg);
      Worklist< size_t, PathingProcessor<I> > worklist(processor, NUM_THREADS - 1);
      
      for (size_t parcel_ID = 0; parcel_ID < n_parcels; parcel_ID++)
        worklist.add(parcel_ID);
      
      // The Worklist destructor (called at the end of scope) 
      // doesn't return until it has done all its processing
    } 
    
    //cout << Tag() << "end " << n_parcels << " parcels ("
    //     << NUM_THREADS << "-way parallelized)" << endl;
    

    // ribeiro: REALLY important to destroy the parcel store.
    // memory builds up pretty fast if you don't.
    // But really, the parcel store should be declared outside this function
    // and this deletion would be unecessary.  
    // (see my comment before the declaration of parcels)
    delete parcels;
    



  
    if (CHECKPOINT_HEAD != "") {
      cout << Tag() << "CHECKPOINT: keeping kmer parcels." << endl;
      cout << Tag() << "CHECKPOINT: saving Balanced Mutmer Graph 'bmg'."<< endl;
      bmg->ToFile(bmg_fn);
      cout << Tag() << "CHECKPOINT: saved " << bmg->size() << " BMG nodes." << endl;

      cout << Tag() << "CHECKPOINT: saving 'kmer_chosen'." << endl;
      kmer_chosen->WriteAll(kmer_chosen_fn);
      cout << Tag() << "CHECKPOINT: saved " << kmer_chosen->size() << " kmer_chosen entries." << endl;
    }
    
    //cout << Date( ) << ": MemUsage after parallel runs: " << MemUsage( )/1000.0 << "Mb" << endl;
  }

}







template void MakeAlignsPathsParallelX<1>(const size_t, 
                                          const BaseVecVec &,
                                          BitVecVec *,
                                          BMG<1> *,
                                          const String &,
                                          const size_t,
                                          const String CHECKPOINT_HEAD);

template void MakeAlignsPathsParallelX<2>(const size_t, 
                                          const BaseVecVec &,
                                          BitVecVec *,
                                          BMG<2> *,
                                          const String &,
                                          const size_t,
                                          const String CHECKPOINT_HEAD);











