/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// CloseUnipathGaps.
//
// Purpose: Align fragment reads to the unipath graph with the specific goal of 
// aiding in the closure of "gaps" in the unipath graph.  Therefore we have to 
// allow for alignments that go off the end of a unipath.  For computational 
// efficiency we may also decide to keep only the alignments that are relevant 
// to the problem at hand.
//
// Seeding: We seed on perfect 20-mer matches, using the Kmer Parcel paradigm.
// We use only the first 20-mer in each read.
//
// Extension: We allow all possible paths through the unipath graph, including 
// those that "go off an end".  But see below: some searches will be truncated, 
// and the final output is restricted.  Indels are not allowed.
//
// Scoring: A placement of a read of length n is assigned scores s1,..., sn, where 
// si is the sum of the read quality scores at mismatches, up to position i.
//
// Assessment: Partial placements can be rejected according to the following rules:
// (a) if si exceeds a fixed upper bound (say 200);
// (b) if there is another placement s' such that s'i is small relative to si
// (say s'i < si - 100).
//
// Uniqueness: For a given read, we require that only one alignment passes.
//
// Output: We only keep the alignments that start on a terminal unipath and come 
// within 50 bases of its end.
//
// Note: now this does TWO passes.  One the second pass the reads are reversed.

// OpenMP requirements (necessary for the #pragma omp, below)
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include <omp.h>


#include "paths/CloseUnipathGapsCore.h"
#include "VecTemplate.h"
#include "kmers/KmerParcels.h"      // KmerParcelsCoupledBuilder
#include "kmers/KmerParcelsTools.h" // GetNextMatchingKmerBatches
#include "math/Functions.h"






class partial_placement {
private:

  vec<int> score; // score at base K + i
  vec<int> uni;   // sequence of unipaths
  int pos1;       // position on first unipath


public:

  partial_placement() { }

  partial_placement(const int u1, const int pos1)
    : pos1(pos1)
  {    uni.push_back(u1);    }
  
  // Pre-reserve memory for the score vector (since we know how big it is
  // likely to get)
  void ReserveScore(const int N) { score.reserve(N); }

  void AddAgree() { score.push_back(LastScore()); }

  void AddDisagree(const int q) { score.push_back(LastScore() + q); }

  int Score(const int i) const { return score[i]; }

  size_t NScores() const { return score.size(); }

  int LastScore() const
  {    return (score.empty() ? 0 : score.back());    }

  void AddUnipath(const int u) {
    uni.push_back(u);
  }

  Bool IsFirstUni() const { return uni.solo(); }

  int LastUni() const { return uni.back(); }

  int Uni(const int i) const { return uni[i]; }

  int NUni() const { return uni.size(); }

  int FirstPos() const { return pos1; }

  friend ostream& operator<<(ostream& out, const partial_placement& p)
  {    out << p.Uni(0) << "." << p.FirstPos();
  for (int j = 1; j < p.NUni(); j++)
    out << ", " << p.Uni(j);
  return out << ": " << p.NScores()
             << " scores recorded, last score = " 
             << p.LastScore();    }

};





void FindUnipathGapSeeds(const vecbasevector & bases, 
                         const vecbasevector & unibases,
                         const size_t K, 
                         const size_t n_threads, 
                         const String & work_dir, 
                         VecULongVec * seeds)
{
  
  // Create coupled KmerParcel files for the reference unibases and the read bases.
  // We combine the datasets with GetMatchingKmerParcelBatches
  cout << Date() << ": Creating KmerParcel files for FindUnipathGapSeeds" << endl;
  KmerParcelsDiskStore parcels_unibases(K, work_dir + "/FindUnipathGapSeeds.unibases");
  KmerParcelsDiskStore parcels_bases_trunc(K, work_dir + "/FindUnipathGapSeeds.bases_trunc");

  // Create truncated version of the reads.  Idiotic.
  size_t n_reads = bases.size();
  cout << Date() << ": n_reads = " << n_reads << endl;
  {
    vecbasevector bases_trunc(bases);
    for (size_t i = 0; i < n_reads; i++)
      if (bases_trunc[i].size() >= K) bases_trunc[i].resize(K);
    
    
    // build identical number of parcels for unibases and bases_trunc
    KmerParcelsCoupledBuilder builder(unibases, bases_trunc, 
				      parcels_unibases, parcels_bases_trunc,
				      n_threads);
    builder.Build();
  }


  
  // We will ignore any read that has more than 20 seed kmers; these
  // reads are probably repetitive.  This cutoff is necessary to prevent
  // the runtime from exploding.
  const size_t n_seeds_max = 20;
  vec<uint8_t> n_seeds(n_reads, 0);
  
  cout << Date() << ": Creating seeds..." << endl;
  // We will ignore any kmer that appears more than 20 times in the unibases;
  // these kmers are probably repetitive.  (N.B.: all non-palindromic kmers
  // appear more than once in the unibases, because the unibases contain fw/rc.)
  const size_t kmer_freq_max = 20;

  // Find seeds.
  // Loop over all KmerParcels, and parallelize over the loops.

  cout << Date() << ": compute seeds[query_ID].size() for each query_ID." << endl;
  
  seeds->clear();
  seeds->resize(n_reads);

  size_t n_parcels = parcels_unibases.GetNumParcels();

  // ---- pass 1: 
  //      Figure out the amount of space needed.
  //      We need this step to know how much memory to reserve for 
  //      each seeds[query_ID].  
  //      If we just do seeds[query_ID].push_back(seed), we impact dramatically 
  //      both time and memory because of repeated reallocations.
  //              
#pragma omp parallel for
  for (size_t parcel_ID = 0; parcel_ID < n_parcels; parcel_ID++) {
    
    // Open the KmerParcels for reference and reads.
    KmerParcelReader reader_unibases   (parcels_unibases,    parcel_ID);
    KmerParcelReader reader_bases_trunc(parcels_bases_trunc, parcel_ID);

    while (GetNextMatchingKmerBatches(reader_unibases, reader_bases_trunc)) {
      
      const vec<KmerLoc> & kmer_locs_unibases    = reader_unibases.CurrentKmerLocs();
      const vec<KmerLoc> & kmer_locs_bases_trunc = reader_bases_trunc.CurrentKmerLocs();
      const size_t n_kmers_unibases    = kmer_locs_unibases.size();
      const size_t n_kmers_bases_trunc = kmer_locs_bases_trunc.size();
      
      // keep going if this is NOT a repetitive kmer in the unibases
      if (n_kmers_unibases <= kmer_freq_max) {
      
        // Loop over every instance of this kmer in the bases_trunc 
        // and over every instance in the unibases.
        for (size_t i = 0; i != n_kmers_bases_trunc; i++) {

          const size_t query_ID = kmer_locs_bases_trunc[i].GetReadID();
          if (n_seeds[query_ID] < n_seeds_max) {
            
            for (size_t j = 0; j != n_kmers_unibases; j++) {
              
              if (kmer_locs_unibases[j].IsRC() == kmer_locs_bases_trunc[i].IsRC() && 
                  n_seeds[query_ID] < n_seeds_max)

                  n_seeds[query_ID]++;

            }

          }
        }

      }
    }
  }

  // ---- allocate space for seeds
  cout << Date() << ": reserve space for seeds[query_ID]." << endl;
  for (size_t i = 0; i != n_reads; i++)
    if (n_seeds[i] < n_seeds_max)
      (*seeds)[i].reserve(n_seeds[i]);


  // ---- pass 2: actually store the seeds.
  cout << Date() << ": store seeds in seeds[query_ID]." << endl;
  size_t n_seeds_total = 0;
#pragma omp parallel for
  for (size_t parcel_ID = 0; parcel_ID < n_parcels; parcel_ID++) {
    
    // Open the KmerParcels for reference and reads.
    KmerParcelReader reader_unibases   (parcels_unibases,    parcel_ID);
    KmerParcelReader reader_bases_trunc(parcels_bases_trunc, parcel_ID);

    while (GetNextMatchingKmerBatches(reader_unibases, reader_bases_trunc)) {
      
      const vec<KmerLoc> & kmer_locs_unibases    = reader_unibases.CurrentKmerLocs();
      const vec<KmerLoc> & kmer_locs_bases_trunc = reader_bases_trunc.CurrentKmerLocs();
      const size_t n_kmers_unibases    = kmer_locs_unibases.size();
      const size_t n_kmers_bases_trunc = kmer_locs_bases_trunc.size();
      
      // keep going if this is NOT a repetitive kmer in the unibases
      if (n_kmers_unibases <= kmer_freq_max) {
      
        // Loop over every instance of this kmer in the bases_trunc 
        // and over every instance in the unibases.
        for (size_t i = 0; i != n_kmers_bases_trunc; i++) {

          const size_t query_ID = kmer_locs_bases_trunc[i].GetReadID();
          if (n_seeds[query_ID] < n_seeds_max) {
            
            for (size_t j = 0; j != n_kmers_unibases; j++) {
              
              if (kmer_locs_unibases[j].IsRC() == kmer_locs_bases_trunc[i].IsRC()) {
                
                if ((*seeds)[i].size() < n_seeds_max) {

                  uint64_t seed = ((kmer_locs_unibases[j].GetReadID() << 32) |
                                   kmer_locs_unibases[j].GetUnsignedPos());
                  
                  (*seeds)[query_ID].push_back(seed);
                  n_seeds_total++;
                }
                else {
                  cout << "huh?" << endl;
                }
              }
            }

          }
        }

      }
    }
  }

  cout << Date() << ": " << n_reads << " number of reads processed." << endl;
  cout << Date() << ": " << n_seeds_total << " seeds created." << endl;
}



// Explanation of "extenders" in the following function:
// - first = index of unipath on which read is placed
// - second = start position of read on unipath
// - third = identifier of read.

void CloseUnipathGapsCore( const vecbasevector & bases, 
			    const vecqualvector & quals,
			    const vecbasevector & unibases, 
			    const vec< vec<int> > & nexts, 
			    const vec<int> & to_rc, 
			    const size_t K, 
                             const size_t UNIBASES_K,
			    const VecULongVec & seeds,
			    const VecULongVec & seeds_rc,
			    vec< triple<int,int,longlong> > & extenders, 
			    const int VERBOSITY, ostream& log)
{
  // Go through two passes, one for each orientation of the reads.

  size_t n_reads = bases.size();
  extenders.clear();
  extenders.reserve(n_reads);
  vec<Bool> placed(n_reads, False);
  vec<partial_placement> places, places_active;

  for (int pass = 1; pass <= 2; pass++) {
    log << Date() << ": starting pass " << pass << endl;
    const VecULongVec& SEEDS = (pass == 1 ? seeds : seeds_rc);

    // Go through the reads one by one.

#pragma omp parallel for private(places, places_active)
    for (size_t id = 0; id < n_reads; id++) {
      if (pass == 2 && placed[id]) continue;

      const basevector& b = bases[id];
      const qualvector& q = quals[id];
      const size_t read_length = b.size();
      const size_t n_seeds = SEEDS[id].size();
      if ( read_length < K ) continue;
      
      // For each position on the read, after the first K positions, we
      // keep track of the best observed score thus far.

      vec<int> best(read_length - K, 1e9);

      // Define rules.

      const int max_qual_err = 200;
      const int max_qual_err_diff = 100;

      // Go through the seeds.

      places.clear();
      places.reserve(n_seeds);
	       
      for (size_t si = 0; si < n_seeds; si++) {
        
        // Find all places arising from the seed.

        places_active.clear();
        places_active.reserve(n_seeds);
        // Unpack the longlong into two ints.
        int tig = SEEDS[id][si] >> 32;
        int pos = SEEDS[id][si] & 0xffffffff;
        ForceAssertGe(tig, 0); // trying to catch bug
        places_active.push(tig, pos);
        while(places_active.nonempty()) {
          partial_placement p = places_active.back(); // override p def
          places_active.pop_back();
	  p.ReserveScore(read_length - K);
          int u = p.LastUni();
          ForceAssertGe(u, 0); // trying to catch bug
          const basevector& U = unibases[u];
          Bool reject = False;
          int start = (p.IsFirstUni() ? pos + K : UNIBASES_K - 1);
          for (size_t j = start; j < U.size(); j++) {
	    int read_pos = p.NScores() + K;
	    if (read_pos == (int)read_length) break;
	    
	    // On second pass, use the inverted basevector/qualvector.
	    if (pass == 2) read_pos = read_length - 1 - read_pos;
            if (U[j] == (pass == 1 ? b[read_pos] : GetComplementaryBase(b[read_pos]))) p.AddAgree();
            else p.AddDisagree(q[read_pos]);
	    
            if (p.LastScore() > max_qual_err || 
                 p.LastScore() > best[p.NScores() - 1] + max_qual_err_diff) {
              reject = True;
              break;
            }
            best[ p.NScores() - 1 ] = Min(best[p.NScores() - 1], p.LastScore());
          }
          if (reject) continue;
          if (p.NScores() + K == read_length || nexts[u].empty()) {
            places.push_back(p);
          }
          else {
            for (size_t j = 0; j < nexts[u].size(); j++) {
              partial_placement pn(p);
	      pn.ReserveScore(read_length - K);
              pn.AddUnipath(nexts[u][j]);
              places_active.push_back(pn);
            }
          }
        }
      }

      // Screen placements.

      vec<unsigned int> reject_places;
      reject_places.reserve(places.size());
      for (size_t j = 0; j < places.size(); j++) {
        const partial_placement& p = places[j];
        for (size_t l = 0; l < p.NScores(); l++)
          if (places[j].Score(l) > best[l] + max_qual_err_diff) {
	    reject_places.push_back(j);
	    break;
	  }
      }
      EraseTheseIndices(places, reject_places);
      Destroy(reject_places);

      // Eliminate some equivalent placements.

      vec<Bool> to_delete(places.size(), False);
      for (size_t i1 = 0; i1 < places.size(); i1++) {
	const partial_placement &p1 = places[i1];
        for (size_t i2 = 0; i2 < places.size(); i2++) {
          const partial_placement &p2 = places[i2];
          if (p1.LastUni() != p2.LastUni()) continue;
          if (p1.NUni() != 1 || p2.NUni() == 1) continue;
          if (p1.NScores() != p2.NScores()) continue;
          to_delete[i2] = True;
        }
      }
      EraseIf(places, to_delete);

      // Announce results if requested.

      if (VERBOSITY == 2) {
        log << "\nread " << id << (pass == 1 ? " fw" : " rc") << "\n";
        log << "found " << places.size() << " placements\n";
        for (size_t i = 0; i < places.size(); i++) {
          log << "[" << i+1 << "] " 
               << (pass == 1 ? "fw" : "rc") << " "
               << places[i] << "\n";
        }
      }
      
      // Reject placements that are not solo, in an appropriate sense.

      if (places.empty() || places.size() > 2) continue;
      if (places.size() == 2) {
        Bool OK = False;
        const partial_placement &p1 = places[0], &p2 = places[1];
        if (p1.NUni() != 1 || p2.NUni() != 1) continue;
        int u1 = p1.Uni(0), u2 = p2.Uni(0);
        if (nexts[u1].empty() && nexts[ to_rc[u2] ].empty())
          OK = True;
        if (nexts[u2].empty() && nexts[ to_rc[u1] ].empty())
          OK = True;
        if (!OK) continue;
      }

      // Only keep unique alignments that go off the end of a unipath.

      for (size_t j = 0; j < places.size(); j++) {
        const partial_placement& p = places[j];
        if (p.NUni() != 1) continue;
        const int min_to_end = 100;
        Bool close_to_end = False;
        int u = p.Uni(0);
        int nu = unibases[u].size();
        if (nexts[u].empty()) {
          int dist_to_end = nu - (p.FirstPos() + read_length);
          if (dist_to_end <= min_to_end) close_to_end = True;
        }
        int ru = to_rc[u];
        if (nexts[ru].empty()) {
          int dist_to_end = p.FirstPos();
          if (dist_to_end <= min_to_end) close_to_end = True;
        }

        if (!close_to_end) continue;
        if (VERBOSITY == 1) 
          log << "\nread " << id << "\n" << p << "\n";
        if (pass == 1) {
#pragma omp critical 
	  {
	    extenders.push(u, p.FirstPos(), id);    
	    placed[id] = True;
	  }
        }
        else {
#pragma omp critical
	  {
	    extenders.push(ru, nu - (p.FirstPos() + read_length), id);
	  }
        }
      }
    }
  }
  Sort(extenders);
}
