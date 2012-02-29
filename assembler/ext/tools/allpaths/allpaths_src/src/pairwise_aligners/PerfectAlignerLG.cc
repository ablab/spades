///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "pairwise_aligners/PerfectAlignerLG.h"

#include "kmers/KmerParcels.h"
#include "pairwise_aligners/MatchList.h"


// For general documentation, see PerfectAlignerLG.h
void PerfectAlignerLG::Align(const vecbasevector& sequences, 
                              vec<alignment_plus>& aligns,
                              const size_t n_processes,
                              const int partition,
                              const size_t max_kmer_freq)
{
  aligns.clear();
  
  ForceAssertLe(partition, (int)sequences.size());
  
  if (partition == 0 || partition == (int)sequences.size()) return;
  
  /// If partition greater than or equal to zero, only create aligns between
  /// sequences on either side of the partition index.
  bool using_partition = (partition >= 0);
  
  vec<int> seqIds(sequences.size());
  for (size_t i = 0; i < seqIds.size(); ++i) 
    seqIds[i] = i;
  
  MatchList perfectMatchList(sequences.size());
  
  // Our method of finding "matches" is to create a set of KmerParcels (i.e.,
  // a list of all kmers in the bases, sorted by their appearance) and then
  // look through the KmerBatches.  Each KmerBatch consists of a kmer and its
  // set of appearances in the bases; hence it implies one or more matches of
  // length K.
  
  // Create KmerParcels.
  if (m_pLog)
    *m_pLog << Date() << ": Creating KmerParcels..." << endl;
  
  // Build KmerParcels.  The parcels are stored in memory, not written to disk.
  KmerParcelsMemStore parcels(m_K);
  KmerParcelsBuilder builder(sequences, parcels, n_processes);
  builder.Build();
  size_t n_parcels = parcels.GetNumParcels();

  
  // Loop over all KmerBatches.  (Each KmerBatch appears in only one KmerParcel.)
  for (size_t parcel_ID = 0; parcel_ID < n_parcels; parcel_ID++) {
    
    KmerParcelReader parcel_reader(parcels, parcel_ID);
    
    while (parcel_reader.GetNextKmerBatch()) {

      const vec<KmerLoc> & kmer_locs = parcel_reader.CurrentKmerLocs();
      const size_t kmer_freq = kmer_locs.size();

      if (max_kmer_freq == 0 || kmer_freq <= max_kmer_freq) {
      
        // Loop over all pairs of kmer occurrences in this KmerBatch.
        for (size_t i = 0; i < kmer_freq; i++) {
          int64_t  id1 = kmer_locs[i].GetReadID();
          int64_t pos1 = kmer_locs[i].GetSignedPos();
          if (pos1 >= 0) pos1++; // required for MatchList
	
          for (size_t j = i+1; j < kmer_freq; j++) {
            int64_t  id2 = kmer_locs[j].GetReadID();
            int64_t pos2 = kmer_locs[j].GetSignedPos();
            if (pos2 >= 0) pos2++; // required for MatchList
	  
	  
            // Make sure these ids are on opposite sides of the partition, if using it.
            if (using_partition && ((id1 < partition) == (id2 < partition))) continue;
	  
            perfectMatchList.ProcessMatchingKmers(pos1, pos2, m_K, 
                                                  id1, id2, 
                                                  &sequences[id1], &sequences[id2]);
	  
          }
        } 
      }
    }
  }
  
  if (m_pLog)
    *m_pLog << Date() << ": Building alignments... " << endl;
  
  vec<Match> matches;
  
  for (size_t id = 0; id < sequences.size(); ++id) {
    perfectMatchList.GetMatches(id, matches);
    Sort(matches);
    
    for (unsigned int m = 0; m < matches.size(); ++m) {
      int64_t id1 = id;
      int64_t id2 = matches[m].GetId2();
      bool rc1 = false;
      bool rc2 = matches[m].GetRc();
      int64_t pos1 = matches[m].GetPos1();
      int64_t pos2 = matches[m].GetPos2();
      int64_t len  = matches[m].GetLen();
      int64_t len1 = sequences[id1].size();
      int64_t len2 = sequences[id2].size();
      
      // Count the number of ends reached.
      int n_ends = 0;
      if (pos1 == 0 || pos2 == 0) n_ends++;
      if (pos1 + len == len1 || pos2 + len == len2) n_ends++;
      
      // Apply the matching filter (which may be one of many behaviors)
      // to determine whether this alignment is acceptable.
      
      // Proper alignment: Two ends must be reached.
      if (m_behavior == PerfectAlignerLG::findProperOnly) {
	if (n_ends < 2) continue;
      }
      // Semiproper alignment: One end must be reached.
      else if (m_behavior == PerfectAlignerLG::findSemiproper) {
	if (n_ends < 1) continue;
      }
      // Improper alignment: All alignments accepted.
      else if (m_behavior == PerfectAlignerLG::findImproper) {}
      
      
      if (using_partition) {
	// Adjust id1, id2 so that they indicate offsets from the partition.
	if (id1 >= partition) {
	  id1 -= partition;
	  swap(id1, id2);
	  swap(len1, len2);
	  swap(rc1, rc2);
	  swap(pos1, pos2);
	}
	else {
	  ForceAssertGe(id2, partition);
	  id2 -= partition;
	}
      }
      
      // Create the alignment object and add it to the list.
      avector<int> gaps(1), lens(1);
      gaps(0) = 0;
      lens(0) = len;
      int errs = 0;
      alignment theAlignment(pos1, pos2, errs, gaps, lens);
      if (rc1) {
	theAlignment.ReverseThis(len1, len2);
	rc1 = !rc1;
	rc2 = !rc2;
      }
      
      alignment_plus theAlignmentPlus(id1, id2, len1, len2,
				       rc2, theAlignment, Float(0.0));
      
      aligns.push_back(theAlignmentPlus);
      
      if (!using_partition && id1 != id2) {
	theAlignmentPlus.Swap(len1, len2);
	aligns.push_back(theAlignmentPlus);
      }
    }
  }
  
  if (m_pLog)
    *m_pLog << Date() << ": Done building alignments." << endl;
}
