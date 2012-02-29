///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOOKUP__FIRST_LOOK_ALIGN_H
#define LOOKUP__FIRST_LOOK_ALIGN_H

#include "lookup/LookAlign.h"
#include "lookup/LookupTab.h"

/**
 * struct first_look_align
 *
 * A simplified alignment class representing a gapless alignment found
 * via various FirstLookupFinder* classes.
 */
struct first_look_align
{
  uint64_t query_ID;     // could have more than 2^32 queries
  Location target_loc;   // 8 bytes
  int n_mismatches;      // sign indicates orientation of alignment
  

public:
  uint64_t get_query_ID() const { return query_ID; }
  uint64_t get_target_ID() const { return target_loc.getContig(); }
  int      get_num_mismatches() const { return n_mismatches; }
  bool     is_FW() const { return n_mismatches >= 0; }
  bool     is_RC() const { return n_mismatches < 0; }

  uint64_t get_start_on_target(const size_t query_len, const size_t K) const 
  { 
    if (this->is_FW())     
      return target_loc.getOffset(); 

    else {
      uint64_t n_bases_preceding_kmer_end = this->target_loc.getOffset() + K;
      uint64_t alignment_len = min(query_len, n_bases_preceding_kmer_end);
      uint64_t offset = n_bases_preceding_kmer_end - alignment_len;
      return offset;
    }
  
  }
  
  void PrintReadableBrief( ostream& out, const size_t query_len, const size_t K ) const {
    out << query_ID << ( this->is_RC() ? "rc" : "fw" ) << " vs "
	<< this->get_target_ID() << ", " << (this->is_FW() ? n_mismatches : -n_mismatches + 1) << " mismatches/(of " << query_len << "), starting at " << this->get_start_on_target( query_len, K ) << " on target\n";
  }
  
  // Convert to a gapless LookAlign.
  void convert_to_look_align(look_align &la,
                             const BaseVecVec &query,
                             const BaseVecVec &target,
                             const size_t K) const
  {
    la.a.SetNblocks(1);
    la.a.SetGap(0, 0);
    la.query_id = this->query_ID;
    la.query_length = query[la.query_id].size();
    la.target_id = this->target_loc.getContig();
    la.target_length = target[la.target_id].size();
    
    // Determine alignment length, orientation, etc.
    // The logic here is split into fw/rc cases and mirrors alignF & alignR.
    if (this->n_mismatches >= 0) { // i.e., if this alignment is forward
      la.rc1 = False;
      la.mutations = this->n_mismatches;
      
      size_t nBasesFollowingKmerStart = la.target_length - this->target_loc.getOffset();
      size_t alignmentLen = min((size_t)la.query_length, nBasesFollowingKmerStart);
      la.a.SetStartOnQuery(0);
      la.a.SetStartOnTarget(this->target_loc.getOffset());
      la.a.SetLength(0, alignmentLen);
    }
    else { // i.e., if this alignment is RC
      la.rc1 = True;
      la.mutations = -this->n_mismatches - 1;
      
      size_t nBasesPrecedingKmerEnd = this->target_loc.getOffset() + K;
      size_t alignmentLen = min((size_t)la.query_length, nBasesPrecedingKmerEnd);
      size_t offset = nBasesPrecedingKmerEnd - alignmentLen;
      la.a.SetStartOnQuery(la.query_length - alignmentLen);
      la.a.SetStartOnTarget(offset);
      la.a.SetLength(0, alignmentLen);
    }
    
    ForceAssert(la.IsProper());
  }
  
};

#endif
