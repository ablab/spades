/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef __INCLUDE_paths_PseudoUnipathGapCloser_h
#define __INCLUDE_paths_PseudoUnipathGapCloser_h


#include "Basevector.h"
#include "StdMethods.h"
#include "CommonSemanticTypes.h"
#include "paths/KmerPath.h"

/**
   Class: gapsequence
  
   Represents the base sequence that can be used to bridge the gap between
   two unipaths. When the two unipaths overlap the gap size is zero or
   negative and there is no base sequence. When the gap size is positive the
   base sequence to fill the gap is stored in seq.

   A <gapcloser> is a set of <gapsequences> representing the possible variants
   for what's between two CN1 unipaths.
*/
class gapsequence {
public:
  gapsequence() {};
  
  gapsequence(nbases_t size) : size(size) { seq.resize(0); }

  gapsequence(nbases_t size, const basevector& seq) : size(size), seq(seq) {};

  gapsequence(const basevector& seq) : size(seq.size()), seq(seq) {};

  STD_METHODS2( gapsequence, size, seq );

  friend ostream& operator<<(ostream& out, const gapsequence& gs) {
    return out << gs.size << " [" << gs.seq.ToString() << "]";
  }
  
  /// gap size - zero or negative for overlapping unipaths
  nbases_t size;
  /// sequence of bases that fills the gap - empty if gap size is negative
  basevector seq;
};

/**
   Class: gapcloser

   Represents the junction between two copy-number-one unipaths:

   >           gap
   >    u1    closer   u2  
   > ---------xxxxxx-------
   >   CN1     CN?    CN1

   The two unipaths may overlap, be immediately adjacent, or be separated
   by a non-zero-length gap.  In the latter case, we represent the sequence
   between the two unipaths, that we were able to infer.  Currently we have two
   methods of inferring this sequence: 1) bridge the gap by reads that align
   to the two unipaths, and of all such bridging reads take the consensus as to
   each base in the gap; 2) find paths in kmer space from the end of one unipath
   to the beginning of the other.
   The gapcloser can represent *multiple* valid paths between the two unipaths and
   the paths may be of mixed size; this represents our uncertainty about what actually
   goes on in the gap between the two CN1 unipaths.
*/
class gapcloser {

public:
  gapcloser( ) { }
    
  gapcloser( const unipath_id_t uid1, const unipath_id_t uid2, const int gap )
    : uid1(uid1), uid2(uid2)
  { gaps.push(gap); }
  
  gapcloser( const unipath_id_t uid1, const unipath_id_t uid2, const int gap, 
	    const basevector& gapseq )
    : uid1(uid1), uid2(uid2) 
  { gaps.push(gap, gapseq); }
  
  gapcloser( const unipath_id_t uid1, const unipath_id_t uid2, 
	    const vecbasevector& gapseqs )
    : uid1(uid1), uid2(uid2) {
    this->addGaps(gapseqs);
  }

  gapcloser( const unipath_id_t uid1, const unipath_id_t uid2, 
	    const vec<gapsequence>& gapseqs )
    : uid1(uid1), uid2(uid2), gaps(gapseqs) {}
 
  void addGaps(const vecbasevector& gapseqs) {
    for (size_t i = 0; i < gapseqs.size(); ++i)
      gaps.push(gapseqs[i]);
  }

  void addGap(const basevector& gapseq) {
    gaps.push(gapseq);
  }
  
  void addGap(const int gap) {
    ForceAssertLe(gap, 0);
    gaps.push(gap);
  }
    
  STD_METHODS3( gapcloser, uid1, uid2, gaps );

  /// Returns the number of valid gap sequences for this closer
  inline int validSequenceCount() {
    return gaps.isize();
  }

  /// Returns true if there is only one gap sequence for this closer
  inline Bool solo() const {
    return gaps.solo();
  }

  /// Returns the mean size of the gap sequences for this closer
  int meanGapSize() const;

  /// Returns the id of the unipath to the left of the closer
  inline unipath_id_t getUid1() const {
    return uid1;
  }

  /// Returns the id of the unipath to the right of the closer
  inline unipath_id_t getUid2() const {
    return uid2;
  }
  
  /// Returns the first gap filling sequence
  inline basevector getFirstSeq() const {
    return gaps[0].seq;
  }

  /// Returns the size of the first gap sequence (may be zero or negative)
  inline int getFirstSize() const {
    return gaps[0].size;
  }

  const vec<gapsequence>& getGapSequences() const {
    return gaps;
  }

  friend ostream& operator<<(ostream& out, const gapcloser& gc);

private:
  ///   uid1 - <unipath id> of the unipath to the left of the closer
  ///   uid2 - <unipath id> of the unipath to the right of the closer
  unipath_id_t uid1, uid2;
  ///   gaps - vec of gaps sequences
  vec<gapsequence> gaps;

};  // class gapcloser

#endif
// #ifndef __INCLUDE_paths_PseudoUnipathGapCloser_h


