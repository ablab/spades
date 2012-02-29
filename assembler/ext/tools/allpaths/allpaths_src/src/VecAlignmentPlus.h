///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef VEC_ALIGNMENT_PLUS
#define VEC_ALIGNMENT_PLUS

#include "Alignment.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "system/file/FileReader.h"
#include <cstddef>

/*
 * vec_alignment_plus
 *
 * It allows to view at alignments, keeping in memory only those for
 * which read_id1 < read_id2.
 *
 * Warning. The alignments must be sorted by read_id1.
 *
 * How to use.
 *
 * To fix the ideas, call all_aligns the vector of alignments you would
 * get the old way; and all_aligns_index the vector which returns the
 * first occurrence in all_aligns in which a given read_id is the first
 * read of the alignment (or -1). Then:
 *
 *  GetAlignment( al_plus, ii ); where al_plus = all_aligns[ii];
 *
 *  GetAlignsIndex( ii ) = all_aligns_index[ii].
 */
class vec_alignment_plus {

public:

  vec_alignment_plus( const String &alignments_file,
		      const vec<int> &read_lengths );
  
  vec_alignment_plus( const vec<alignment_plus> &orig_aligns,
		      const vec<int> &read_lengths );
  
  unsigned int GetNumberAlignments() const { return all_aligns_ids_.size(); }
  
  void GetAlignment( alignment_plus &al_plus, int align_id ) const;
  
  int GetAlignsIndex( int read_id ) const;

  int GetAlignmentId1( int align_id ) const;

  int GetAlignmentId2( int align_id ) const;

  float GetAlignmentScore( int align_id ) const;

  int GetAlignmentLength( int align_id, bool use_pos2 = false ) const;

  int GetSymmetricAlignmentLength( int align_id ) const;

  Bool Rc2( int align_id ) const;

  void SaveAlignments( const String &filename );
  
  // Set the alignment in all_aligns[ align_id ].
  void SetPlainAlignment( int align_id,
			  int pos1,
			  int pos2,
			  int errors,
			  const avector<int>& gaps,
			  const avector<int>& lengths,
			  int nblocks );
  
  // Set the alignment in all_aligns[ align_id ].
  void SetPlainAlignment( int align_id, alignment &plain_al );
  
  // Set the align in all_aligns[ align_id ] (beware: align != alignment).
  void SetPlainAlign( int align_id, align &plain_al );

  // Kill aligns that fail the RequireProper test (return killed count).
  int KillImproperAligns( );

  void SetAlignmentScore( int align_id, float score );
  

private:
  
  void Load( const String &alignments_file );
  
  
private:
  /*
   * Convention: call all_al the hypothetical vector with all the
   * alignments, sorted by rd_id1 ("all" means that if rd1<rd2 align,
   * then all_al contains two alignments, the one for rd1-rd2, and the
   * one for rd2-rd1).
   *
   * alignments_: all and only the entries of all_al with rd_id1<rd_id2;
   * all_aligns_ids_[ii] = jj, plus 
   * all_aligns_flip_[ii] = true (resp. false) means:
   *  all_al[ii] is equal to flipped (resp. not flipped) alignments_[jj];
   * all_aligns_index_[read_jj] = ii means:
   *  all_al[ii] is the first alignment with read_id1 = read_jj;
   */
  const vec<int> *read_lengths_;
  vec<alignment_plus> alignments_;
  mutable vec<Bool> all_aligns_flip_;  
  mutable vec<int> all_aligns_ids_;
  mutable vec<int> all_aligns_index_;

};

// Alignment indexing tools.  Note that these will not access
// preFindSeeds reads unless you modify the run_dir argument
// appropriately.

// BuildAlignsIndex: from all_aligns generate file aligns.index.  Then the
// AlignsIndexReader (below) can random-access this file to rapidly
// determine which reads a given read is known to align to.

void BuildAlignsIndex( const String& run_dir,
		       const vec_alignment_plus& all_aligns,
		       int n_reads );
void BuildAlignsIndex( const String& run_dir,
		       const vec<alignment_plus>& all_aligns,
		       int n_reads );

// Note that there is a main program util/BuildAlignmentsMain.cc which
// calls BuildAlignsIndex.  Ultimately, we should instead incorporate
// calls to BuildAlignsIndex into individual assembly executables.

class AlignsIndexReader
{
public:
    // assuming that BuildAlignsIndex has been run for this run_dir
    AlignsIndexReader( String const& run_dir, size_t nReads )
    : mFR( (run_dir + "/aligns.index").c_str() ), mNReads(nReads) {}

    // compiler-supplied destructor is OK
    // copying OK if mFR can be copied (which it currently can't be)

    // return as "to" all the read ids which read id1 is known to align to
    void readIndex( int id1, vec<int>& to ) const;

private:
    FileReader mFR;
    size_t mNReads;
};

#endif // VEC_ALIGNMENT_PLUS
