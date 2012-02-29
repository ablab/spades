///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "Basevector.h"

/* PerfectAlignerLG
 *
 * A class to generate all the perfect alignments between sequences in some
 * vecbasevector.
 *
 * PerfectAlignerLG::Align() is the main function.  It returns all perfect matches
 * of length >= K.  The aligns can also be filtered, as follows:
 *
 * behavior = PerfectAlignerLG::findImproper:
 *     No filtering.
 * behavior = PerfectAlignerLG::findSemiproper:
 *     One end of the alignment must coincide with the end of a target sequence.
 * behavior = PerfectAlignerLG::findProperOnly:
 *     Both ends of the alignment must coincide with the ends of target sequences.
 *
 *
 * This class is adapted from the very old PerfectAligner class.  It is better
 * for the following reasons:
 * -- PerfectAlignerLG uses the ParcelKmers paradigm instead of the older
 *    SortKmers.  ParcelKmers is fast, parallelizable, and safe for very large
 *    datasets ( sequences.size() > 2^31 )
 * -- PerfectAlignerLG returns stably-sorted results.
 *
 *
 * Josh Burton
 * January 2010
 *
 *********************************************************************************/



class PerfectAlignerLG {
public:
  
  enum Behavior {
    findProperOnly,
    findSemiproper,
    findImproper
  };
  
  PerfectAlignerLG( int K, 
		    Behavior behavior,
		    ostream* pLog = 0 )
    : m_K( K ), m_behavior( behavior ), m_pLog( pLog ) {}
  ~PerfectAlignerLG() {}
  
  void SetK( const int K ) { m_K = K; }
  void SetBehavior( const Behavior behavior ) { m_behavior = behavior; }
  
  /* Calculate all the perfect alignments in <sequences>.  Output goes to <aligns>.
   * 
   * n_processes:
   *      Parallelization level.  Used in ParcelKmers.
   * partition:
   *      If partition < 0, compare all sequences to all sequences.
   *      Otherwise, compare all the sequences with ids less than <partition> to
   *      all the sequences with ids greater than or equal to <partition>.
   * max_kmer_freq:
   *      Ignore K-mer matches that appear more than <max_kmer_freq> times.
   *
   *******************************************************************************/
  void Align( const vecbasevector& sequences, 
              vec<alignment_plus>& aligns,
	      const size_t n_processes = 1,
              const int partition = -1,
	      const size_t max_kmer_freq = 0 );
  
private:
  int m_K;
  Behavior m_behavior;
  ostream* m_pLog;
  
};

