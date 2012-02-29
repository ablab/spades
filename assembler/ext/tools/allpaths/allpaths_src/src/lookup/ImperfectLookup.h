/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

//// ImperfectLookup:  Find good unique alignments with no indels quickly.
///
/// \fn ImperfectLookup
///
/// Given a vector "query" of sequences Q,
/// (a) find all gap-free alignments A of Q to target (defined by "lookup_file",
/// from MakeLookupTable with the given "K"), that subsume a perfect match of length
/// K and that are unique in the sense that there is no other such alignment B for
/// which errors(B) <= errors(A) + best_prox;
/// (b) find the minimum number of mismatches in a gap-free alignment of Q to target
/// that subsumes a perfect match of length K.
///
/// If no alignment of Q is found, the value for (b) is infinitely_many.
///
/// The variable AlignDir controls whether we look for forward alignments, or for
/// both forward and reverse alignments.
///
/// Note that this code will not work for finding indel polymorphisms!
///
/// Note that for certain applications, the code could be speeded up by keeping the
/// lookup table in memory, rather than reading it from disk on each call.
///
/// I'm not sure how this handles Ns in the target.
///
/// If unique_aligned is specified, do not modify aligns, and instead return a
/// boolean value as to whether each read aligns uniquely.
///
/// If quality scores are passed in, they are used to weight the errors.
/// Details of the algorithm are in the doc for CountMismatches in .cc file.
///
/// \sa PerfectLookup()
/// \ingroup grp_eval

#ifndef IMPERFECT_LOOKUP_H
#define IMPERFECT_LOOKUP_H

#include "Basevector.h"
#include "CoreTools.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"
#include "BasevectorTools.h"
#include "lookup/LookupTools.h"
#include "lookup/AlignCollector.h"
#include "solid/Solid.h"

class TaskTimer;


void ImperfectLookup( lookup_table &look, const vecbasevector& reads,
		      AlignCollectorBase & aligns,
		      AlignDir direction,
		      vec<TaskTimer *> * timers = 0,
		      const vecqualvector * quals = 0,
		      bool useReverse=false, const unsigned int max_freq=0);



/// Wrapper overload that reads the lookup file from disk and
/// passes execution to ImperfectLookup(lookup_table &, vecbasevector &, AlignCollector &,...)
/// (see docs for that method).
void ImperfectLookup(
     // inputs:
     const String & lookup_file, vecbasevector& reads,
     // inputs-outputs (the logic of accepting/rejecting aligns is passed in with the collector):
     AlignCollectorBase & aligns,
     // parameters:
     AlignDir direction,
     // other:
     vec<TaskTimer *> * timers = 0,
     const vecqualvector * quals = 0,
     bool useReverse = false, const unsigned int max_freq=0);


/// Overload that reads in the lookup file from disk and passes execution
/// to ImperfectLookup(lookup_table &, vecbasevector &, vec<look_align> &,...).
/// This is a legacy method with flat list of parameters which is too long and
/// with more rigid logic,
/// see more flexible overload templatized by AlignCollector.
void ImperfectLookup(
     // inputs:
     unsigned int K, const vecbasevector& query, const String& lookup_file,
     // outputs:
     vec<look_align>& aligns, vec<int>& min_errors,
     // parameters:
     AlignDir direction, int best_prox, int max_errors = -1,
     // other:
     vec<Bool>* unique_aligned = 0,
     vec<TaskTimer *> * timers = 0,
     const vecqualvector * quals = 0,
     bool useReverse = false, const unsigned int max_freq = 0);

/// Searches for gapless alignment with mismatches between each of the query
/// sequences (reads) and the reference passed as a lookup table object. This is a legacy
/// interface with long flat list of parameters and rigid logic. This overload
/// can only search for best "unique" alignments satisfying the following criteria:
/// 1) the best alignment has max_errors or fewer errors,
/// 2) the next best alignment has at least best_prox + 1 errors more than the
///    best alignment.
/// This method saves and returns only these best unique alignments (in \c aligns)
/// if unique_aligned==0. Otherwise, \c aligns is not populated and only vector
/// of boolean flags \c unique_aligned is populated (with true at position i if
/// the query reads[i] aligned uniquely, false otherwise. In both cases, the
/// number of errors in the best alignment for query i is returned in best_errors[i].
/// NOTE: the value returned in best_errors is correct only for uniquely aligned reads;
/// for the reads that do not align uniquely, this method gives up
/// early and reported error count may be greater than in the best possible alignment.
/// See ImperfectLookup(lookup_table &, vecbasevector &, AlignCollector &,...) for
/// more flexible implementation templatized by the "align acceptance" policy.
void ImperfectLookup(
     // inputs:
     lookup_table &look, const vecbasevector& reads,
     // outputs:
     vec<look_align>& aligns, vec<int>& best_errors,
     // parameters:
     AlignDir direction, int best_prox, int max_errors = -1,
     // other:
     vec<Bool>* unique_aligned = 0,
     vec<TaskTimer *> * timers = 0,
     const vecqualvector * quals = 0,
     bool useReverse = false, const unsigned int max_freq = 0 );


/// A "modern" overload of ImperfectLookup that takes align collector as an argument.
/// This method achieves better concept separation: the only job of the implementation
/// is to look for all putative ungapped alignments of the reads and pass them
/// to the collector. It is the job of the collector implementation to decide whether
/// the alignment is "good enough" (i.e. unique enough, has fewer mismatches than some
/// max cutoff etc); collector can also decide to accept a few alignments per read if needed.
///
/// The template parameter 'AlignCollector' must be a model of the type concept AlignCollector
/// defined in AlignCollector.h .
///
/// @param[in] look lookup table (reference genome) to generate alignments against
/// @param[in] reads sequences, for which the alignments are sought
/// @param[out] aligns Align collector; any previous content will be \em conserved up to the number of passed reads
///              or previous size whichever is smaller -
///             upon return the size of the collector is equal to the number of passed reads or previous size, whichever
///             is larger
/// @param direction for each passed read consider only alignment of the passed sequence as is(FW),
///                  or of both the sequence and its reverse complement (FW_OR_RC), see also useReverse.
/// @param timers mast be a vector of length 2 containing pointers to 2 timer objects; if passed,
///               instruments the method with separate time counts for I/O and for alignment itself
/// @param[in] quals if specified, should point to the array of per-base quality scores for all reads
/// @param useReverse if true, and direction is FW_OR_RC, then in addition to forward alignment, an
///                    attempt will be made to align simple *reverse* of each sequence instead of
///                    its reverse complement (useful for alignments in color space - SOLiD technology).
/// @param max_freq ignore altogether Kmers with frequencies >= max_freq

/*
void ImperfectLookup( lookup_table &look, const vecbasevector& reads,
		      AlignCollectorBase & aligns,
		      AlignDir direction,
		      vec<TaskTimer *> * timers,
		      const vecqualvector * quals = 0,
		      bool useReverse=false, const unsigned int max_freq=0)
{
  // PRINT3( aligns.GetCollectorName(), int( direction ), useReverse );

  const unsigned int K = look.K();
  // Deal with quality scores: verify that the quals match the bases
  // and calculate the mean quality
  double meanQual=0;
  if (quals) {
    vec<int> qs, bs;
    quals->ElementSizes(qs);
    reads.ElementSizes(bs);

    // length of each query sequence must match the length of the
    // corresponding quality vector:
    ForceAssert(bs == qs);
    longlong qsum=0, qcount=0;
    for (int i=0; i != quals->size(); ++i) {
      qsum = accumulate((*quals)[i].begin(), (*quals)[i].end(), qsum);
    }
    qcount=accumulate(qs.begin(), qs.end(), qcount);

    // average quality, per base across all the query seqs:
    meanQual = double(qsum)/qcount;
    PRINT3(meanQual, qsum, qcount);//exit(0);
  }

  // For each query (or its rc), find the indices of each of its kmers.

  // two arrays here, one for all fw, the other for all rc queries;
  // the query object keeps information about the direction and we could
  // go generic here using Query -> Hit transformation that would
  // take care of everything; however this would be a little slower and
  // this particular case is too simple (and important), so we go after
  // performance.
  VecQueryVec all_queries[2];
  BasesToQueries(reads, all_queries[0], all_queries[1], look, max_freq, useReverse);

    //######################
    //#unsigned int TEST_ID = 2;
    //#####################

  const int npasses = ( direction == FW ? 1 : 2 );
  const unsigned int nqueries = (unsigned int)reads.size( );

  //  aligns.clear();
  if ( nqueries > aligns.size() ) {
    aligns.resize( nqueries );
  }

  // Set up for alignment generation. The fields set up below are the same
  // for all the alignments we are going to generate here
  static look_align la;
  la.nhits = la.indels = 0;
  la.a.SetNblocks(1);
  la.a.SetGap( 0, 0 );
  la.a.SetStartOnQuery(0);

  // Go through the lookup table chunks.

  for ( unsigned int i = 0; i < look.NChunks( ); i++ ) {
    if (timers) (*timers)[0]->Start();
    look.ReadChunk(i);
    if (timers) (*timers)[0]->Stop();
    if (timers) (*timers)[1]->Start();

    unsigned int ChunkSize = look.NBasesInChunk(i);
    unsigned int FirstBaseInChunk = look.BasesStart();

    // Go through the query sequences trying to align them in the current chunk.

    for ( la.query_id = 0; (unsigned int)la.query_id < nqueries ; ++la.query_id ) {
      if ( ! aligns.AlignsWanted(la.query_id) ) continue;
      la.query_length = reads[la.query_id].size(); // preset length of this query seq. into look_align
      if ( la.query_length < K ) continue;

      // Go through the orientations:
      for ( int pass = 0; pass < npasses && aligns.AlignsWanted(la.query_id) ; pass++ ) {

	static basevector src;
	if ( pass == 1 ) {
	  if ( useReverse ) src.Reverse(reads[la.query_id]);
	  else src.ReverseComplement(reads[la.query_id]);
	}
	// this is what we are going to align for this query and this pass:
	const basevector& S = ( pass == 0 ? reads[la.query_id] : src );

	static qualvector rq;
	const qualvector * q = 0;
	if (quals) {
	  if (0 == pass) q = &((*quals)[la.query_id]);
	  else {
	    rq = (*quals)[la.query_id];
	    rq.ReverseMe();
	    q = &rq;
	  }
	}

	// this vector is going to hold the offsets of all potential
	// alignments of the current query sequence/orientation based on
	// each individual Kmer. For instance if the 5th Kmer
	// (zero-based numbering) in this query is found at position X in the
	// reference, this results in a potential alignment (offset) of the
	// whole query sequence at position X-5 in the reference.
	static vec<unsigned int> offsets;
 	offsets.clear( );
	QueriesToSeqOffsets(look, (all_queries[pass])[la.query_id].begin(),
			          (all_queries[pass])[la.query_id].end(), offsets);
	UniqueSort(offsets);

	for ( vec<unsigned int>::iterator u = offsets.begin();
	      u != offsets.end( ) && aligns.AlignsWanted( la.query_id ) ; ++u ) {

	  unsigned int offset = *u;
	  // catch wrap-arounds:
	  if ( offset >= (unsigned int) ( K - la.query_length ) ) continue;

	  unsigned int pos_on_tig;
	  // convert the absolute alignment position into contig:position_on_contig:
	  look.GetContigPos( offset, la.target_id, pos_on_tig );

	  // The valid alignment region is bounded by this contig or
	  // the chunk, whichever ends first.
	  unsigned int first = max( look.ContigStart(la.target_id), FirstBaseInChunk );
	  unsigned int last = min( look.ContigStop(la.target_id), look.BasesStop());

	  // We want subsumed alignments only.
	  if ( offset < first || offset + la.query_length > last )
	    continue;

	  // Validate alignment.

	  // stop counting mismatches as soon as we find at least the
	  // same number of mismatches as in the second-best alignment
	  // already observed before:
	  int max_mismatches = aligns.ErrorThreshold(la.query_id);

	  // count mismatches:
	  int mismatches;
	  if ( q == 0 ) {
	    mismatches = MismatchCount(ExactBaseMatch(),
					 S.Begin(),
					 S.End(),
					 look.Bases().Begin( offset - FirstBaseInChunk ),
					 max_mismatches);
	  } else {
	    mismatches = MismatchScore(ExactBaseMatch(),
				       S.Begin(),
				       S.End(),
      				       look.Bases().Begin( offset - FirstBaseInChunk ),
				       q->begin(),
				       max_mismatches*meanQual,
				       meanQual);
	  }

	  //	  cerr << i << ":" << la.query_id << ":" << mismatches << endl;

	  if ( mismatches > max_mismatches ) continue;

	  // Create alignment.
	  la.a.Setpos2(pos_on_tig);
	  la.target_length = look.ContigSize(la.target_id);
	  la.rc1 = ( pass == 1 );
	  la.a.SetLength( 0, la.query_length );
	  la.mutations = (0 == quals)
	    ? mismatches
	    : MismatchCount(ExactBaseMatch(),
			    S.Begin(),
			    S.End(),
			    look.Bases().Begin( offset- FirstBaseInChunk ) );

	  // Try to insert. Collector will figure out the details -
	  // whether the current alignment is better than the best,
	  // or only better then the second best, or not good at all etc

    //###################################
    //#		    if ( (unsigned int)la.query_id == TEST_ID ) {
    //#		      la.PrintReadableBrief(cout);
    //#		      cout << endl;
    //#		    }
    //###################################
	  aligns.Insert(la);
	}
	// End of for ( i=unique_offset_cnt-1 ) loop over all
	// putative alignments (offests) for given id/orientation
	// on the current chunk
      } // end of for (pass ) - loop over forward/reverse alignment passes

    }   // End of for ( id=...) loop over all query sequences;
    if (timers) (*timers)[1]->Stop();
    aligns.Consolidate();
  } // End of for ( int i... ) loop over all chunks

}
*/


class ImperfectLookupAligner {
 public:


  enum {
    IO_TIMER = 0x1,
    CPU_TIMER = 0x02
  };

  ImperfectLookupAligner() :
    m_direction(FW_OR_RC),
    m_instr_level(0),
    m_plookup(0),
    m_owns_lookup(false),
    m_use_reverse(false),
    m_reference_seq_fname(""),
    m_target_seq(0) {}

  ~ImperfectLookupAligner() { if ( m_owns_lookup ) delete m_plookup; }

  /// Set direction(s), in which alignments for every query sequence
  /// will be condsidered: currently FW (forward only) or FW_OR_RC (both fw and rc)
  void SetDirection(AlignDir d) { m_direction = d; }

  /// Getter method: shows what direction(s) this aligner is going to
  /// try when aligning every query seq.
  AlignDir GetDirection(AlignDir d) { return m_direction; }

  /// Tells aligner to set instrument flag(s)
  void InstrumentOn(unsigned int flag) { m_instr_level |= flag; }

  /// Tells aligner to clear instrument flag(s)
  void InstrumentOff(unsigned int flag) { m_instr_level &= ( ~flag ); }

  /// Sets lookup table (index) to use when searching for alignments;
  /// this method loads lookup table from the file on disk (specified by file name);
  /// if a lookup table was already loaded and owned by this aligner class, it will be
  /// deleted. The lookup table loaded by this method is internal for the aligner class and
  /// owned by it.
  void SetTargetLookup(const String & lookup_name) {
    if ( m_owns_lookup ) delete m_plookup;
    m_plookup = new lookup_table(lookup_name);
    m_K = m_plookup->K();
    m_owns_lookup = True;
  }

  /// Sets lookup table (index) to use when searching for alignments;
  /// this method reuses lookup table that was already loaded elsewhere.
  /// Aligner does \em not take the ownership of the lookup table passed
  /// via this method. The caller should take care of destroying it when
  /// the alignments are done.
  /// NOTE: unsafe method! lookup tables can not be copied, so it is the
  /// passed instance of the lookup table that will be used; do not destroy
  /// externally loaded lookup table prematurely or the Aligner will crash!!
  void SetTargetLookup(lookup_table & l ) {
    if ( m_owns_lookup ) delete m_plookup;
    m_plookup = & l;
    m_K = m_plookup->K();
    m_owns_lookup = False;
  }

  /// Sets the name of the file (fastb) where the target reference
  /// sequence resides. Target reference sequence is used by few
  /// less-trivial alignment strategies (bisulfite, SOLiD), see the
  /// corresponding documentation pages. The "standard" alignment
  /// strategy reads target sequence directly from the lookup table and
  /// makes no use of the "original" target reference passed via this method.
  void SetTargetSequence(const String & seq_name) {
    m_reference_seq_fname = seq_name;
  }



  /// If 'reverse' is set to <True>, the queries ("reads") will be
  /// simply reverted rather than reverse complemented when looking
  /// for the "other strand" alignments ('reverse' mode has any effect
  /// at all only if direction is FW_OR_RC, i.e. alignments are sought
  /// on both strands). This is a trick and a workaround for SOLiD alignments.
  void SetReverse(Bool rev) { m_use_reverse = rev; }

  /// Returnst the 'reverse' status (<True> if a query is going
  /// to be reversed rather than reverse complemented when looking
  /// for the "other strand" match).
  Bool GetReverse() { return m_use_reverse; }

  template <typename AlignCollector>
  void ComputeAlignments(const vecbasevector & queries, AlignCollector & aligns) {
      RunAlignmentMain(queries, 0, aligns, STANDARD_ALIGN);
  }

  template <typename AlignCollector>
  void ComputeColorspaceAlignments(const vecbasevector & queries, AlignCollector & aligns) {
      if ( ! m_reference_seq_fname.empty() ) {
	//	  m_target_seq.ReadAll(m_reference_seq_fname);
      } else {
	  cout << "ImperfectLookupAligner setup is incomplete: target basespace reference must be specified for " << endl
	       << "full colorspace read (SOLiD) alignment strategy" << endl;
	  exit(1);
      }

      SetReverse(True);
      RunAlignmentMain(queries, 0, aligns, SOLID_FULL_ALIGN);
  }

  static double computeAverageQuality( vecbvec const& reads, vecqvec const& quals );

 protected:

  enum AlignmentStrategy {
    STANDARD_ALIGN = 0,
    SOLID_FULL_ALIGN = 1 } ;

  /// Actual alignment logic is implemented here (this method is called by public interfaces)
  template <typename AlignCollector>
  void RunAlignmentMain(const vecbasevector & query_seqs,
			const vecqualvector * query_quals,
			AlignCollector & aligns,
			AlignmentStrategy strategy=STANDARD_ALIGN) {
    PRINT3( aligns.GetCollectorName(), int( m_direction ), m_use_reverse );

    if (m_plookup == 0) {  // lookup table must exist
      cout << "ImperfectLookupAligner setup is incomplete: lookup index is not specified" << endl;
      exit(1);
    }


    if ( strategy == SOLID_FULL_ALIGN && m_reference_seq_fname.empty() ) {
      cout << "ImperfectLookupAligner setup is incomplete: target basespace reference must be specified for " << endl
	   << "full colorspace read (SOLiD) alignment strategy" << endl;
      exit(1);
    }


    // Deal with quality scores: verify that the quals match the bases
    // and calculate the mean quality
    double meanQual = 0.;
    if (query_quals) {
       if ( strategy != STANDARD_ALIGN ) {
	   cout << "Quality scores are not supported in the requested alignment strategy" << endl;
	   exit(1);
       }
       meanQual = computeAverageQuality(query_seqs,*query_quals);
    }

    // For each query (or its rc), find the indices of each of its kmers.

    // two arrays here, one for all fw, the other for all rc queries;
    // the query object keeps information about the direction and we could
    // go generic here using Query -> Hit transformation that would
    // take care of evrything; however this would be a little slower and
    // this particular case is too simple (and important), so we go after
    // performance.
    VecQueryVec all_queries[2];
    const unsigned int nqueries = (unsigned int)query_seqs.size( );

    basevector first_bases; // used for full colorspace alignment only;
                            // will keep the first base after the primer
                            // as suggested by the primer base + 1st color
    unsigned int start_pos = 0;
    if ( strategy == SOLID_FULL_ALIGN ) {
      start_pos = 2; // colorspace reads with primer!
      first_bases.resize(nqueries);
      for ( unsigned int i = 0 ; i < nqueries ; i++ ) {
	first_bases.Set(i, SOLiD_toBase(query_seqs[i][0],query_seqs[i][1]) );
      }
    }

    // transform query sequences into queries (kmers), separately for fw and rc:
    BasesToQueries(query_seqs, all_queries[0], all_queries[1], *m_plookup, 0, m_use_reverse,start_pos);

    const int npasses = ( m_direction == FW ? 1 : 2 );

    if ( nqueries > aligns.size() ) {
      aligns.resize( nqueries );
    }

    auto_ptr<TaskTimer> t1, t2;

    if ( m_instr_level & IO_TIMER != 0 ) t1.reset( new TaskTimer() );
    if ( m_instr_level & CPU_TIMER != 0 ) t2.reset( new TaskTimer() );

    // Set up for alignment generation. The fields set up below are the same
    // for all the alignments we are going to generate here
    static look_align la;
    la.nhits = la.indels = 0;
    la.a.SetNblocks(1);
    la.a.SetGap( 0, 0 );
    la.a.SetStartOnQuery(0);

    //######################
    // ## unsigned int TEST_ID = 2;
    //#####################

    // Go through the lookup table chunks.

    if ( m_plookup->NChunks() > 1 ) cout << "lookup chunks:" ;
    for ( unsigned int i = 0; i < m_plookup->NChunks( ); i++ ) {
        if ( m_instr_level & IO_TIMER != 0 ) t1->Start();
	m_plookup->ReadChunk(i);
	if ( m_plookup->NChunks() > 1 ) {
	  cout << "." ;
	  flush(cout);
	}

	unsigned int first_contig_in_chunk = m_plookup->FirstContigInChunk();
	unsigned int contigs_in_chunk = m_plookup->ContigsInChunk();
	if ( strategy == SOLID_FULL_ALIGN ) {
	  m_target_seq.clear();
	  m_target_seq.ReadRange(m_reference_seq_fname,
				 first_contig_in_chunk,
				 first_contig_in_chunk+contigs_in_chunk);
	}

        if ( m_instr_level & IO_TIMER != 0 ) t1->Stop();

	if ( m_instr_level & CPU_TIMER != 0 ) t2->Start();

	unsigned int ChunkSize = m_plookup->NBasesInChunk(i);
	unsigned int FirstBaseInChunk = m_plookup->BasesStart();

	// Go through the query sequences trying to align them in the current chunk.

	for ( la.query_id = 0; (unsigned int)la.query_id < nqueries ; ++la.query_id ) {
	    if ( ! aligns.AlignsWanted(la.query_id) ) continue;
	    la.query_length = query_seqs[la.query_id].size() - start_pos; // preset length of this query seq. into look_align
	    if ( (unsigned int) la.query_length < m_K ) continue;

	    // Go through the orientations:
	    for ( int pass = 0; pass < npasses && aligns.AlignsWanted(la.query_id) ; pass++ ) {

	        static basevector src;
	        if ( pass == 1 ) {
		    if ( m_use_reverse ) src.Reverse(query_seqs[la.query_id]);
		    else src.ReverseComplement(query_seqs[la.query_id]);
		}

		// this is what we are going to align for this query and this pass:
		const basevector& S = ( pass == 0 ? query_seqs[la.query_id] : src );

		static qualvector rq;
		const qualvector * q = 0;
		if (query_quals) {
		    if (0 == pass) q = &((*query_quals)[la.query_id]);
		    else {
		        rq = (*query_quals)[la.query_id];
			rq.ReverseMe();
			q = &rq;
		    }
		}

		// this vector is going to hold the offsets of all potential
		// alignments of the current query sequence/orientation based on
		// each individual Kmer. For instance if the 5th Kmer
		// (zero-based numbering) in this query is found at position X in the
		// reference, this results in a potential alignment (offset) of the
		// whole query sequence at position X-5 in the reference.
		static vec<unsigned int> offsets;

		offsets.clear( );

		// lookup queries (kmers) in the index table and record all offsets:
		QueriesToSeqOffsets(*m_plookup,
				    (all_queries[pass])[la.query_id].begin(),
			            (all_queries[pass])[la.query_id].end(),
				    offsets);
		UniqueSort(offsets);

		for ( int u = 0; u < offsets.isize( ) && aligns.AlignsWanted( la.query_id ) ; u++ ) {

		    unsigned int offset = offsets[u];
		    // catch wrap-arounds:
		    if ( offset >= (unsigned int) ( m_K - la.query_length ) ) continue;

		    unsigned int pos_on_tig;
		    // convert the absolute alignment position into contig:position_on_contig:
		    m_plookup->GetContigPos( offset, la.target_id, pos_on_tig );

		    // The valid alignment region is bounded by this contig or
		    // the chunk, whichever ends first.
		    unsigned int first = max( m_plookup->ContigStart(la.target_id), FirstBaseInChunk );
		    unsigned int last = min( m_plookup->ContigStop(la.target_id), m_plookup->BasesStop());

		    // We want subsumed alignments only.
		    if ( offset < first || offset + la.query_length > last )
		      continue;

		    // Validate alignment.

		    // will iterate along the target sequence
		    basevector::const_iterator target_iter=m_plookup->Bases().begin( offset - FirstBaseInChunk );
		    if ( m_reference_seq_fname.empty() || strategy == SOLID_FULL_ALIGN ) {
		      // if we use target sequence directly as stored in the lookup (we always use that for
		      // colorspace reads!):
		      //ALREADY DONE:		      target_iter = m_plookup->Bases().Begin( offset - FirstBaseInChunk );
		    } else {
		      // if the target reference sequence was explicitly supplied:
		      target_iter = m_target_seq[la.target_id-first_contig_in_chunk].Begin(pos_on_tig);
		    }

		    // stop counting mismatches as soon as we find at least the
		    // same number of mismatches as in the second-best alignment
		    // already observed before:
		    int max_mismatches = aligns.ErrorThreshold(la.query_id);

		    // count mismatches:
		    int mismatches;
		    if ( q == 0 ) {
		      if ( pass == 0 ) {
			mismatches = MismatchCount(ExactBaseMatch(),
						   S.Begin(start_pos),
						   S.End(),
						   target_iter,
						   max_mismatches);
		      } else {
			basevector::const_iterator end = S.End();
			end -= start_pos;
			mismatches = MismatchCount(ExactBaseMatch(),
						 S.Begin(),
						 end,
						 target_iter,
						 max_mismatches);
		      }

		    } else {
		      mismatches = MismatchScore(ExactBaseMatch(),
						 S.Begin(),
						 S.End(),
						 target_iter,
						 q->begin(),
						 max_mismatches*meanQual,
						 meanQual);
		    }
    /* ######################################
	#	    if ( (unsigned int)la.query_id == TEST_ID) PRINT2(mismatches,max_mismatches);
	#
 	########## */
		    if ( strategy == SOLID_FULL_ALIGN ) {
		      basevector & ref_contig = m_target_seq[la.target_id-first_contig_in_chunk];
		      base_t ref_base = ( ( pass==0 ) ?
					  ref_contig[pos_on_tig] :
					  3-ref_contig[pos_on_tig+la.query_length] );
		      if ( first_bases[la.query_id] != ref_base ) {

			//#		cout << "FIRST BASE WRONG!" << endl;
			mismatches++; // first base mismatch!
		      }
		    }

		    //	  cerr << i << ":" << la.query_id << ":" << mismatches << endl;

		    if ( mismatches > max_mismatches ) continue;

		    // Create alignment.
		    la.a.Setpos2(pos_on_tig);
		    la.target_length = m_plookup->ContigSize(la.target_id);
		    la.rc1 = ( pass == 1 );
		    la.a.SetLength( 0, la.query_length );
		    la.mutations = (0 == query_quals)
		      ? mismatches
		      : MismatchCount(ExactBaseMatch(),
				      S.Begin(),
				      S.End(),
				      target_iter);

		    // Try to insert. Collector will figure out the details -
		    // whether the current alignment is better than the best,
		    // or only better then the second best, or not good at all etc
	  /*###################################
	    #          if ( (unsigned int)la.query_id == TEST_ID ) {
	    #	          basevector btmp;
	    #             la.PrintReadableBrief(cout);
	    #             la.target_id-=first_contig_in_chunk;
	    #             btmp.SetToSubOf(query_seqs[la.query_id],2,query_seqs[la.query_id].size()-2);
	    #             PrintVisualColorAlign(cout, btmp,
	    #                    m_target_seq, la,
	    #			 SOLiD_toBase(query_seqs[la.query_id][0],
	    #							query_seqs[la.query_id][1]));
	    #	          cout << endl;
	    #             la.target_id+=first_contig_in_chunk;
	    #	       }
            ###################################	*/

		    aligns.Insert(la);
		}
		// End of for ( i=unique_offset_cnt-1 ) loop over all
		// putative alignments (offests) for given id/orientation
		// on the current chunk
	    } // end of for (pass ) - loop over forward/reverse alignment passes

	}   // End of for ( id=...) loop over all query sequences;
	if ( m_instr_level & CPU_TIMER != 0 ) t2->Stop();
	aligns.Consolidate();
    } // End of for ( int i... ) loop over all chunks
    if ( m_plookup->NChunks() > 1 ) {
      cout << endl ;
      flush(cout);
    }

    if ( m_instr_level & IO_TIMER ) cout << "Loading time: " << *t1.get() << endl;
    if ( m_instr_level & CPU_TIMER ) cout << "Aligning time: " << *t2.get() << endl;
  }

  /*
  template <typename AlignCollector>
  void OffsetsToAligns(unsigned int query_id,
		       const basevector & S,
		       vec<unsigned int> & offsets,
		       AlignCollector & aligns,
		       vecqualvector *q = 0,
		       meanQual = 0 ) {

       look_align la;
       la.nhits = la.indels = 0;
       la.a.SetNblocks(1);
       la.a.SetGap( 0, 0 );
       la.a.SetStartOnQuery(0);
       la.query_id = query_id;
       la.query_length = S.size();
       la.a.SetLength( 0, la.query_length );

       unsigned int first_base_in_chunk = m_plookup->BasesStart();
       unsigned int first_contig_in_chunk = m_plookup->FirstContigInChunk();

       for ( int u = 0; u < offsets.isize( ) && aligns.AlignsWanted( la.query_id ) ; u++ ) {

	   unsigned int offset = offsets[u];
	   // catch wrap-arounds:
	   if ( offset >= (unsigned int) ( m_K - la.query_length ) ) continue;

	   unsigned int pos_on_tig;
	   // convert the absolute alignment position into contig:position_on_contig:
	   m_plookup->GetContigPos( offset, la.target_id, pos_on_tig );

	   // The valid alignment region is bounded by this contig or
	   // the chunk, whichever ends first.
	   unsigned int first = max( m_plookup->ContigStart(la.target_id), first_base_in_chunk );
	   unsigned int last = min( m_plookup->ContigStop(la.target_id), m_plookup->BasesStop());

	   // We want subsumed alignments only.
	   if ( offset < first || offset + la.query_length > last )
	       continue;

	   // Validate alignment.

		    // will iterate along the target sequence
	   basevector::iterator target_iter;
	   if ( m_reference_seq.Empty() ) {
	       // if we use target sequence directly as stored in the lookup:
	       target_iter = m_plookup->Bases().Begin( offset - first_base_in_chunk );
	   } else {
	       // if the target reference sequence was explicitly supplied:
	       target_iter = m_reference_seq[la.target_id-first_contig_in_chunk].Begin(pos_on_tig);
	   }

	   // stop counting mismatches as soon as we find at least the
	   // same number of mismatches as in the second-best alignment
	   // already observed before:
	   int max_mismatches = aligns.ErrorThreshold(la.query_id);

	   // count mismatches:
	   int mismatches;
	   if ( q == 0 ) {
  	       mismatches = MismatchCount(ExactBaseMatch(),
					  S.Begin(),
					  S.End(),
					  target_iter,
					  max_mismatches);
	   } else {
	       mismatches = MismatchScore(ExactBaseMatch(),
					  S.Begin(),
					  S.End(),
					  target_iter,
					  q->begin(),
					  max_mismatches*meanQual,
					  meanQual);
	   }

	   //	  cerr << i << ":" << la.query_id << ":" << mismatches << endl;

	   if ( mismatches > max_mismatches ) continue;

	   // Create alignment.
	   la.a.Setpos2(pos_on_tig);
	   la.target_length = m_plookup->ContigSize(la.target_id);
	   la.rc1 = ( pass == 1 );
	   la.mutations = (0 == quals)
	     ? mismatches
	     : MismatchCount(ExactBaseMatch(),
			     S.Begin(),
			     S.End(),
			     target_iter);

	   // Try to insert. Collector will figure out the details -
	   // whether the current alignment is better than the best,
	   // or only better then the second best, or not good at all etc

	   aligns.Insert(la);
       }
  }
  */
 private:

  AlignDir m_direction; ///< perform alignments in only one or in both directions
  int m_instr_level;  ///< stores instrumentation attributes
  lookup_table *m_plookup; ///< pointer to the lookup table (index)
  Bool m_owns_lookup;      ///< is lookup table pointed to by this aligner owned by it?
  Bool m_use_reverse;      ///< use reverse instead of reverse complement when aligning to the second strand
  String m_reference_seq_fname; ///< if empty (default), the reference seq stored in lookup table will be used

  vecbasevector m_target_seq; ///< stores target sequence (if different from the one stored in lookup index)

  unsigned int m_K; // K-mer size (set automatically when lookup table is set)
};


/// This is a convenience method that breaks large alignment job into parts, sends those parts
/// to the specified LSF queue, waits for all chunks to complete and merges the partial results into
/// the final output file(s). This method does not use LSF package of the codebase and should probably
/// be replaced with methods from that package (eventually).
///
/// force_wrire_level: 0 - do not recompute anything that already exists
///                    1 - recompute full alignment file, but do not recompute
///                        batch-aligned chunks that finished successfully
///                    2 - recompute all (all batch chunks will be resubmitted too)

void AlignInChunksOnFarm_ILT(const String & readfile,
			 const String & reffile,
			 const String & out_head,
			 const String & out_suffix,
			 const String & LSF_QUEUE,
			 const String & align_params_imp,
			 int batch_size=300000,
			     int force_write_level=0 );


#endif
