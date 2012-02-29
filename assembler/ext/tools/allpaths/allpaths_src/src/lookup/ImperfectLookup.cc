/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "lookup/AlignCollector.h"
#include "lookup/ImperfectLookup.h"
#include "TaskTimer.h"


/// Overload that reads in the lookup file from disk and passes execution
/// to ImperfectLookup(lookup_table &, vecbasevector &, vec<look_align> &,...).
/// This is a legacy method with flat list of parameters which is too long and
/// with more rigid logic,
/// see more flexible overload templatized by AlignCollector.
void ImperfectLookup( unsigned int K, const vecbasevector& reads,
		      const String& lookup_file, vec<look_align>& aligns,
		      vec<int>& best_errors,
		      AlignDir direction, int best_prox, int max_errors,
		      vec<Bool>* unique_aligned,
		      vec<TaskTimer *> * timers,
		      const vecqualvector * quals,
		      bool useReverse, const unsigned int max_freq)
{
  // Read header information from lookup table.

  lookup_table look(lookup_file);
  ForceAssertEq( look.K( ), K );
  ImperfectLookup(look, reads, aligns, best_errors, direction, best_prox, max_errors,
		  unique_aligned, timers, quals,useReverse, max_freq);
}


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
     vec<TaskTimer *> * timers ,
     const vecqualvector * quals,
     bool useReverse ,unsigned int max_freq) {

  lookup_table look(lookup_file);
  ImperfectLookup(look, reads, aligns, direction, timers, quals, useReverse, max_freq);
}

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
void ImperfectLookup( lookup_table &look, const vecbasevector& reads,
		      vec<look_align>& aligns, vec<int>& best_errors,
		      AlignDir direction, int best_prox, int max_errors,
		      vec<Bool>* unique_aligned,
		      vec<TaskTimer *> * timers,
		      const vecqualvector * quals,
		      bool useReverse, const unsigned int max_freq )
{
  // Initialize unique_aligned.

  if ( unique_aligned != 0 )
    unique_aligned->resize_and_set( reads.size( ), False );


  if ( max_errors < 0 ) max_errors = infinitely_many;
  UniqueByErrDiffAlignCollector aligns_tmp(best_prox,max_errors);

  ImperfectLookup(look, reads, aligns_tmp, direction, timers, quals, useReverse,
     max_freq);

  // Finish up.

  unsigned int nqueries = reads.size();
  best_errors.resize( nqueries );

  if ( unique_aligned == 0 ) aligns.clear( );
  for ( unsigned int id = 0; id < nqueries ; id++ ) {
    best_errors[id] = aligns_tmp.MinErrors(id);
    if ( best_errors[id] < 0 ) best_errors[id] = infinitely_many;
    if ( aligns_tmp.UniquelyAligned(id) ) {
      if ( unique_aligned == 0 ) aligns.push_back( aligns_tmp.Align(id) );
      else (*unique_aligned)[id] = True;
    }
  }
}

void ImperfectLookup( lookup_table &look, const vecbasevector& reads,
		      AlignCollectorBase & aligns,
		      AlignDir direction,
		      vec<TaskTimer *> * timers,
		      const vecqualvector * quals,
		      bool useReverse, const unsigned int max_freq)
{
  // PRINT3( aligns.GetCollectorName(), int( direction ), useReverse );

  const unsigned int K = look.K();
  // Deal with quality scores: verify that the quals match the bases
  // and calculate the mean quality
  double meanQual=0.;
  if ( quals )
      meanQual = ImperfectLookupAligner::computeAverageQuality(reads,*quals);

  // For each query (or its rc), find the indices of each of its kmers.

  // two arrays here, one for all fw, the other for all rc queries;
  // the query object keeps information about the direction and we could
  // go generic here using Query -> Hit transformation that would
  // take care of evrything; however this would be a little slower and
  // this particular case is too simple (and important), so we go after
  // performance.
  VecQueryVec   all_queries[2];
  BasesToQueries(reads, all_queries[0], all_queries[1], look, max_freq, useReverse);

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
	for ( int u = 0; u < offsets.isize( ) && aligns.AlignsWanted( la.query_id ) ; u++ ) {

	  unsigned int offset = offsets[u];
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



/// This is a convenience method that breaks large alignment job into parts, sends those parts
/// to the specified LSF queue, waits for all chunks to complete and merges the partial results into
/// the final output file(s). This method does not use LSF package of the codebase and should probably
/// be replaced with methods from that package (evsentually).
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
			 int batch_size,
			 int force_write_level ) {
          // Below is so we do not have to keep the .bsub.output files

    String outfile = out_head+out_suffix;

    if ( force_write_level == 0  && IsRegularFile( outfile ) ) return; // output file already exists

    cout << "Alignment file missing. Starting new alignment." << endl;

    // total number of reads to align:
    unsigned int nreads = MastervecFileObjectCount(readfile);
    // total number of batch jobs to spawn:
    int nbatches = ( nreads + batch_size - 1 ) / batch_size;
    int batch = 1; // current batch (counter)

    for ( unsigned int start = 0; start < nreads ; start += batch_size, ++batch )  {
        String BATCH = ToString(batch);
	String bsub_output = out_head + ".bsub.output." + BATCH;

	// don't recompute batch jobs that finished successfully (unless force_write is specified)
	if ( ( force_write_level < 2 ) &&  IsRegularFile(bsub_output) &&
	     ( LineOfOutput( "tail -1 " + bsub_output, True ).Contains( ": done.", -1 ) ||
	       LineOfOutput( "tail -2 " + bsub_output + " | head -1", True ).Contains( ": done.", -1 ) )
	   )  continue;

	int stop = Min( start + batch_size, nreads );
	Remove( out_head + ".bsub.out." + BATCH );
	Remove( out_head + ".bsub.err." + BATCH );
	Remove( out_head + ".bsub.output." + BATCH );
	Remove( outfile+"."+BATCH); // remove partial alignment file
	String command = "bsub -q " + LSF_QUEUE + " -o " + out_head + ".bsub.out."
	                 + BATCH + " " + "-e " + out_head + ".bsub.err." + BATCH + " "
	                 + "-R \"rusage[mem=2048]\" "
	                 + "\"ImperfectLookupTable K=12 SEQS=" + readfile + " "
	                 + "L=" + reffile + " "
	                 + "MIN_ERRORS_PREFIX=." + BATCH + " " + "O=" + out_head + " "
	                 + "OUT_SUFFIX=" + out_suffix +"." + BATCH + " "
	                 + align_params_imp + " "
	                 + "START=" + ToString(start) + " END=" + ToString(stop)
	                 + " > " + bsub_output + "\"" ;
	cout << command << endl;
	flush(cout);
	SystemSucceed(command);
    } // end for ( int start = 0 ; start < nreads; ...)
    // all batch jobs (if force_write_level==2) all all missing batch jobs
    // (if force_write_level < 2 ) are now submitted

    // Now wait for the alignments to finish:

    while(1) {
        sleep(10);
	int j;
	for ( j = 1; j <= nbatches; j++ ) {
	    String bsub_output = out_head + ".bsub.output." + ToString(j);
	    if ( !IsRegularFile(bsub_output) ) break; // no output file yet, keep waiting
	    if ( !LineOfOutput( "tail -1 " + bsub_output, True ).Contains( ": done.", -1 ) &&
		 !LineOfOutput( "tail -2 " + bsub_output + " | head -1", True ).Contains( ": done.", -1 )
		 )    break;  // output file does not contain 'done' line yet; keep waiting
	}
	if ( j > nbatches ) break; // we got all 'nbatches' output files accounted for; all contain 'done'
    }
    // all batch jobs finished!!

    // Merge global alignment files.

    Remove( outfile );
    Remove( out_head + ".minAlignErrors.txt" );

    // merge qltout and minAlignErrors files, and also batch output and error files;
    // after everything is merged, chunk alignments and output/error files are deleted!
    for ( int j = 1; j <= nbatches; j++ ) {
        String aligns = out_head + out_suffix + "." + ToString(j);
	Cp2( aligns, outfile , True );
	String errs = out_head + "." + ToString(j) + ".minAlignErrors.txt";
	if ( IsRegularFile(errs) ) {
	       Cp2( errs, out_head + ".minAlignErrors.txt", True );
	       Remove(errs);
	}
	String output =  out_head + ".bsub.output." + ToString(j);
	String outerrs =  out_head + ".bsub.err." + ToString(j);
	Cp2( output, out_head + out_suffix.SafeBefore(".qltout") + ".bsub.output",True);
	Cp2( outerrs, out_head + out_suffix.SafeBefore(".qltout") + ".bsub.err",True);
	Remove(aligns);
	Remove(output);
	Remove(outerrs);
	Remove(out_head + ".bsub.out." + ToString(j));
    }
}


double ImperfectLookupAligner::computeAverageQuality( vecbvec const& reads,
                                                      vecqvec const& quals )
{
    ForceAssert(reads.size()==quals.size());

    unsigned long qsum = 0;
    unsigned long qcount = 0;
    vecbvec::const_iterator rEnd(reads.end());
    vecbvec::const_iterator rItr(reads.begin());
    vecqvec::const_iterator qItr(quals.begin());
    for ( ; rItr != rEnd; ++rItr, ++qItr )
    {
      qvec const& qs = *qItr;
      ForceAssert(rItr->size()==qs.size());
      qcount += qs.size();
      qsum = accumulate(qs.begin(), qs.end(), qsum);
    }

    // average quality, per base across all the query seqs:
    double meanQual = static_cast<double>(qsum) / qcount;
    PRINT3(meanQual, qsum, qcount);
    return meanQual;
}
