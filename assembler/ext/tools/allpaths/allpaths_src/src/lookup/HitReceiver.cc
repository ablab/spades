/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "lookup/HitReceiver.h"

void UniqueGlobalUngappedHitReceiver::operator()(unsigned int q,
						    unsigned int queryPos, bool isRc,
						    unsigned int offset, unsigned int contig,
						    unsigned int )
{
  // PRINT5(q, queryPos, isRc, offset, contig);
  if ( ProvedAmbiguous(q) ) return;  // Only want unambiguous query alignments
  if ( offset < queryPos + look.BasesStart() ) return; // Only want aligns within current chunk
  unsigned int s0 = 0, r0 = offset - look.BasesStart() - queryPos; // start positions in S, R
  // PRINT(r0);

  unsigned int targetPos = offset - look.ContigStart(contig);
  if (targetPos < queryPos) return; // Only want global aligns of query within contig
  int alignStartPos = int(targetPos) - int(queryPos); 
  // PRINT2(targetPos, alignStartPos);

  const basevector &R = look.Bases(); // The bases available in current chunk
  SelectQuery(q, isRc); // get member basevector S set appropriately
  unsigned int length = S.size();
  if (r0 + length > R.size()) return; // Only want aligns within current chunk
  if (offset - queryPos + length > look.ContigStop(contig)) return; // Only want global query aligns
  // Compute alignment quality.
  int mismatches = 0;
  int max_mismatches = Min( best_errors[q] + best_prox, second_best_errors[q] - 1 );
  for ( unsigned int y = 0; y < length; y++ ) {
    if ( S[s0 + y] != R[ r0 + y ] ) {
      ++mismatches;
      if ( mismatches > max_mismatches ) 
	break; // No need to keep counting...
    }
  }
  // PRINT3(mismatches, best_errors[q], second_best_errors[q]);

  // If the alignment is equally good, we might have found the
  // same alignment again (should only happen in the overlap between
  // chunks.)
  if ( mismatches == best_errors[q]
       && aligns[q].target_id == int(contig)
       && aligns[q].a.pos2() == alignStartPos )
    return;

  // Now that we've ruled that problem out, if this alignment is not
  // an improvement, it need not be generated.
  if ( mismatches >= best_errors[q] ) {
    second_best_errors[q] = Min( second_best_errors[q], mismatches );
    return; 
  }

  // This alignment is an improvement, but if it's not good enough,
  // it need not be generated.
  if ( mismatches > max_errs ) {
    second_best_errors[q] = best_errors[q];
    best_errors[q] = mismatches;
    return;
  }

  // Create alignment.
  la.query_id = q;
  la.target_id = contig;
  la.a.Setpos1(s0);
  la.a.Setpos2(alignStartPos);
  la.query_length = length;
  la.target_length = look.ContigSize(contig);
  la.rc1 = isRc;
  la.a.SetLength( 0, length );
  la.mutations = mismatches;

  // Update. We already know this is an improvement.
  second_best_errors[q] = best_errors[q];
  best_errors[q] = mismatches;
  aligns[q] = la;
  // la.PrintReadableBrief(cout);
  return;
}



