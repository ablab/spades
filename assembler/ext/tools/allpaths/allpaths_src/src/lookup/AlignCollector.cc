/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "lookup/AlignCollector.h"
#include "VecTemplate.h"
BINARY3_DEF(look_align);

// Use this align to test read for ambiguity
void AmbiguousAlignCollector::Insert(const look_align & la) {
  int i = la.query_id;

  if ( counts_per_read_[i] == num_better_aligns_) {
    return;
  }

  if ( (unsigned int)la.Errors() < best_[i] ) {
    best_[i] = la.Errors();
    counts_per_read_[i] = num_better_aligns_;
  }
  else if ( (unsigned int)la.Errors() == best_[i] ) {
    ++counts_per_read_[i];
  }

//   cout << "ambiguous_read: " << la.Errors() <<" "<< counts_per_read_[i] <<  flush;
//   if ( counts_per_read_[i] == num_better_aligns_ )
//     cout << " AMB " <<flush;
//   else
//     cout << " OK " <<flush;
//   cout << " thresh_per_read: " << threshold_per_read_[i] << endl;
//   cout << endl;
}



/// Add in a look_align (or don't if we don't want it).
void MaxErrDiffAlignCollector::Insert(const look_align & la) {
  int i = la.query_id;
  // We need to keep alignments that are worse than maxErrs_ but could be
  // multiple alignments for best_[i] when (if) we find it...  note: best_[i] is
  // always initialized with 'infinitely_many', so that the first "real"
  // alignment always passes the first condition below; after that we are always
  // and only interested in alignments that are not worse than current
  // best[i]+maxErrDiff (best[i] will be also adjusted below if required); if at
  // the end all the alignments (including best_[i]) are worse than maxErr they
  // will be erased later.
  if ( (unsigned int)la.Errors() <= best_[i] + maxErrDiff_ /* || la.Errors() <= maxErrs_ */) {
    consolidated_ = false;
    aligns_[i].push_back(la);
    if ((unsigned int)la.Errors() < best_[i]) {
      best_[i] = la.Errors();
      EraseBad(i);
      ++betterBest_;
    }
  }
}

/// Inspect all the alignments currently stored for query \c i and remove those
/// that have at least maxErrDiff+1 more errors than the best alignment found.
void MaxErrDiffAlignCollector::EraseBad(int i) {
  int bad = best_[i] + maxErrDiff_;
  vec<look_align>::iterator firstBad =
    remove_if(aligns_[i].begin(),aligns_[i].end(), gt_lookalign_Errors(bad));
  aligns_[i].erase(firstBad, aligns_[i].end());
}


void MaxErrDiffAlignCollector::Erase(int i) {
  aligns_[i].clear();
  best_[i] = infinitely_many;
}



/// Clean up: remove any look_aligns that are not good enough.
void MaxErrDiffAlignCollector::Consolidate() {

  if (true == consolidated_) return; // nothing to do

  longlong capacity=0, size=0;
  for (int i=0; i != aligns_.isize(); ++i) {
      capacity += aligns_[i].capacity(); size += aligns_[i].size();
  }
  //  PRINT3("before", capacity, size);

  for (int i=0; i != aligns_.isize(); ++i) {
      // erase all alignments for query i that have more than best[i]+maxErr errors
      EraseBad(i);

      // Clear aligns where the best alignment is worse than maxErrs_
      // Important - this will not remove alignments which have more than maxErrs_
      // but less than best[i] + maxErrDiff + 1.
      if (!aligns_[i].empty() && best_[i] > maxErrs_) {
	vec<look_align>().swap(aligns_[i]); //clear size and capacity.
      }

      // Make sure we have no duplicate alignments
      sort(aligns_[i].begin(),aligns_[i].end(), order_lookalign_TargetQueryStartEnd());
      aligns_[i].erase(unique(aligns_[i].begin(),
			      aligns_[i].end(),
			      equal_lookalign_TargetQueryStartEnd() ),
		       aligns_[i].end());
      // re-sort by error rate:
      sort(aligns_[i].begin(), aligns_[i].end(), order_lookalign_ErrorRate());
  }

  capacity = size=0;
  for (int i=0; i != aligns_.isize(); ++i) {
      capacity += aligns_[i].capacity(); size += aligns_[i].size();
  }
  //  PRINT2(capacity, size);
  //  PRINT(betterBest_);
  betterBest_=0;
  consolidated_ = true;
}

/// Add in a look_align (or don't if we don't want it).

void UniqueByErrDiffAlignCollector::Insert(const look_align & la) {
  int q = la.query_id;

  if (!AlignsWanted(q)) {
    //      if ( la.Errors() < best_[q] ) {
    //	  second_best_[q] = best_[q];
    //          best_[q] = la.Errors();
    //      }
    return;
  }

  // If the alignment is equally good, we might have found the
  // same alignment again
  if ( (unsigned int)la.Errors() == best_[q] && aligns_[q] == la )
    return;

  // Now that we've ruled that problem out, if this alignment is not
  // an improvement, it need not be saved. it still can be better than the
  // current second best, so we take care of that:
  if ( (unsigned int)la.Errors() >= best_[q] ) {
    second_best_[q] = Min( second_best_[q], (unsigned int)la.Errors() );
    return;
  }

  // This alignment is an improvement, but if it's not good enough,
  // it need not be saved.
  if ( (unsigned int)la.Errors() > maxErrs_ ) {
    second_best_[q] = best_[q];
    best_[q] = la.Errors();
    return;
  }

  // Save the alignment.
  second_best_[q] = best_[q];
  best_[q] = la.Errors();
  aligns_[q] = la;

  ++betterBest_;
}



/// Add in a look_align (or don't if we don't want it).

void BestNextBestAlignCollector::Insert(const look_align & la) {
  int q = la.query_id;
  unsigned int e = la.Errors();


  // if this alignment is not
  // an improvement of at least second_best, it need not be saved.
  if ( e > second_best_[q] ) return;

  // same alignment, no need to save; this line guarantees that
  // we never allow a duplication of the *best* alignment we keep;
  // note that the *second* best alignments might contain duplicates
  if ( e == best_[q] && la == aligns_[q][0] ) return;

  // alignment as good as current second best; if it's worth saving, add it and bail out
  // (note that here and only here we can add duplicate alignments as
  // second best hits). Note that we always keep best align(s) (because they may
  // become second best if we find a better one), but we do not bother storing
  // second best ones if best is still above maxErr threshold!
  if ( e == second_best_[q] && ( best_[q] <= maxErrs_ || best_[q] == second_best_[q] ) ) {
    aligns_[q].push_back(la);
    return;
  }

  // at this point we are guaranteed that our alignment is strictly better
  // than the second best.

  if ( e >= best_[q] ) {
    // This alignment is an improvement of second_best, but not best;

    // all current second best alignments should be removed if we are at the point when we store second best aligns at all
    if ( best_[q] <= maxErrs_  || e==best_[q] ) { // do we need to store second best alignment?
      if ( best_[q] < second_best_[q] ) { // if the old second best was strictly worse than best, it should go
	vec<look_align>::iterator start = aligns_[q].begin();
	aligns_[q].erase(++start, aligns_[q].end());
      }
      aligns_[q].push_back(la); // add new second best (which can be as good as best!)
    }
    // under any circumstances, we do need to save the error count as new second best:
    second_best_[q] = e;

    return;
  }

  // at this point we are guaranteed that our alignment is better
  // than the current best. Two possibilities exist:
  // 1) current best is strictly better than the current second best-
  //    in this case second best should go and current best becomes second best (we definitely store
  //    new best and maybe new second best - former best - depending on how good the new best is);
  // 2) current best is the same as second best, in this case both best and second best
  //    become new second best (if new best is good enough), or they all should go:

  if ( e <= maxErrs_ ) {
    // if e is good enough so that we'll need to keep both best and seconf best:

    if ( second_best_[q] > best_[q] ) {
      // best will be new second best, old second best should go;
      // otherwise (second_best==best), both best and second best are new second_best
      vec<look_align>::iterator start = aligns_[q].begin();
      aligns_[q].erase(++start, aligns_[q].end());
    }
    aligns_[q].insert(aligns_[q].begin(), la);
  } else {
    // e is better than the current best, but still not good; we need to
    // store only new best:
    aligns_[q].clear();
    aligns_[q].push_back(la);
  }
  // update best and second best:
  second_best_[q] = best_[q];
  best_[q] = e;

  ++betterBest_;
}


/// Clean up: remove any look_aligns that are not good enough.
void BestNextBestAlignCollector::Consolidate() {


  for (int i=0; i != aligns_.isize(); ++i) {
    if ( aligns_[i].empty() ) continue;

      // Make sure we have no duplicate alignments
      vec<look_align>::iterator second_best_start = aligns_[i].begin();
      if ( second_best_[i] > best_[i] ) ++second_best_start; // now points to the start of the list of second-best aligns

      // sort second best only, don't touch the best (unless second_best is as good as best, in which case we sort all):
      sort(second_best_start,aligns_[i].end(), order_lookalign_TargetQueryStartEnd());

      // make sure all second bests are unique, erase duplicates:
      aligns_[i].erase(unique(second_best_start,
			      aligns_[i].end(),
			      equal_lookalign_TargetQueryStartEnd() ),
		       aligns_[i].end());
  }
}


/// Add in a look_align (or don't if we don't want it).

void BestAlignCollector::Insert(const look_align & la) {
  int q = la.query_id;
  unsigned int e = la.Errors();


  // if this alignment is not
  // an improvement of at least second_best, nothing needs to be saved.
  if ( e >= second_best_[q] ) return;

  // alignment is worse than best but better than second_best: update second_best
  if ( e > best_[q] && e < second_best_[q] ) {
    second_best_[q] = e;
    return;
  }

  // at this point we are guaranteed that alignment is as good or better than current best;

  if ( e == best_[q] ) {
    // it's another "best" alignment; save it if needed, update second_best and bail out.

    if ( e <= maxErrs_ ) { // if we are already storing the best align:
      // we never allow a duplication of the *best* alignment we keep as long as there is only one;
      if ( la == aligns_[q][0] ) return;

      // important: we are not going to get second_best all wrong despite the possibility
      // to store multiples: at first, there is only *one* best alignment, and the check
      // above guarantees that we are not going to add another instance of the same alignment
      // again. Hence we get to this point only when we find *another* alignment with the
      // same number of errors and setting second_best to e is correct. After this point,
      // adding more "best" aligns with the same error count may result in duplicates,
      // but second_best will remain correct.

      aligns_[q].push_back(la); // add another "best" alignment
    }
    // if we indeed added another best alignment, the second best is indeed the same as best;
    // if e was actually > maxErrs and we did not add anything, the value of second best
    // is irrelevant: we will have first to find better alignment with e < maxErrs which will
    // become best, the current best will become second best, and the current second best will be lost
    second_best_[q] = e;
    return;
  }

  // at this point we are guaranteed that our alignment is strictly better
  // than the best.

  second_best_[q] = best_[q];
  if ( best_[q] <= maxErrs_ ) aligns_[q].clear(); // erase old best aligns if we were keeping any
  best_[q] = e; // reset new best

  if ( e <= maxErrs_ ) aligns_[q].push_back(la); // save the new best align if it's time to save at all

  ++betterBest_;
}


/// Clean up: remove any look_aligns that are not good enough.
void BestAlignCollector::Consolidate() {


  for (int i=0; i != aligns_.isize(); ++i) {
    if ( aligns_[i].empty() ) continue;

      // sort best aligns:
      sort(aligns_[i].begin(),aligns_[i].end(), order_lookalign_TargetQueryStartEnd());

      // make sure all bests are unique, erase duplicates:
      aligns_[i].erase(unique(aligns_[i].begin(),
			      aligns_[i].end(),
			      equal_lookalign_TargetQueryStartEnd() ),
		       aligns_[i].end());
  }
}
