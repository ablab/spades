///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// Class to store all the information we want to keep from a lookalign
/// and access it easily.

#include "system/Types.h"

class look_align;

struct FlowAlignSummary {
  FlowAlignSummary(int name = 0, int target= 0, int mis= 0, int ins= 0, 
		   int del = 0, int heur= 0, int len= 0, int cover = 0,
		   int start=0, int end=0, bool bad = false):
    name(name),target(target),mis(mis),ins(ins),del(del),heur(heur),len(len),
    cover(cover), start(start), end(end), bad(bad)
  {}

  ///Constructor from a look_align.
  /// This does not have info to set heur or bad, but leaves them as 0.
  FlowAlignSummary(const look_align & la);

  int name;
  int target;
  int mis;
  int ins;
  int del;
  int heur;
  int len;
  int cover;
  int start; // start of alignment in target
  int end; // end of alignment in target
  Bool bad; //does not pass the heuristics test.
  int goodBases() const { return cover - mis - ins - del; }
  float score() const { return len ? (cover - mis) / float(len) : 0.0; }
  int coverScorePercent() const { 
    return cover ? (cover - mis) * 100 /cover  : 0; 
  }
};

inline bool 
operator<(const FlowAlignSummary & lhs, const FlowAlignSummary &rhs) {
  return lhs.score() < rhs.score();
}

inline bool 
operator>(const FlowAlignSummary & lhs, const FlowAlignSummary &rhs) {
  return lhs.score() > rhs.score();
}
