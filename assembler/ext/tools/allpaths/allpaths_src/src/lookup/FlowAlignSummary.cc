///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "lookup/LookAlign.h"
#include "lookup/FlowAlignSummary.h"

FlowAlignSummary::FlowAlignSummary(const look_align & la):
  name(la.query_id),
  target(la.target_id),
  mis(la.mutations),
  heur(0),
  len(la.query_length),
  start(la.a.pos2()),
  end(la.a.Pos2()),
  bad(false)
{
    pair<int,int> indels = la.a.Gap1Gap2();
    cover = 0;
    for (int i = 0; i != la.a.Nblocks(); ++i) {
      cover += la.a.Lengths(i);
    }
    //ForceAssert(cover != 0);
    ins = indels.second;
    del = indels.first;
}
