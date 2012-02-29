///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PostPatcherBridge.  Utilities for briding gaps for PostPatcher.

#include "Basevector.h"
#include "Fastavector.h"
#include "MainTools.h"
#include "paths/HyperFastavector.h"
#include "paths/PostPatcherBridgeGap.h"

// Choose a bridge.  For a given gap, if there is more than one patch, we
// require that the patches have only mismatches between each other, and the
// total mismatch fraction is at most 5%.  If there are multiple gaps, choose 
// the one with the most possible bridges, or if equal, the one closest to 
// predicted gap.  Return true if we have something.  Experimental.  --bruce

bool
BridgeGap(const vec<basevector> & bridges,
	  const vec<int> &gaps,
	  int pgap,	// predicted gap
	  int pdev,	// sigma of predicted gap
	  float max_sigma, // within how many sigma must we be?
	  int flags,		// flags to control algorithm
	  const bool verbose,
	  fastavector &bridge)
{
  bridge.clear();
  u_int nbridges = bridges.size();
  vec<int> sizes(nbridges);
  vec<int> deltas(nbridges);
  vec<int> id(nbridges);
  for (u_int i = 0; i < nbridges; ++i) {
    sizes[i] = bridges[i].size();
    deltas[i] = gaps[i] - pgap;
    id[i] = i;
  }

  if (flags & BRIDGEGAP_ONESIZE) {
    // ensure all candidates are same size
    UniqueSort(sizes);
    if (sizes.size() > 1) {
      if (verbose)
	cout << "Rejecting patches because more than one size." << endl;
      return false;
    }
  }

  // Select a winning gap (delta): if there is a gap value which
  // occurs more than any other, then use it.  If there are several
  // ties, use the one which is closest to the predicted gap.
  int delta = deltas[0];
  if (deltas.size() > 1) {
    SortSync(deltas, id);
    // Count how many of each delta value
    vec<int> udeltas, ucounts;
    for (u_int i = 0; i < deltas.size(); ++i) {
      if (i == 0 || udeltas.back() != deltas[i]) {
	udeltas.push_back(deltas[i]);
	ucounts.push_back(1);
      } else {
	++ucounts.back();
      }
    }

    if (flags & BRIDGEGAP_MAXCOUNT) {
      // Trim to those which share largest count
      ReverseSortSync(ucounts, udeltas);
      for (u_int i = 1; i < ucounts.size(); ++i) {
	if (ucounts[i] != ucounts[0]) {
	  udeltas.resize(i);
	  ucounts.resize(i);
	  break;
	}
      }
    }

    if (flags & BRIDGEGAP_BESTFIT) {
      // If there is more than one, use one closest to predicted gap
      if (ucounts.size() > 1) {
	vec<int> adeltas(udeltas.size());
	for (u_int i = 0; i < udeltas.size(); ++i)
	  adeltas[i] = abs(udeltas[i]);
	SortSync(adeltas, udeltas);
      }
    }
    // The winner
    delta = udeltas[0];
  }

  float zscore = float(abs(delta)) / float(pdev);
  if (max_sigma > 0 && zscore > max_sigma) {
    if (verbose)
      cout << "Rejecting patches because gap is " << ToString(zscore,2) << " sigma from predicted gap." << endl;
    return false;
  }

  vec<int> used;
  for (u_int i = 0; i < deltas.size(); ++i) {
    int j = id[i];
    if (deltas[i] == delta) {
      if (used.size() == 0)
	bridge = fastavector(bridges[j]);
      else
	bridge.combine(fastavector(bridges[j]));
      used.push_back(j);
    }
  }

  if (verbose) {
    cout << endl;
    String label = "winning_bridge";
    for (u_int i = 0; i < used.size(); ++i)
      label += "_" + ToString(used[i]+1);
    bridge.Print(cout, label);
  }

  int amb = bridge.AmbCount( );
  const double max_amb_frac = 0.05;
  const int max_amb = 15;
  if ( amb > max_amb && amb > int( floor( double( bridge.size( ) ) * max_amb_frac ) ) )
    {    if (verbose)
	{    cout << "Rejecting patches because too much ambiguity."
		  << endl;    }
      return false;    }

  return true;
}










