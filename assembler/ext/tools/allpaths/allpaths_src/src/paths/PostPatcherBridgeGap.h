///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PostPatcherBridgeGap.  Utilities for briding gaps for PostPatcher.

bool
BridgeGap(const vec<basevector> & bridges,
	   const vec<int> &gaps,
	   int pgap,	// predicted gap
	   int pdev,	// sigma of predicted gap
	   float max_sigma, // within how many sigma must we be?
	   int flags,	// flags to control method
	   const bool verbose,
	   fastavector &bridge);

// Flags
#define BRIDGEGAP_ONESIZE  1	/* only bridge if all candidates are same size */
#define BRIDGEGAP_BESTFIT  2	/* use those which most closely match predicted gap */
#define BRIDGEGAP_MAXCOUNT 4	/* use bridge length with most candidates */

