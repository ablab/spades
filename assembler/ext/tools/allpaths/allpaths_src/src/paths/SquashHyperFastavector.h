///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SQUASH_HYPER_FASTAVECTOR_H
#define SQUASH_HYPER_FASTAVECTOR_H

#include "paths/HyperFastavector.h"

/**
 * SquashHyperFastavector
 *
 * Look for and destroy scaffolds that either are duplicates of, or
 * are strictly contained in, other scaffolds.
 *
 * work_dir: full path name (several temp files will be saved in here)
 * num_threads: used by ConvertToSuperbs
 * max_error_rate: max error to accept the alignmnt between two scaffolds
 * hfv: input and output
 * plog: optional log
 */
void SquashHyperFastavector( const String &work_dir,
			     const int num_threads,
			     const float max_error_rate,
			     HyperFastavector &hfv,
			     ostream *plog = 0 );

#endif
