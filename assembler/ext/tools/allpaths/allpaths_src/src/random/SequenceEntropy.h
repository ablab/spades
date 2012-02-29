/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SEQUENCE_ENTROPY_H
#define SEQUENCE_ENTROPY_H

#include "Basevector.h"
#include "CoreTools.h"

// DinukeEntropy: calculate the entropy of the distribution of dinucleotides in a
// given sequence.  This is Sum( p_i log(p_i) ) as p_i varies over the probabilities
// occurring in the dinucleotide distribution.

double DinukeEntropy( const basevector& b );

#endif
