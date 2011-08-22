/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

#ifndef _TRANSLATE_H_
#define _TRANSLATE_H_

#include "io_lib/scf.h"
#include "io_lib/Read.h"
#include "io_lib/expFileIO.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Translates an Scf structure into a Read structure.
 * The Scf structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Read structure upon success.
 *    NULLRead upon failure.
 */
Read *scf2read(Scf *scf);

/*
 * Translates a Read structure into a Scf structure.
 * The Read structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Scf structure upon success.
 *    NULL upon failure.
 */
Scf *read2scf(Read *read);

/*
 * Translates a Read structure and an Experiment file.
 * The Read structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Exp_info structure upon success.
 *    NULL upon failure (FIXME: need to free memory here)
 */
Exp_info *read2exp(Read *read, char *EN);

/*
 * Controls the use of the SQ and ON lines when loading an experiment file.
 * The default (value&1 == 1) is to load these into the Read structure.
 * With value&1 == 0 we load the sequence directly from the trace file
 * (LT line).
 * value&2 controls whether to use the SL/SR fields when setting the cutoff.
 * value&2 == 0 implies to do so, and value&2 == 2 implies to not.
 *
 * Returns:
 *    The old value.
 */
int read_experiment_redirect(int value);

/*
 * Takes an original read structure and a set of edit change arrays and
 * produces a new base position array incorporating all the edits. For
 * insertions, interpolation is used to derive a suitable sample position.
 *
 * INPUTS:
 *
 * Read   *r       = The original unedited read structure
 * int     Comp    = 0=Normal sequence, 1=Complemented sequence
 * int     Ned     = Length of edited arrays to follow
 * char   *edBases = Sequence of base characters incorporating ins/del edits
 * uint_2 *edPos   = Corresponding original base numbers, 0 indicates an
 *		     insertion. Base numbers start at 1.
 *
 * OUTPUTS:
 *
 * This array is assumed to be empty with an allocated length of Ned elements.
 *
 * uint_2* basePos = Base positions in samples
 */

void read_update_base_positions( Read *r, int Comp, int Ned, char *edBases,
				 int_2 *edPos, uint_2 *basePos );

/*
 * Takes a set of edit change arrays and produces a new set of confidence
 * arrays incorporating all the edits.
 *
 * INPUTS:
 *
 * int    Ned     = Length of edited arrays to follow
 * char*  edBases = Sequence of base characters incorporating ins/del edits
 * int1*  edConf  = Corresponding confidence values, 100 for insertions
 *
 *
 * OUTPUTS:
 *
 * These output arrays are assumed to be empty with an allocated length
 * of Ned elements each. The names and types are identical to the same
 * elements in the Read structure.
 *
 * char*  prob_A  = Base confidence A
 * char*  prob_C  = Base confidence C
 * char*  prob_G  = Base confidence G
 * char*  prob_T  = Base confidence T
 *
 */

void read_update_confidence_values( int Ned, char* edBases, int1* edConf,
                                    char* prob_A, char* prob_C, char* prob_G, char* prob_T );


/*
 * Translates an experiment file to a Read structure.
 * The Exp_Info structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Read structure upon success.
 *    NULLRead upon failure.
 */
Read *exp2read(Exp_info *e, char *fn);

#ifdef __cplusplus
}
#endif

#endif
