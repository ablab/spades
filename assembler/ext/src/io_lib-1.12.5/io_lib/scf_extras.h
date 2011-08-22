/*
 * Copyright (c) Medical Research Council 1998. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies and that credit is given
 * where due.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * Kathryn Beal, as part of the Staden Package at the MRC Laboratory of
 * Molecular Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

#ifndef _SCFBITS_H_
#define _SCFBITS_H_

#include "io_lib/expFileIO.h"

#ifdef __cplusplus
extern "C" {
#endif

int get_read_conf(Exp_info *e, int length, int2 *opos, int1 *conf);

#ifdef __cplusplus
}
#endif

#endif /* _SCFBITS_H_ */
