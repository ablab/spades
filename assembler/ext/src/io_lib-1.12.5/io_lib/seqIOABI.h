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

#ifndef _seqIOABI_h_
#define _seqIOABI_h_

#include <sys/types.h> /* off_t */
#include "io_lib/os.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * The ABI magic number - "ABIF"
 */
#define ABI_MAGIC	((int_4) ((((('A'<<8)+'B')<<8)+'I')<<8)+'F')

/*
 * The index is located towards the end of the ABI trace file.
 * It's location is given by a longword at a fixed place.
 */
#define IndexPO ((off_t)26)

#define IndexEntryLength 28

/*
 * Here are some labels we will be looking for, four chars packed
 * into an int_4
 */
#define LABEL(a) ((int_4) ((((((a)[0]<<8)+(a)[1])<<8)+(a)[2])<<8)+(a)[3])
#define DataEntryLabel    LABEL("DATA")
#define BaseEntryLabel    LABEL("PBAS")
#define BasePosEntryLabel LABEL("PLOC")
#define SpacingEntryLabel LABEL("SPAC")
#define SignalEntryLabel  LABEL("S/N%")
#define FWO_Label         LABEL("FWO_")
#define MCHNLabel         LABEL("MCHN")
#define PDMFLabel         LABEL("PDMF")
#define SMPLLabel         LABEL("SMPL")
#define PPOSLabel         LABEL("PPOS")
#define CMNTLabel         LABEL("CMNT")
#define GelNameLabel      LABEL("GELN")
#define LANELabel         LABEL("LANE")
#define RUNDLabel         LABEL("RUND")
#define RUNTLabel         LABEL("RUNT")
#define MTXFLabel         LABEL("MTXF")
#define SPACLabel         LABEL("SPAC")
#define SVERLabel         LABEL("SVER")
#define MODLLabel         LABEL("MODL")
#define BaseConfLabel	  LABEL("PCON")


/*
 * From the ABI results file connected to `fp' whose index starts
 * at byte offset `indexO', return in `val' the `lw'th long word
 * from the `count'th entry labelled `label'.
 * The result is 0 for failure, or index offset for success.
 */
int getABIIndexEntryLW(FILE *fp, off_t indexO,
		       uint_4 label, uint_4 count, int lw,
		       uint_4 *val);

/*
 * From the ABI results file connected to `fp' whose index starts
 * at byte offset `indexO', return in `val' the `sw'th short word
 * from the `count'th entry labelled `label'.
 * The result is 0 for failure, or index offset for success.
 */
int getABIIndexEntrySW(FILE *fp, off_t indexO,
		       uint_4 label, uint_4 count, int sw,
		       uint_2 *val);

/*
 * Gets the offset of the ABI index.
 * Returns -1 for failure, 0 for success.
 */
int getABIIndexOffset(FILE *fp, uint_4 *indexO);

/*
 * Get an "ABI String". These strings are either pointed to by the index
 * offset, or held in the offset itself when the string is <= 4 characters.
 * The first byte of the string determines its length.
 * 'string' is a buffer 256 characters long.
 *
 * Returns -1 for failure, string length for success.
 */
int getABIString(FILE *fp, off_t indexO, uint_4 label, uint_4 count,
		 char *string);

/*
 * Get an "ABI Int_1". This is raw 1-byte integer data pointed to by the
 * offset, or held in the offset itself when the data is <= 4 characters.
 *
 * If indexO is 0 then we do not search for (or indeed use) label and count,
 * but simply assume that we are already at the correct offset and read from
 * here. (NB: This negates the length <= 4 check.)
 *
 * Returns -1 for failure, length desired for success (it'll only fill out
 * up to max_data_len elements, but it gives an indication of whether there
 * was more to come).
 */
int getABIint1(FILE *fp, off_t indexO, uint_4 label, uint_4 count,
	       uint_1 *data, int max_data_len);

/*
 * Get an "ABI Int_2". This is raw 2-byte integer data pointed to by the
 * offset, or held in the offset itself when the data is <= 4 characters.
 *
 * Returns -1 for failure, length desired for success (it'll only fill out
 * up to max_data_len elements, but it gives an indication of whether there
 * was more to come).
 */
int getABIint2(FILE *fp, off_t indexO, uint_4 label, uint_4 count,
	       uint_2 *data, int max_data_len);

/*
 * Get an "ABI Int_4". This is raw 4-byte integer data pointed to by the
 * offset, or held in the offset itself when the data is <= 4 characters.
 *
 * Returns -1 for failure, length desired for success (it'll only fill out
 * up to max_data_len elements, but it gives an indication of whether there
 * was more to come).
 */
int getABIint4(FILE *fp, off_t indexO, uint_4 label, uint_4 count,
	       uint_4 *data, int max_data_len);

int dump_labels(FILE *fp, off_t indexO);

/*
 * Change the DATA counts for fetching traces
 */
void abi_set_data_counts(int f, int w, int o, int _);

/*
 * Put the DATA counts back to their defaults.
 */
void abi_reset_data_counts(void);

#ifdef __cplusplus
}
#endif

#endif /* _seqIOABI_h_ */
