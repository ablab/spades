/*
Copyright 2007, 2008, 2009 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#ifndef _KMER_H_
#define _KMER_H_

#include <stdint.h>

#include "globals.h"

void copyKmers(Kmer* k1, Kmer* k2);

void pushNucleotide(Kmer * kmer, Nucleotide nucleotide);
Nucleotide popNucleotide(Kmer * kmer);

int compareKmers(Kmer* k1, Kmer* k2);

void reversePushNucleotide(Kmer * kmer, Nucleotide nucleotide);

KmerKey getKmerKey(Kmer * kmer);

void printKmer(Kmer * kmer);

void clearKmer(Kmer * kmer);

void resetWordFilter(int wordLength);
void resetKeyFilter(int keyLength);

#define KMER_QUOTIENT (MAXKMERLENGTH / 4)
#define KMER_REMAINDER (MAXKMERLENGTH % 4)
#if KMER_REMAINDER
#define KMER_BYTE_SIZE (KMER_QUOTIENT + 1)
#else
#define KMER_BYTE_SIZE KMER_QUOTIENT
#endif
#define KMER_LONGLONGS (KMER_BYTE_SIZE / 8)
#define KMER_LONGS ((KMER_BYTE_SIZE % 8) / 4)
#define KMER_INTS ((KMER_BYTE_SIZE % 4) / 2)
#define KMER_CHARS (KMER_BYTE_SIZE % 2)

struct kmer_st {
#if KMER_LONGLONGS
	uint64_t longlongs[KMER_LONGLONGS];
#endif
#if KMER_LONGS
	uint32_t longs;
#endif
#if KMER_INTS
	uint16_t ints;
#endif
#if KMER_CHARS
	uint8_t chars;
#endif
}  ATTRIBUTE_PACKED;

#endif 
