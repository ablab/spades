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
#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "kmer.h"
#include "utility.h"

static const uint64_t longLongLeftFilter = (uint64_t) 3 << 62; 
static const uint32_t longLeftFilter = (uint32_t) 3 << 30; 
static const uint16_t intLeftFilter = (uint16_t) 3 << 14; 
static const uint8_t charLeftFilter = (uint8_t) 3 << 6; 

#if KMER_LONGLONGS
static uint64_t longLongWordFilter = (uint64_t) ((int64_t) -1); 
#endif
static uint32_t longWordFilter = (uint32_t) ((int32_t) -1); 
static uint16_t intWordFilter = (uint16_t) ((int16_t) -1); 
static uint8_t charWordFilter = (uint8_t) ((int8_t) -1); 

#define UNDEFINED 0
#define CHARS 1
#define INTS 2
#define LONGS 3
#define LONGLONGS 4
static int kmerFilterIndex = UNDEFINED;
static int kmerFilterOffset = 0;
static int kmerFilterLength = 0;
#if KMER_LONGLONGS
static int longLongKmerFilterIndex = KMER_LONGLONGS - 1;
static uint64_t longLongKmerFilter = (uint64_t) ((int64_t) -1); 
#endif

static uint64_t keyFilter = 0;
static int keyFilterIndex = UNDEFINED;
static int keyFilterOffset = 0;
static int keyFilterLength = 0;
#if KMER_LONGLONGS
static int longLongKeyFilterIndex = KMER_LONGLONGS - 1;
#endif

void resetWordFilter(int wordLength) {
	int kmer_bit_size = wordLength * 2;

	if (wordLength > MAXKMERLENGTH) 
		exitErrorf(EXIT_FAILURE, true, "Word length %i greater than max allowed value (%i).\nRecompile Velvet to deal with this word length.", wordLength, MAXKMERLENGTH);

#if KMER_LONGLONGS
	int i;
	for (i = 0; i < KMER_LONGLONGS; i++) {
		if (kmer_bit_size > 64) {
			kmer_bit_size -= 64;
			continue;
		} else if (kmer_bit_size == 64) {
			longLongKmerFilterIndex = i;
			longLongKmerFilter = longLongWordFilter; 
			kmerFilterIndex = LONGLONGS;
			kmerFilterOffset = kmer_bit_size - 2;
			kmerFilterLength = kmer_bit_size;
			longWordFilter = 0;
			intWordFilter = 0;
			charWordFilter = 0;
			return;
		} else {
			longLongKmerFilterIndex = i;
			longLongKmerFilter = (((uint64_t) 1) << kmer_bit_size) - 1;	
			kmerFilterIndex = LONGLONGS;
			kmerFilterOffset = kmer_bit_size - 2;
			kmerFilterLength = kmer_bit_size;
			longWordFilter = 0;
			intWordFilter = 0;
			charWordFilter = 0;
			return;
		}
	}
#endif
#if KMER_LONGS
	if (kmer_bit_size > 32)
		kmer_bit_size -= 32;
	else if (kmer_bit_size == 32) {
		kmerFilterIndex = LONGS;
		kmerFilterOffset = kmer_bit_size - 2;
		kmerFilterLength = kmer_bit_size;
		intWordFilter = 0;
		charWordFilter = 0;
		return;
	} else {
		longWordFilter = (((uint32_t) 1) << kmer_bit_size) - 1;	
		kmerFilterIndex = LONGS;
		kmerFilterOffset = kmer_bit_size - 2;
		kmerFilterLength = kmer_bit_size;
		intWordFilter = 0;
		charWordFilter = 0;
		return;
	}
#endif
#if KMER_INTS
	if (kmer_bit_size > 16)
		kmer_bit_size -= 16;
	else if (kmer_bit_size == 16) {
		kmerFilterIndex = INTS;
		kmerFilterOffset = kmer_bit_size - 2;
		kmerFilterLength = kmer_bit_size;
		charWordFilter = 0;
		return;
	} else {
		intWordFilter = (((uint16_t) 1) << kmer_bit_size) - 1;	
		kmerFilterIndex = INTS;
		kmerFilterOffset = kmer_bit_size - 2;
		kmerFilterLength = kmer_bit_size;
		charWordFilter = 0;
		return;
	}

#endif
#if KMER_CHARS
	if (kmer_bit_size < 8) 
		charWordFilter = (((uint8_t) 1) << kmer_bit_size) - 1;	

	kmerFilterIndex = CHARS;
	kmerFilterOffset = kmer_bit_size - 2;
	kmerFilterLength = kmer_bit_size;
#endif

}

void resetKeyFilter(int keyLength) {
	int kmer_bit_size = keyLength;

	if (keyLength > MAXKMERLENGTH) 
		exitErrorf(EXIT_FAILURE, true, "Key length %i greater than max allowed value (%i).\nRecompile Velvet to deal with this word length.", keyLength, MAXKMERLENGTH);

#if KMER_CHARS
	if (kmerFilterIndex == CHARS) {
		if (kmer_bit_size > kmerFilterLength) {
			kmer_bit_size -= kmerFilterLength;
		} else {
			keyFilterIndex = CHARS;
			keyFilterLength = kmer_bit_size;
			keyFilterOffset = kmerFilterLength - keyFilterLength;
			keyFilter = ((((uint64_t) 1) << keyFilterLength) - 1) << keyFilterOffset;	
			return;
		}
	}
#endif
#if KMER_INTS
	if (kmerFilterIndex == INTS) {
		if (kmer_bit_size > kmerFilterLength) {
			kmer_bit_size -= kmerFilterLength;
		} else {
			keyFilterIndex = INTS;
			keyFilterLength = kmer_bit_size;
			keyFilterOffset = kmerFilterLength - keyFilterLength;
			keyFilter = ((((uint64_t) 1) << keyFilterLength) - 1) << keyFilterOffset;	
			return;
		}
	} else if (intWordFilter) {
		if (kmer_bit_size > 16)
			kmer_bit_size -= 16;
		else {
			keyFilterIndex = INTS;
			keyFilterLength = kmer_bit_size;
			keyFilterOffset = 16 - keyFilterLength;
			keyFilter = ((((uint64_t) 1) << keyFilterLength) - 1) << keyFilterOffset;	
			return;
		}
	}
#endif
#if KMER_LONGS
	if (kmerFilterIndex == LONGS) {
		if (kmer_bit_size > kmerFilterLength) {
			kmer_bit_size -= kmerFilterLength;
		} else {
			keyFilterIndex = LONGS;
			keyFilterLength = kmer_bit_size;
			keyFilterOffset = kmerFilterLength - keyFilterLength;
			keyFilter = ((((uint64_t) 1) << keyFilterLength) - 1) << keyFilterOffset;	
			return;
		}
	} else if (longWordFilter) {
		if (kmer_bit_size > 32)
			kmer_bit_size -= 32;
		else {
			keyFilterIndex = LONGS;
			keyFilterLength = kmer_bit_size;
			keyFilterOffset = 32 - keyFilterLength;
			keyFilter = ((((uint64_t) 1) << keyFilterLength) - 1) << keyFilterOffset;	
			return;
		}
	}
#endif
#if KMER_LONGLONGS
	if (kmerFilterIndex == LONGLONGS) {
		if (kmer_bit_size > kmerFilterLength) {
			kmer_bit_size -= kmerFilterLength;
		} else {
			keyFilterIndex = LONGLONGS;
			longLongKeyFilterIndex = longLongKmerFilterIndex;
			keyFilterLength = kmer_bit_size;
			keyFilterOffset = kmerFilterLength - keyFilterLength;
			keyFilter = ((((uint64_t) 1) << keyFilterLength) - 1) << keyFilterOffset;	
			return;
		}
	} else {
		if (kmer_bit_size > 64)
			abort();
		else {
			keyFilterIndex = LONGLONGS;
			longLongKeyFilterIndex = longLongKmerFilterIndex;
			keyFilterLength = kmer_bit_size;
			keyFilterOffset = 64 - keyFilterLength;
			keyFilter = ((((uint64_t) 1) << keyFilterLength) - 1) << keyFilterOffset;	
			return;
		}
	} 

	if (kmer_bit_size > 64)
		abort();
	else {
		keyFilterIndex = LONGLONGS;
		longLongKeyFilterIndex = longLongKmerFilterIndex - 1;
		keyFilterLength = kmer_bit_size;
		keyFilterOffset = 64 - keyFilterLength;
		keyFilter = ((((uint64_t) 1) << keyFilterLength) - 1) << keyFilterOffset;	
		return;
	}
#endif

}

static void shiftRight(Kmer * kmer) {
	uint64_t leftBits = 0;
	uint64_t rightBits;

#if KMER_CHARS

#if KMER_INTS | KMER_LONGS | KMER_LONGLONGS
	rightBits = kmer->chars & 3;
#endif 

	kmer->chars >>= 2;
	kmer->chars += (uint8_t) leftBits;

#if KMER_INTS | KMER_LONGS | KMER_LONGLONGS
	leftBits = rightBits;
#endif 
#endif

#if KMER_INTS

#if KMER_LONGS | KMER_LONGLONGS
	rightBits = kmer->ints & 3;
#endif

	leftBits <<= 14;
	kmer->ints >>= 2;
	kmer->ints += (uint16_t) leftBits;

#if KMER_LONGS | KMER_LONGLONGS
	leftBits = rightBits;
#endif
#endif

#if KMER_LONGS

#if KMER_LONGLONGS
	rightBits = kmer->longs & 3;
#endif

	leftBits <<= 30;
	kmer->longs >>= 2;
	kmer->longs += (uint32_t) leftBits;

#if KMER_LONGLONGS
	leftBits = rightBits;
#endif
#endif 

#if KMER_LONGLONGS
	int i;
	for (i = KMER_LONGLONGS - 1; i >= 0; i--) {
		rightBits = kmer->longlongs[i] & 3;
		leftBits <<= 62;
		kmer->longlongs[i] >>= 2;
		kmer->longlongs[i] += leftBits;
		leftBits = rightBits;
	}
#endif
}

void copyKmers(Kmer* k1, Kmer* k2) {

#if KMER_LONGLONGS
	int i;
	for (i = 0; i < KMER_LONGLONGS; i++)
		k1->longlongs[i] = k2->longlongs[i];
#endif
#if KMER_LONGS
	k1->longs = k2->longs;
#endif
#if KMER_INTS
	k1->ints = k2->ints;
#endif
#if KMER_CHARS
	k1->chars = k2->chars;
#endif
}

int compareKmers(Kmer* k1, Kmer* k2) {
#if KMER_LONGLONGS
	int i;
#endif

#if KMER_CHARS
	if (k1->chars == k2->chars)
		;
	else if (k1->chars > k2->chars)
		return 1;
	else 
		return -1;
#endif
#if KMER_INTS
	if (k1->ints == k2->ints)
		;
	else if (k1->ints > k2->ints)
		return 1;
	else 
		return -1;
#endif
#if KMER_LONGS
	if (k1->longs == k2->longs)
		;
	else if (k1->longs > k2->longs)
		return 1;
	else 
		return -1;
#endif
#if KMER_LONGLONGS
	for (i = KMER_LONGLONGS - 1; i >= 0; i--) {
		if (k1->longlongs[i] == k2->longlongs[i])
			continue;
		else if (k1->longlongs[i] > k2->longlongs[i])
			return 1;
		else 
			return -1;
	}
#endif

	return 0;
}

void clearKmer(Kmer * kmer) {
#if KMER_LONGLONGS
	int i;
	for (i = 0; i < KMER_LONGLONGS; i++)
		kmer->longlongs[i] = 0;
#endif
#if KMER_LONGS
	kmer->longs = 0;
#endif
#if KMER_INTS
	kmer->ints = 0;
#endif
#if KMER_CHARS
	kmer->chars = 0;
#endif
}

void printKmer(Kmer * kmer) {

#if KMER_CHARS
	velvetLog("%hx\t", kmer->chars);
#endif
#if KMER_INTS
	velvetLog("%x\t", kmer->ints);
#endif
#if KMER_LONGS
	velvetLog("%x\t", kmer->longs);
#endif
#if KMER_LONGLONGS
	int i;
	for (i = KMER_LONGLONGS - 1; i >= 0; i--)
		velvetLog("%llx\t", (long long) kmer->longlongs[i]);
#endif
	velvetLog("\n");
}

void pushNucleotide(Kmer * kmer, Nucleotide nucleotide) {

#if KMER_LONGLONGS
	uint64_t * ptr;
	int i;
#endif
#if KMER_LONGLONGS > 1 | KMER_LONGS | KMER_INTS | KMER_CHARS
	uint64_t leftBits; 
#endif
	uint64_t rightBits = 0;

#if KMER_LONGLONGS
	ptr = kmer->longlongs;

#if KMER_LONGLONGS > 1
	for (i = 0; i < longLongKmerFilterIndex; i++) {
		leftBits = (*ptr & longLongLeftFilter);
		leftBits >>= 62;
		*ptr <<= 2;
		*ptr += rightBits;
		*ptr &= longLongWordFilter;
		rightBits = leftBits;
		ptr++;
	}
#endif

#if KMER_LONGS | KMER_INTS | KMER_CHARS
	leftBits = (*ptr & longLongLeftFilter);
	leftBits >>= 62;
#endif

	*ptr <<= 2;
	*ptr += rightBits;
	*ptr &= longLongKmerFilter;

#if KMER_LONGS | KMER_INTS | KMER_CHARS
	rightBits = leftBits;
#endif
#endif

#if KMER_LONGS

#if KMER_INTS | KMER_CHARS
	leftBits = kmer->longs & longLeftFilter;
	leftBits >>= 30;
#endif
	kmer->longs <<= 2;
	kmer->longs += rightBits;
	kmer->longs &= longWordFilter;

#if KMER_INTS | KMER_CHARS
	rightBits = leftBits;
#endif

#endif

#if KMER_INTS

#if KMER_CHARS
	leftBits = kmer->ints & intLeftFilter;
	leftBits >>= 14;
#endif 
	kmer->ints <<= 2;
	kmer->ints += rightBits;
	kmer->ints &= intWordFilter;

#if KMER_CHARS
	rightBits = leftBits;
#endif 

#endif

#if KMER_CHARS
	kmer->chars <<= 2;
	kmer->chars += rightBits;
	kmer->chars &= charWordFilter;
#endif

#if KMER_LONGLONGS
	kmer->longlongs[0] += nucleotide;
	if (kmer->longlongs[0] >= nucleotide)
		return;

	for (i = 1; i < KMER_LONGLONGS; i++)
		if (++kmer->longlongs[i])
			return;
#if KMER_LONGS
	if (++kmer->longs)
		return;
#endif
#if KMER_INTS
	if (++kmer->ints)
		return;
#endif
#if KMER_CHARS
	++kmer->chars;
#endif

#else

#if KMER_LONGS
	kmer->longs += nucleotide;
	if (kmer->longs >= nucleotide)
		return;
#if KMER_INTS
	if (++kmer->ints)
		return;
#endif
#if KMER_CHARS
	++kmer->chars;
#endif

#else

#if KMER_INTS
	kmer->ints += nucleotide;
	if (kmer->ints >= nucleotide)
		return;
#if KMER_CHARS
	++kmer->chars;
#endif

#else 

#if KMER_CHARS
	kmer->chars += nucleotide;
#endif

#endif
#endif
#endif
}

Nucleotide popNucleotide(Kmer * kmer) {
	Nucleotide nucl;

#if KMER_LONGLONGS
	nucl = kmer->longlongs[0] & 3;
#elif KMER_LONGS
	nucl = kmer->longs & 3;
#elif KMER_INTS
	nucl = kmer->ints & 3;
#elif KMER_CHARS
	nucl = kmer->chars & 3;
#endif

	shiftRight(kmer);
	return nucl;
}

void reversePushNucleotide(Kmer * kmer, Nucleotide nucleotide) {
	uint64_t templongLong = nucleotide;

	shiftRight(kmer);

	switch(kmerFilterIndex) {
		case UNDEFINED:
			abort();
#if KMER_LONGLONGS
		case LONGLONGS:
			kmer->longlongs[longLongKmerFilterIndex] += templongLong << kmerFilterOffset; 			
			return;
#endif
#if KMER_LONGS
		case LONGS:
			kmer->longs += templongLong << kmerFilterOffset; 			
			return;
#endif
#if KMER_INTS
		case INTS:
			kmer->ints += templongLong << kmerFilterOffset; 			
			return;
#endif
#if KMER_CHARS
		case CHARS:
			kmer->chars += templongLong << kmerFilterOffset; 			
			return;
#endif
	}		

	exitErrorf(EXIT_FAILURE, true, "Anomaly in k-mer filering");
}

KmerKey getKmerKey(Kmer * kmer) {
	KmerKey key = 0;	

#if KMER_CHARS
	if (keyFilterIndex == CHARS) {
		key = kmer->chars & keyFilter;
		key >>= keyFilterOffset;
		return key;
	}
#endif

#if KMER_INTS
	if (keyFilterIndex == INTS) {
#if KMER_CHARS
		key += kmer->chars;
#endif
		key <<= keyFilterLength;
		key += (kmer->ints & keyFilter) >> keyFilterOffset;
		return key;
	}
#endif

#if KMER_LONGS
	if (keyFilterIndex == LONGS) {
#if KMER_CHARS
		key += kmer->chars;
#endif
#if KMER_INTS
		key <<= 16;
		key += kmer->ints;
#endif
		key <<= keyFilterLength;
		key += (kmer->longs & keyFilter) >> keyFilterOffset;
		return key;
	}
#endif

#if KMER_LONGLONGS
	if (keyFilterIndex == LONGLONGS) {
#if KMER_CHARS
		key += kmer->chars;
#endif
#if KMER_INTS
		key <<= 16;
		key += kmer->ints;
#endif
#if KMER_LONGS
		key <<= 32;
		key += kmer->longs;
#endif
		if (longLongKeyFilterIndex != longLongKmerFilterIndex) {
			key <<= kmerFilterLength;
			key += kmer->longlongs[longLongKmerFilterIndex];
		}
		key <<= keyFilterLength;
		key += (kmer->longlongs[longLongKeyFilterIndex] & keyFilter) >> keyFilterOffset;

		return key;
	}
#endif
	abort();
	return 0;
} 
