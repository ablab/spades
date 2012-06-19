/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

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
#ifndef _TIGHTSTRING_H_
#define _TIGHTSTRING_H_

#include <stdio.h>

typedef unsigned char Codon;

struct tString_st {
	Codon *sequence;
	IDnum length;
}  ATTRIBUTE_PACKED;

//////////////////////////////////////////////////////////////
// Creators and destructors
/////////////////////////////////////////////////////////////

TightString *newTightString(Coordinate size);

void destroyTightString(TightString * tString);

///////////////////////////////////////////////////////////////
// Reading
//////////////////////////////////////////////////////////////

IDnum getLength(TightString * tString);

char *readTightString(TightString * tString);

Nucleotide getNucleotide(Coordinate nucleotideIndex,
			 TightString * tString);

char getNucleotideChar(Coordinate nucleotideIndex, TightString * tString);

void readTightStringFragment(TightString * tString, Coordinate start,
			     Coordinate finish, char *string);

TightString *getTightStringInArray(TightString * tString,
				   IDnum	 position);

///////////////////////////////////////////////////////////////
// Writing
///////////////////////////////////////////////////////////////

void setTightStringLength(TightString * tString, Coordinate length);

void writeNucleotideAtPosition(Nucleotide nucleotide, Coordinate position,
			       TightString * tString);

///////////////////////////////////////////////////////////////
// Array wide operations
///////////////////////////////////////////////////////////////

TightString *newTightStringArrayFromStringArray(char **sequences,
						IDnum sequenceCount,
						char **tSeqMem);

///////////////////////////////////////////////////////////////
// Misc
///////////////////////////////////////////////////////////////

void exportTightString(FILE * outfile, TightString * sequence, IDnum index);
#endif
