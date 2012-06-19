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

//////////////////////////////////////////////////////////////
// Creators and destructors
/////////////////////////////////////////////////////////////

TightString *newTightString(Coordinate size);

TightString *newTightStringFromString(char *sequence);

void destroyTightString(TightString * tString);

///////////////////////////////////////////////////////////////
// Reading
//////////////////////////////////////////////////////////////

Coordinate getLength(TightString * tString);

char *readTightString(TightString * tString);

Nucleotide getNucleotide(Coordinate nucleotideIndex,
			 TightString * tString);

char getNucleotideChar(Coordinate nucleotideIndex, TightString * tString);

Nucleotide getNucleotideFromString(Coordinate nucleotideIndex,
				   char *string);

char getInverseNucleotideChar(Coordinate nucleotideIndex,
			      TightString * tString);

void readTightStringFragment(TightString * tString, Coordinate start,
			     Coordinate finish, char *string);

Nucleotide charToNucleotide(char c);

///////////////////////////////////////////////////////////////
// Writing
///////////////////////////////////////////////////////////////

void setTightStringLength(TightString * tString, Coordinate length);

void writeNucleotideAtPosition(Nucleotide nucleotide, Coordinate position,
			       TightString * tString);

///////////////////////////////////////////////////////////////
// Array wide operations
///////////////////////////////////////////////////////////////

TightString **newTightStringArrayFromStringArray(char **sequences,
						 IDnum sequenceCount);

TightString **concatenateTightStringArrays(TightString ** array1,
					   TightString ** array2,
					   IDnum size1, IDnum size2);

void destroyTightStringArray(TightString ** array, IDnum arrayLength);

///////////////////////////////////////////////////////////////
// Misc
///////////////////////////////////////////////////////////////

void trimTightString(TightString * tString, Coordinate length);

void trimTightStringArray(TightString ** tStringArray, IDnum arrayLength,
			  Coordinate maxLength);

void trimTightStringArraySanger(TightString ** tStringArray,
				IDnum arrayLength, Coordinate minLength,
				Coordinate maxLength);

void clipTightString(TightString * sequence, Coordinate start,
		     Coordinate finish);

// Exports an array of sequences under FastA format
void exportSequenceArray(char *filename, TightString ** array,
			 IDnum sequenceCount);
void exportTightString(FILE * outfile, TightString * sequence, IDnum index);
#endif
