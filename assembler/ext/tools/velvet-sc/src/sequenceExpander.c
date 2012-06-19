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
#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "tightString.h"
#include "passageMarker.h"

static TightString *readPositivePassageMarker(PassageMarker * marker,
					      TightString ** seqs,
					      int WORDLENGTH)
{
	Coordinate index;
	Nucleotide nucleotide;
	TightString *tString =
	    seqs[getPassageMarkerSequenceID(marker) - 1];
	TightString *res = newTightString(getPassageMarkerLength(marker));

	for (index = 0; index < getLength(tString); index++) {
		nucleotide =
		    getNucleotide(getPassageMarkerStart(marker) + index +
				  WORDLENGTH - 1, tString);
		writeNucleotideAtPosition(nucleotide, index, res);
	}

	return res;
}

static TightString *readNegativePassageMarker(PassageMarker * marker,
					      TightString ** seqs)
{
	Coordinate index;
	Nucleotide nucleotide;
	TightString *tString =
	    seqs[getAbsolutePassMarkerSeqID(marker) - 1];
	TightString *res = newTightString(getPassageMarkerLength(marker));

	for (index = 0; index < getPassageMarkerLength(marker); index++) {
		nucleotide =
		    getNucleotide(getPassageMarkerStart(marker) - index,
				  tString);
#ifndef COLOR
		writeNucleotideAtPosition(3 - nucleotide, index, res);
#else
		writeNucleotideAtPosition(nucleotide, index, res);
#endif
	}

	return res;
}

TightString *expandPassageMarker(PassageMarker * marker,
				 TightString ** sequences, int WORDLENGTH)
{
	if (getPassageMarkerSequenceID(marker) > 0)
		return readPositivePassageMarker(marker, sequences,
						 WORDLENGTH);
	else
		return readNegativePassageMarker(marker, sequences);
}
