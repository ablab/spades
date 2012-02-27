/***************************************************************************
 * Title:          GFFEntry.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef GFF_ENTRY_H
#define GFF_ENTRY_H


/*
  Field number	Description	Comments

Description from: http://www.drive5.com/pals/pals_userguide.htm

chr3  pals  hit  10863665  10864065  310  +  .  Target  chr1  17255 17657;  maxe  0.058

1 Query sequence name	PALS truncates the FASTA annotation at the first blank character.
2 Source                Is always set to pals.
3 Feature               Is always set to hit.
4 Start                 Position within the sequence of hit start. 
                        The first base is at position 1 (GFF uses 1-based coordinates).
5 End                   Position of hit end.
6 Score                 Alignment score. Matches score +1, differences
                        (substitutions, deletions and insertions) score ?Ð3.
7 Strand                Is + if the aligned regions are on the same
                        strand, Ð if the two regions are reverse
                        complemented with respect to each other. 
8 Frame                 Is always set to ".", a period.
9 Attributes            Formatted as follows:  Target name start end ; maxe e

Fields are separated by spaces (not tabs). Name is the name of the
target sequence, start and end are the coordinates of the matching
region, and e is an upper bound on the error ratio, i.e. the number of
differences between the two regions divided by the length of the
region. A difference is a substitution, deletion or insertion of
length one. E.g., e=0.058 means identity ³ 94.2%. 
*/

#include <string>
#include <istream>

class GFFEntry {
public:
  std::string refSeq;
  std::string source;
  std::string feature;
  ssize_t refStart, refEnd;
  double score;
  char strand;
  char frame;
  std::string qrySeq;
  ssize_t qryStart, qryEnd;
  double errorRatio;
  bool ParsePALSLine(std::istream& in);
};

#endif
