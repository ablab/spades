/***************************************************************************
 * Title:          SAP.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SAP.h"


//  Each possible variant is represented

void CreateGrid(ssize_t dim, Grid& grid) {
  grid = new Column[dim+1];
  ssize_t d;
  for (d = 0; d <= dim; d++ ) {
    grid[d] = new Cell[4];
  }
}

void DeleteGrid(ssize_t dim, Grid &grid) {
  ssize_t d;
  for (d = 0; d <= dim; d++ ) {
    delete[] grid[d];
  }
  delete[] grid;
}



ssize_t FindSolidPosition(DNASequence &seq, 
											//CountedIntegralTupleList &spectrum, 
											CountedIntegralTupleDict &spectrum,
											ssize_t span,
											ssize_t &lastSolidPos, ssize_t& pos) {

	CountedIntegralTuple tuple;
  ssize_t solidSpanFound = 0;
  lastSolidPos = -1;
  for (pos = 0; pos < seq.length - IntegralTuple::tupleSize + 1; pos++) {
		
		ssize_t isValid = tuple.StringToTuple(&seq.seq[pos]);
	//		tuple.assign((char*) &(seq.seq[pos]));

		if (isValid) {
      if (spectrum.DictLookupBinaryTuple(tuple) != -1) {
				// found a solid tuple
				if (lastSolidPos == -1) {
					lastSolidPos = pos;
				}
				else {
					// if the span is long enough, return this position
					if (pos - lastSolidPos >= span) {
						solidSpanFound = 1;
					}
				}
      }
      else {
				// found an erroneous tuple

				if (lastSolidPos > -1 and pos - lastSolidPos - 1 >= span) {
					// If this tuple is erroneous, but we are currently in a valid span
					// of solid tuples, return the previous position, as this does not have an
					// error.
					pos -= 1;
					return 1;
				}
				else {
					//	  std::cout << "found an error at " << pos << " and last solid at: " << lastSolidPos << std::endl;
					lastSolidPos = -1;
				}
      }
    } else {
      // found an invalid (masked, etc.) tuple, keep looking
      if (lastSolidPos > -1 and pos - lastSolidPos -1 >= span) {
				pos -=1;
				return 1;
      }
      else 
				lastSolidPos = -1;
    }
  }
	pos--;
  return solidSpanFound;
}


ssize_t SolidifyUntilFixed(DNASequence &seq, 
											 CountedIntegralTupleList &origSeqEnum, ssize_t origEnumPos,
											 //CountedIntegralTupleList &spectrum, 
											 CountedIntegralTupleDict &spectrum, 
											 IntMatrix &scoreMat,
											 FixParams &params, DNASequence &fixedSeq, 
											 ssize_t &replaceLength, Stats &stats, ssize_t &firstEdit, ssize_t &lastSolidifiedPos) {

  stats.Reset();
  // the reads starts having errors after solidpos. fix them.
  ssize_t fixed;
  ssize_t pos; // for iterating over positions
  ssize_t k; // for iterating over gaps
  ssize_t n; // for iterating over nucleotides
  fixed = 0;
  ssize_t i;
  replaceLength = 0;
  Grid grid;

  // start at 'pos' and fix until a tuple is found that corresponds
  // to the read
  ssize_t init = 0;
  SEdge edge;
	CountedIntegralTuple tuple;
	tuple.StringToTuple(seq.seq);

  pos = 0 ;
  Cube cube;
  ssize_t solidExtensionFound = 1;
  CountedIntegralTuple newTuple;
	//  newTuple.reserve(spectrum.tupleSize);
  while (pos < seq.length - IntegralTuple::tupleSize + 1 and 
				 fixed == 0 and
				 solidExtensionFound) {

    solidExtensionFound = 0;
    CreateGrid(params.maxGap, grid);
    cube.push_back(grid);

		// Try and extend each gap level in the grid.
		for (k = 0; k < params.maxGap; k++) {
			// foreach nucleotide
			if (init == 0) {
				// we do not want to compute values for all 4 nucleotides
				// on the first iteration.

				// this is the first iteration, so we need to initialize
				// the only edge.  We assume that this points to nothing.
				solidExtensionFound    = 1;

				edge.prevNuc      = -1;
				edge.prevLevel    = -1;
				// TODO: sparc and mac compilers warn prevLevel is unsigned but is
				// set to -1 here; it appears safe to do this because
				// (i) an edge is ignored when prevNuc==-1 (w/o even checking
				// other fields), and
				// (ii) prevNuc and prevLevel are set to values != -1 simultaneously.
				//  But may need to check this further.
				edge.prevPosition = -1;
				edge.score = 0;

				//if (!GetHashValue(seq, pos, spectrum.tupleSize, tuple)) {
				ssize_t isValid = tuple.StringToTuple(&seq.seq[pos]);
				if (!isValid) {
					std::cout << "Error, can only solidify starting on valid tuples (internal error - please" 
										<< std::endl << "contact the authors)." << std::endl;
					exit(1);
				}
				edge.tuple = tuple;
				init = 1;

				// store the edge here.
				//				grid[0][numeric_nuc_index[(unsigned char) (edge.tuple.tuple & 3)]][edge.tuple] = edge;
				grid[0][numeric_nuc_index[(unsigned char) edge.tuple.GetNucleotide0()]][edge.tuple] = edge;

			}
			else {
				// assign cell values for each nucleotide
				for ( n = 0; n < 4; n++ ) {
					// try to find incoming edges to each nucleotide.

					EdgeIterator edgeIt, edgeEnd;
					ssize_t m;
					unsigned char prevNuc;
					ssize_t score;


					// Search for paths that use deletions.  Deletions are found by 
					// searching for transitions more than 1 layer back in the grid.
	   
					for (m = 1; m <= std::min(params.maxGap, pos); m++ ) {
	     
						// reference cell at position [pos -solidPos]
						//                   level k
						//                   nucleotide n
						// and find all edges going into it

						// Look through all tuples at position pos - m
						// for a transition involving a deletion or deletion followed by a mutation.
						
						// m == 1 implies mutation
						// m > 1  searches deletions.
						for (prevNuc = 0; prevNuc < 4 ; prevNuc++) {
							// Look at all edges stored at the previous nucleotide 
							// to see if the k-1 tuple overlaps with a 
							for (edgeIt = cube[pos  - m][k][prevNuc].begin();
									 edgeIt != cube[pos - m][k][prevNuc].end();
									 edgeIt++) {
								ForwardNuc(edgeIt->first, n, newTuple);
								if (spectrum.DictLookupBinaryTuple(newTuple) != -1) {
									// create a new edge here.
									// compare the nucleotide at 'pos' with what we are 
									// changing to at 'n'
									//									if (unmasked_nuc_index[seq.seq[pos + IntegralTuple::tupleSize-1]] < 4) 
									score = scoreMat[n][unmasked_nuc_index[seq.seq[pos+IntegralTuple::tupleSize-1]]] + 
										edgeIt->second.score + params.gapExtend * (m-1);
									/*									else {
										score = 1 + edgeIt->second.score + params.gapExtend * (m-1);
									}
									*/
									// Store a new edge if this score is better than
									// all other edges
									// and edge is defined as a 'prev tuple' 'cur tuple' edge
									// for each nucleotide.
		   
									// so we need to pass in the set at the current position
									// the current level 'm', the current nucleotide 'n'
									// and try to connect it to the tuple 'newTuple'
									if (score < params.scoreThreshold) {
										solidExtensionFound = 1;
										
										if (StoreEdge(pos  - m, k, prevNuc, edgeIt->first,
																	cube[pos][k][n], newTuple, score, 
																	edgeIt->second.score, edgeIt->second.solidStretch)) {
											// A new edge was stored.
											/*											std::string st; newTuple.ToString(st);
											std::cout << "new ext: " << st << std::endl;
											*/
										}
										if (cube[pos][k][n][newTuple].solidStretch >= params.span) {
											// A stretch of nucleotides has been found where the 
											// no modifications have been made.
											if (
													newTuple == origSeqEnum[pos + origEnumPos] // compare tuples
													//													newTuple.tuple == origSeqEnum[pos + origEnumPos].tuple
													and
													origSeqEnum[pos + origEnumPos].count > 0) {
												// No changes have been made, and the fixed sequence is the same
												// as the original sequence at this position, so we can stop trying to fix
												// it for now.
												fixed = 1;
											}
										}
										stats.numEdge++;
									}
								}
							}
						}
					} // done fixing deletions
					// Look for insertions.  These are transitions from the previous level
					// in the grid.
					if (k > 0) {
						for (prevNuc = 0; prevNuc < 4; prevNuc++ ) {
							for (edgeIt = cube[pos ][k-1][prevNuc].begin();
									 edgeIt != cube[pos][k-1][prevNuc].end();
									 ++edgeIt) {
								ForwardNuc(edgeIt->first, n, newTuple);
								if (spectrum.DictLookupBinaryTuple(newTuple) != -1) {
									// Only allow insertions before regular nucleotides
									if (unmasked_nuc_index[seq.seq[pos+ IntegralTuple::tupleSize-1]] < 4) 
										score = edgeIt->second.score + params.gapExtend;
									/*									score = ScoreMat[n][unmasked_nuc_index[seq.seq[pos+IntegralTuple::tupleSize-1]]] + 
										edgeIt->second.score + params.gapExtend;
									*/
									else
										score = 10000 + edgeIt->second.score + params.gapExtend;

									if (score < params.scoreThreshold) {
										solidExtensionFound = 1;
										StoreEdge(pos , k-1, prevNuc, edgeIt->first,
															cube[pos][k][n], newTuple, score, 
															edgeIt->second.score, edgeIt->second.solidStretch );
										stats.numEdge++;
									}
								}
							}
						}
					}
				}
      }
    }
    pos++;
  }

  //std::cout << std::endl;
  EdgeIterator minEdgeIt;
  ssize_t minN, minK;
  ssize_t numMin, minScore;
  minScore = 0;
  numMin   = 0;
  ssize_t success = 0;
  ssize_t trim;
  //MultTuple minTuple;
	CountedIntegralTuple minTuple;
  ssize_t minTrim;

	// This is a valid fix if the fix ended on a good note (solidExtensionFound)
  if (solidExtensionFound == 1 or 
      (solidExtensionFound == 0 and 
			 cube.size() + IntegralTuple::tupleSize - 1 > seq.length - params.maxTrim)) {
    // solidExtensionFound == 1 means that this ended 
		//                          with an edge that 'explains' the last nucleotide in 'seq'.
    // solidExtensionFound == 0 means that the method had to bail out since no 
		//                          possible explanations were available for seq at 
		//                          (pos + tupleSize - 1).  However it is possible that some amount
		//                          of the sequence was fixed, so we want to keep that, and discard 
    //                          the remainder of the read.  If the last position fixed, cube.size(),
    //                          is within params.maxTrim of the end then we can just use that.

		// No valid extensions were possible on the last spot.  We cannot trace back from here
		// so get rid of this.
    if (solidExtensionFound == 0) {
      DeleteGrid(params.maxGap, cube[cube.size()-1]);
      cube.pop_back();
    }
    ssize_t numMin;
    if ( (numMin = FindMinimumScoreEdge(cube, params, minEdgeIt, 
																				minN, minK, numMin, minScore, minTuple, 
																				minTrim, trim)) == 1 ) {
      // make sure we finished by creating an edge, otherwise 
      // no good edges were found.
      
      Backtrack(cube, cube.size()-1-trim-minK, 
								minK, minN, minTuple, fixedSeq, stats, params, firstEdit, lastSolidifiedPos );
      success       = 1;
      replaceLength = cube.size() - 1 - trim - minK;

    }
    else {
      stats.numMultiplePaths++;
    }
  }
  else {
    if (solidExtensionFound == 0) {
      stats.numNoPathFound++;
    }
		else {
			stats.numErrorAtEnd++;
		}
  }
  if (success == 1) {
    // std::cout << "fixed: " << seq.namestr << std::endl;
  }
  else {
    // std::cout << "no fix for " << seq.namestr << std::endl;
  }
  for (i = 0; i < cube.size(); i++ ) {
    DeleteGrid(params.maxGap, cube[i]);
  }
  return success;
}

ssize_t FindMinimumScoreEdge(Cube &matrix, FixParams &params, 
												 EdgeIterator &minEdge, 
												 ssize_t &minN, ssize_t &minK, ssize_t& numMin, ssize_t &minScore,
												 CountedIntegralTuple &minTuple, ssize_t &minTrim,
												 ssize_t &trim) {
  numMin = 0;
  minScore = 999999999;
  ssize_t last = matrix.size() - 1;
  // look over all nucleotides
  ssize_t n, k;
  EdgeIterator edgeIt, edgeEnd;
  minN = -1; minK = -1;
  numMin = 0;
  trim = 0;
  while (numMin != 1 and 
				 trim < params.maxTrim) {
    for (k = 0; k < params.maxGap; k++) {
      for (n = 0; n < 4; n++ ) {
				if (last - trim - k >= 0) {
					for (edgeIt = matrix[last-trim-k][k][n].begin();
							 edgeIt != matrix[last-trim-k][k][n].end();
							 ++edgeIt) {
						if (edgeIt->second.score < minScore) {
							minScore = edgeIt->second.score;
							minN     = n;
							minK     = k;
							minEdge  = edgeIt;
							numMin   = 1;
							minTuple = edgeIt->first;
							minTrim  = trim;
						}
						else if (edgeIt->second.score == minScore) {
							++numMin;
						}
					}
				}
      }
    }
    if (numMin != 1)
      ++trim;
  }
  return numMin;
}

void Backtrack(Cube &matrix, ssize_t pos, ssize_t level, ssize_t nuc, CountedIntegralTuple &tuple, 
							 DNASequence &seq, Stats &stats, FixParams &params, 
							 ssize_t &lastEdit, ssize_t &lastSolidPos) {

  // Trace a path in 'matrix' starting at the tuple referenced
  // at the cell matrix[pos][level][nuc].

  if (pos < 0) 
    return; // 0;
  std::vector<unsigned char> newSeq;
  assert(matrix[pos][level][nuc].find(tuple) != matrix[pos][level][nuc].end());
  SEdge edge;
  ssize_t length = 0;
  ssize_t trim = 0;
	//UNUSED// ssize_t firstOk= -1;
	lastEdit = -1;
	lastSolidPos = pos;
  while (matrix[pos][level][nuc][tuple].prevNuc != -1) {
    length++;
    newSeq.push_back(nuc_char[nuc]);
    edge = matrix[pos][level][nuc][tuple];

    // there must be a link to the previous cell
    assert(matrix[edge.prevPosition][edge.prevLevel][(unsigned char) edge.prevNuc].find(edge.tuple) !=
					 matrix[edge.prevPosition][edge.prevLevel][(unsigned char) edge.prevNuc].end());

    // store some statistics about the fix
    ssize_t score = matrix[pos][level][nuc][tuple].score;
    ssize_t prevScore = matrix[edge.prevPosition][edge.prevLevel][(unsigned char) edge.prevNuc][edge.tuple].score;
    
    if (score > prevScore and length < params.maxTrim) {
      trim = length;
    }
		if (score > prevScore and lastEdit == -1 ){
			lastEdit = pos;
		}
    
    if (edge.prevPosition == pos - 1 and
				prevScore < score)
      ++stats.numMut;
    else if (edge.prevLevel < level) {
//			cout << "traced an ins " <<pos << endl;
      ++stats.numIns;
    }
    else if ( edge.prevPosition < pos - 1 ) {
//			cout << "traced a del " << pos << endl;
      ++stats.numDel;
		}
    
    pos   = edge.prevPosition;;
    level = edge.prevLevel;
    nuc   = edge.prevNuc;
    tuple = edge.tuple;

  }
  if (trim > 0) {
    //    std::cout << "trimmed " << trim << std::endl;
  }

  seq.Reset(newSeq.size() - trim);
  ssize_t i;
  for (i = newSeq.size()-1 ; i >= trim; i--) {
    seq.seq[newSeq.size() - i - 1] = newSeq[i];
  }
  
}

ssize_t StoreEdge(ssize_t prevPos, ssize_t prevLevel, ssize_t prevNuc, CountedIntegralTuple prevTuple, 
							Cell& cell, CountedIntegralTuple tuple, ssize_t score, ssize_t prevScore, ssize_t prevSolid) {

  // A cell is part of the grid [pos,ins], where pos is the position in the read, and
	// ins is the number of insertions.
	
  // I do this so that I can have a fixed number of cells at each iteration of error correction.
  // There may be many different tuples that end in the same iteration
  // First try and locate the tuple in the cell.  If it is not 
  // in the cell, append it.
  // If it is in the cell, if the score to reach this tuple is better than
  // the previous score, use the current path.  Otherwise, use the previous path.

  Cell::iterator cellIt;
  cellIt = cell.find(tuple);
  if (cellIt == cell.end() ) {
    cell[tuple].prevNuc      = prevNuc;
    cell[tuple].prevLevel    = prevLevel;
    cell[tuple].prevPosition = prevPos;
    cell[tuple].tuple        = prevTuple;
    cell[tuple].score        = score;
		if (score == prevScore) {
			if (prevSolid < 255)
				cell[tuple].solidStretch = prevSolid + 1;
		}
		else {
			cell[tuple].solidStretch = 0;
		}
		return 1;
  }
  else {
    if (cellIt->second.score > score) {
      cell[tuple].prevNuc      = prevNuc;
      cell[tuple].prevLevel    = prevLevel;
      cell[tuple].prevPosition = prevPos;
      cell[tuple].tuple        = prevTuple;
      cell[tuple].score        = score;
			if (score == prevScore) {
				if (prevSolid < 255)
					cell[tuple].solidStretch = prevSolid + 1;
			}
			else {
				cell[tuple].solidStretch = 0;
			}
			return 1;
    }
		else {
			// A tuple already exists with a better score than this one.
			return 0;
		}
  }
}


void ReadToTuples(DNASequence &seq, CountedIntegralTupleList &readAsTuples) {
	CountedIntegralTuple readTuple;
		
	readAsTuples.resize(seq.length - IntegralTuple::tupleSize + 1);
	ssize_t p;
	for (p = 0; p < seq.length - IntegralTuple::tupleSize + 1; p++) {
		if (readTuple.StringToTuple(&(seq.seq[p]))) {
			readAsTuples[p] = readTuple;
			readAsTuples[p].count = 1;
		}
		else {
			readAsTuples[p].count = 0;
		}
	}
}

ssize_t SolidifyRead(DNASequence &seq, 
								 //								 CountedIntegralTupleList &spectrum,
								 CountedIntegralTupleDict &spectrum,
								 IntMatrix &scoreMat,
								 FixParams &params, Stats &stats, ssize_t &readWasModified) {

  DNASequence toFix, solidSeq;
  ssize_t replacedLength;
  ssize_t lastSolidPos;
  ssize_t firstSolidPos;
  Stats readStats;
	ssize_t changeMade;
	//UNUSED// ssize_t iter = 0;
	if (seq.length < CountedIntegralTuple::tupleSize + 1)
		return 0;
	do {

		// Assume no changes were made to the sequence
		changeMade = 0;
		/*
		 * Create a list of the tuples stored in the read at every
		 * position in the read.  That way, if the read is fixed
		 * and matches a solid pos, a short part of the read is fixed.
		 */
		CountedIntegralTupleList readAsTuples;
		ReadToTuples(seq, readAsTuples);
		
		// The read length may have changed at each iteration due to fixing 
		// an indel, so resize to fit that.


		std::string namestr= seq.namestr;
		if (! FindSolidPosition(seq, spectrum, 
														params.span, firstSolidPos, lastSolidPos)) {
			// No solid position found
			stats.numNotSolid++;
			return 0;
		}

		//
		// Otherwise, firstSolidPos and lastSolidPos are set to values between 0 
		// and the read length - 1.


		// make sure there are errors to correct
		if (firstSolidPos == 0 and lastSolidPos + IntegralTuple::tupleSize == seq.length) {
			// The full sequence is ok, just return success
			return 1;
		}

		//
		// Attempt to make the read solid from 'lastSolidPos' until the end of 
		// the read.
		// 
		if (lastSolidPos < seq.length - IntegralTuple::tupleSize) {

			// Create a reference to the original seq.
			toFix.seq     = &seq.seq[lastSolidPos];
			toFix.length  = seq.length - lastSolidPos ;
			toFix._ascii  = seq._ascii;
			ssize_t lastEdit, lastSolidifiedPos;
			if (SolidifyUntilFixed(toFix, readAsTuples, lastSolidPos, spectrum, scoreMat, params, 
														 solidSeq, replacedLength, 
														 readStats, lastEdit, lastSolidifiedPos) == 0) {
				stats.Append(readStats);
				return 0;
			}
			if (lastSolidifiedPos - lastEdit < params.edgeLimit) {
				stats.numErrorAtEnd++;
				return 0;
			}
			changeMade = 1;
			readWasModified = 1;

			PatchSeq(seq, lastSolidPos + IntegralTuple::tupleSize, solidSeq, replacedLength);
			solidSeq.Reset();
			seq._ascii = 1;
			continue;
		}
		if (firstSolidPos > 0) {
			// Now attempt to solidify the beginning of the sequence.
			DNASequence reverse;
			MakeRC(seq, reverse);
			ssize_t reversePos = reverse.length - firstSolidPos - IntegralTuple::tupleSize;
			toFix.seq = &reverse.seq[reversePos];
			toFix.length = firstSolidPos + IntegralTuple::tupleSize;
			toFix._ascii = 1;
				
			ssize_t lastEdit, lastSolidifiedPos;
			if (SolidifyUntilFixed(toFix, readAsTuples, reversePos, spectrum, scoreMat, params, 
														 solidSeq, replacedLength, readStats, lastEdit, lastSolidifiedPos) == 0) {
				stats += readStats;
				return 0;
			}
			if (lastSolidifiedPos - lastEdit < params.edgeLimit) {
				stats.numErrorAtEnd++;
				return 0;
			}
			changeMade = 1;
			readWasModified = 1;

			solidSeq._ascii = 1;
			PatchSeq(reverse, reverse.length - firstSolidPos, solidSeq, replacedLength);
			solidSeq.Reset();
			reverse._ascii = 1;

			MakeRC(reverse, seq);
			continue;
		}
		else {
			// 
			// The read was solidified until the very end of the
			// read.  This means solidifying worked.
			//
			stats += readStats;
			seq._ascii = 1;
			return 1;
		}

	} while (changeMade);
	return 0;
}

void PatchSeq(DNASequence &seq, ssize_t pos, DNASequence &patch, ssize_t replaceLength) {

  DNASequence patchedSeq;
  patchedSeq.Reset(seq.length - replaceLength + patch.length);

  patchedSeq.namestr = seq.namestr;
  // copy the unchanged segments
  patchedSeq._ascii = 0;
  ssize_t i;

  // Write nonsense into the sequence to make sure I patch every position.
  for (i = 0; i < patchedSeq.length; i++ )
    patchedSeq.seq[i] = 255;

	ssize_t patchStart = pos;
	if (pos > 0) {
		memcpy(&patchedSeq.seq[0], &seq.seq[0], pos);
  }

  // copy the new segment
  for (i = 0; i < patch.length; i++) 
    patchedSeq.seq[pos + i] = patch.seq[i];

	if (replaceLength > 0) {
		memcpy(&patchedSeq.seq[patchStart], &patch.seq[0], patch.length);
	}

	//UNUSED// ssize_t patchEnd = patchStart + replaceLength;
	memcpy(&patchedSeq.seq[patchStart + patch.length],
				 &seq.seq[patchStart + replaceLength], 
				 seq.length - (replaceLength + patchStart));
	
  seq = patchedSeq;
  patchedSeq.Reset();
}



void ReverseSeq(DNASequence &seq, DNASequence &rev) {
  // not reverse complement, just reverse
  ssize_t i;
  rev.Reset(seq.length);
  for (i = 0; i < seq.length; i++ ) {
    rev.seq[rev.length - i - 1] = seq.seq[i];
  }
}
