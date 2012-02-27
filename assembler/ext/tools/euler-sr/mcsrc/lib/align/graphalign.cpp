/***************************************************************************
 * Title:          graphalign.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "graphalign.h"
#include "alignutils.h"
#include <limits.h>
#include <queue>

_INT_ AlignVertex::none	 = -1;
_INT_ AlignVertex::match = 0;
_INT_ AlignVertex::refGap= 1;
_INT_ AlignVertex::qryGap= 2;

_INT_ AlignVertex::infty = INT_MAX;

char nuc[4] = {'a', 'c', 'g', 't'};

std::ostream & operator<<(std::ostream & out, const AlignVertex &v) {
	out << "ms: " << v.matchScore << " rs: " << v.refGapScore << " qs: " << v.qryGapScore;
	return out;
}	 

// Function AlignRegion. 
// Given the starting and ending positions in a reference and query sequence,
// extend alignments in each direction.
// Note that the starting and ending positions may be overlapping.
ssize_t AlignRegion(DNASequence &refSeq, ssize_t refStartPos, ssize_t refEndPos,
									DNASequence &qrySeq, ssize_t qryStartPos, ssize_t qryEndPos,
									IntMatrix &scoreMat, ssize_t  gapOpen, ssize_t gapExtend,
									ssize_t &refAlignStart, ssize_t &qryAlignStart, 
									ssize_t *&alignment, ssize_t &length) {

	//UNUSED// ssize_t forwardAlign;
	//UNUSED// ssize_t reverseAlign;
	
	// The initial score is the score of the overlapping regions from start to end.

	ssize_t initialScore = 0;

	ssize_t refPos, qryPos;
	ssize_t refIndex, qryIndex;
	assert((refEndPos - refStartPos) == (qryEndPos - qryStartPos));

	for (refPos = refStartPos, qryPos = qryStartPos; refPos <= refEndPos; ++refPos, ++qryPos) {
		refIndex = nuc_index[(unsigned char) refSeq.seq[refPos]];
		qryIndex = nuc_index[(unsigned char) qrySeq.seq[qryPos]];
		assert(refIndex >= 0 and refIndex <= 4 and qryIndex >= 0 and qryIndex <= 4);
		initialScore += scoreMat[refIndex][qryIndex];
	}

	ssize_t *forAlignment, forAlignLength;
	ssize_t *revAlignment, revAlignLength;

	ssize_t forRefStart, forQryStart, revRefStart, revQryStart;
	ssize_t revScore, forScore;
	forScore = ScoreBandedAffineAlign(refSeq, qrySeq, initialScore,
																		scoreMat, gapOpen, gapExtend,
																		initialScore + 2000,
																		forAlignment, forAlignLength,
																		forRefStart, forQryStart,
																		refEndPos + 1, qryEndPos + 1,
																		FORWARD_DIR);

	revScore = ScoreBandedAffineAlign(refSeq, qrySeq, initialScore,
																		scoreMat, gapOpen, gapExtend,
																		initialScore + 2000,
																		revAlignment, revAlignLength,
																		revRefStart, revQryStart,
																		refStartPos -1, qryStartPos -1,
																		REVERSE_DIR);
	refAlignStart = revRefStart;
	qryAlignStart = revQryStart;
	// Now join together the result.
	ssize_t overlapLength = refEndPos - refStartPos + 1;
	length = revAlignLength + overlapLength + forAlignLength;

	alignment = new ssize_t[length];
	
	/*
		std::cout << refStartPos + 1
		<< " revScore: " << revScore << " length: " << revAlignLength << std::endl;
		std::cout << refEndPos + 1 
		<< " forScore: " << forScore << " length: " << forAlignLength << std::endl;
	*/
	ssize_t pos, alignPos;
	alignPos = 0;
	for (pos = 0; pos < revAlignLength; pos++, alignPos++) {
		if (revAlignment[alignPos] != -1) 
			alignment[alignPos] = revAlignment[pos] + revQryStart;
		else 
			alignment[alignPos] = -1;
	}
	
	for (pos = refStartPos, qryPos = qryStartPos; 
			 pos <= refEndPos; 
			 pos++, qryPos++, alignPos++ ) {
		alignment[alignPos] = qryPos;
	}

	for (pos = 0; pos < forAlignLength; pos++, alignPos++) {
		if (forAlignment[pos] != -1)
			alignment[alignPos] = forQryStart + forAlignment[pos];
		else 
			alignment[alignPos] = -1;
	}

	delete[] revAlignment;
	delete[] forAlignment;

	
	return revScore + forScore + initialScore;
}




ssize_t ScoreBandedAffineAlign(// sequences to compare 
														 DNASequence &refSeq, DNASequence &qrySeq, 
														 // score of seed
														 ssize_t initialScore,	
														 // scoring parameters
														 IntMatrix &scoreMat, ssize_t gapOpen, ssize_t gapExtend, 
														 // negative score is good.
														 ssize_t maxScore,	
														 ssize_t *&locations, ssize_t &length, 
														 ssize_t &refAlignStart, ssize_t &qryAlignStart, 
	 // used for storing the resulting map
	 // refStartPos is the starting position of the alignment.	For alignments 
	 // done in the forward position, this is simply the starting position of the 
	 // alignment.	For alignments done in the reverse direction, this is the ending 
	 // position of the alignment.
														 ssize_t refStartPos, ssize_t qryStartPos, 
														 ssize_t extDir
														 ) { 
	ssize_t refPos, qryPos;
	// reference is on vertical axis, query on horizontal
	qryPos = 0; refPos = 0;
	//UNUSED+// AlignVertex *prevVertex;
	AlignVertex *curVertex,  *newVertex;
	AlignVertex *root;
	AlignVertex *minScoringVertex;
	ssize_t minScore;
	ssize_t minScoringType;
	std::queue<AlignVertex*> newVertices;
	std::vector<AlignVertex*> graph;
	ssize_t refExtLength, qryExtLength;	 
	if (extDir == 1) {
		refExtLength = refSeq.length - refStartPos;
		qryExtLength = qrySeq.length - qryStartPos;
	}
	else {
		refExtLength = refStartPos - 1;
		qryExtLength = qryStartPos - 1;
	}
	if (refExtLength	== 0 || qryExtLength == 0) {
		locations = NULL;
		length = 0;
		return -1;
	}
	
	refPos = -1;
	qryPos = -1;
	root = new AlignVertex(refPos, qryPos);
	root->matchScore = initialScore;
	minScoringVertex = root;
	minScore = initialScore;
	minScoringType	 = AlignVertex::match;

	// Initialize the reference gap.
	curVertex = root;

	// Create the first aligned position

	newVertices.push(root);
	while (refPos < refSeq.length && // possibly iterate over all positions in ref seq
				 curVertex != NULL && // not sure if this is necessary, is the case when new
				 // vertices->size() == 0
				 newVertices.size() > 0	 // when there were no starting positions in the query seq, 
				 // the alignment is done.
				 ) {
		curVertex = newVertices.front();
		newVertices.pop();
		graph.push_back(curVertex);
		if (curVertex != NULL && curVertex->Capable(minScore + 1500)) {
			refPos = curVertex->refPos;
			qryPos = curVertex->qryPos;

			// Keep track of best-scoring vertex, and position in that vertex
			// that is the best scoring (match, ins/del).
			if (curVertex->matchScore < minScoringVertex->matchScore &&
					curVertex->matchScore < minScore) {
				minScoringVertex = curVertex;
				minScoringType	 = AlignVertex::match;
				minScore = curVertex->matchScore;
			}
			else if (curVertex->refGapScore < minScoringVertex->refGapScore &&
							 curVertex->refGapScore < minScore) {
				minScoringVertex = curVertex;
				minScoringType	 = AlignVertex::refGap;
				minScore				 = curVertex->refGapScore;
			}
			else if (curVertex->qryGapScore < minScoringVertex->qryGapScore &&
							 curVertex->qryGapScore < minScore) {
				minScoringVertex = curVertex;
				minScoringType	 = AlignVertex::qryGap;
				minScore				 = curVertex->qryGapScore;
			}

			if (curVertex->GetEast() == NULL && 
					qryPos < qryExtLength -1 && 
					curVertex->ValidQryGap(gapOpen, minScore + 1500)) {
				newVertex = new AlignVertex(refPos, qryPos+1);
				curVertex->InitializeEast(newVertex);
				newVertex->Score(refSeq, qrySeq, refPos, qryPos+1, 
												 refStartPos, qryStartPos, extDir,
												 scoreMat, gapOpen, gapExtend, minScore + 1500, minScore);	
				newVertices.push(newVertex);
			}
			if (curVertex->GetSouth() == NULL && 
					refPos < refExtLength -1 && 
					curVertex->ValidRefGap(refPos, minScore + 1500)) {
				newVertex = new AlignVertex(refPos+1, qryPos);
				curVertex->InitializeSouth(newVertex);
				newVertex->Score(refSeq, qrySeq, refPos+1, qryPos, 
												 refStartPos, qryStartPos, extDir,
												 scoreMat, gapOpen, gapExtend, minScore + 1500, minScore);
				newVertices.push(newVertex);
			}
			if (curVertex->GetSoutheast() == NULL && 
					refPos < refExtLength-1 && 
					qryPos < qryExtLength-1 && 
					curVertex->ValidMatch(minScore + 1500)) {
				newVertex = new AlignVertex(refPos+1, qryPos+1);
				curVertex->InitializeSouthEast(newVertex);
				newVertex->Score(refSeq, qrySeq, refPos+1, qryPos+1, refStartPos, qryStartPos, extDir,
												 scoreMat, gapOpen, gapExtend, minScore	 + 1500, minScore);
				newVertices.push(newVertex);
			}
		}
	}
	// Record where in each sequence the alignment starts.
	if (extDir == 1) {
		refAlignStart = refStartPos;
		qryAlignStart = qryStartPos;
	}	 
	else {
		refAlignStart = refStartPos - minScoringVertex->refPos;
		qryAlignStart = qryStartPos - minScoringVertex->qryPos;
	}
	GetGraphAlignment(minScoringVertex, locations, length, refStartPos, qryStartPos, extDir);
	minScore = minScoringVertex->matchScore;
	DeleteGraph(graph);
	return minScore;
}

void DeleteGraph(std::vector<AlignVertex*> &graph) {
	ssize_t i;
	for ( i= 0; i < graph.size(); i++) {
		delete graph[i];
	}
}

void GetGraphAlignment(AlignVertex *vertex, ssize_t *&locations, 
											 ssize_t &length, 
											 ssize_t refStartPos, ssize_t qryStartPos, ssize_t extDir
											 ) {
	assert(vertex != NULL);
	_INT_ alignType;
	// optimal alignment ends at vertex->refPos, and starts according to alignment specified by the parameters.
	length = vertex->refPos + 1;
	ssize_t qryLength = vertex->qryPos;
	locations = new ssize_t[length];
	ssize_t i;
	
	for (i = 0; i < length; i++) locations[i] = -1;

	// start on the minimum score for this vertex
	if (vertex->matchScore <= vertex->refGapScore &&
			vertex->matchScore <= vertex->qryGapScore)
		alignType = AlignVertex::match;
	else if (vertex->refGapScore <= vertex->qryGapScore)
		alignType = AlignVertex::refGap;
	else
		alignType = AlignVertex::qryGap;
	while (vertex->refPos != -1 && vertex->qryPos != -1) {
		if (alignType == AlignVertex::match) {
			alignType = vertex->matchIndex;
			assert(vertex->refPos < length);
			if (alignType == AlignVertex::match) {
				if (extDir == 1)
					locations[vertex->refPos] = vertex->qryPos;
				else
					locations[length - vertex->refPos - 1] = qryLength - vertex->qryPos;
			}
			vertex = vertex->matchPrev; 
		}
		else if (alignType == AlignVertex::refGap) {
			alignType = vertex->refGapIndex;
			vertex = vertex->refGapPrev;
		}
		else if (alignType == AlignVertex::qryGap) {
			alignType = vertex->qryGapIndex;
			vertex = vertex->qryGapPrev;
		}
		else {
			assert(printf("align type %d not properly set\n", alignType) == 0);
		}
	}
}


AlignVertex* AlignVertex::GetEast() {
	// try various paths to get east
	if (east == NULL)
		if (north != NULL) {
			if (north->southeast != NULL) 
				east = north->southeast;
			else if (north->east != NULL)
				east = north->east->south;
		}

	return east;
}

AlignVertex* AlignVertex::GetSouth() {
	// try	various paths to get south
	if (south == NULL) {
		if (west != NULL && west->southeast != NULL)
			south = west->southeast;
		else if (west != NULL && west->south != NULL && west->south->east != NULL)
			south = west->south->east;
	}
	return south;
}

AlignVertex* AlignVertex::GetSoutheast() {
	if (southeast == NULL) {
		if (south != NULL && south->east != NULL)
			southeast = south->east;
		if (east != NULL && east->south != NULL)
			southeast = east->south;
	}
	return southeast;
}

AlignVertex* AlignVertex::GetNorthWest() {
	if (northwest != NULL)
		return northwest;

	if (north != NULL) {
		northwest = north->west;
		return north->west;
	}

	if (west != NULL) {
		northwest = west->north;
		return west->north;
	}

	return northwest; // TODO: check
}

AlignVertex*	AlignVertex::GetNorth() {
	if (north != NULL)
		return north;

	if (northwest != NULL) {
		north = northwest->east;
		return northwest->east;
	}

	return north; // TODO: check
}

AlignVertex* AlignVertex::GetWest() {
	if (west != NULL)
		return west;

	if (northwest != NULL) {
		west = northwest->south;
		return northwest->south;
	}

	return west; // TODO: check
}

AlignVertex* AlignVertex::InitializeEast(AlignVertex *vertex) {
	assert(vertex != NULL);

	// Initialize parameter
	vertex->west = this;
	this->east	 = vertex;

	vertex->northwest = this->north;
	if (this->north != NULL) 
		vertex->north = this->north->east;
	else
		vertex->north = NULL;

	return vertex;
}


AlignVertex* AlignVertex::InitializeSouth(AlignVertex *vertex) {
	assert(vertex != NULL);

	vertex->north = this;
	this->south		= vertex;

	vertex->northwest = this->west;
	if (this->west != NULL)
		vertex->west	= this->west->south;
	else
		vertex->west = NULL;
	
	return vertex;
}

AlignVertex* AlignVertex::InitializeSouthEast(AlignVertex *vertex) {
	assert(vertex != NULL);
	
	vertex->northwest = this;
	this->southeast = vertex;

	vertex->west			= this->south;
	if (vertex->west != NULL)
		vertex->west->east = vertex;

	vertex->north			= this->east;
	if (vertex->north != NULL)
		vertex->north->south = vertex;

	return vertex;
}


ssize_t AlignVertex::ScoreRefGap(ssize_t gapOpenCost, ssize_t gapExtendCost, ssize_t maxScore, ssize_t minScore) {

	ssize_t gapOpenScore, gapExtendScore;

	GetNorth();

	if (north == NULL) {
		refGapScore = infty;
	}
	else {
		gapOpenScore = north->matchScore + gapOpenCost;
		gapExtendScore = north->refGapScore + gapExtendCost;

		if (gapOpenScore < gapExtendScore) {
			refGapScore = gapOpenScore;
			refGapPrev	= north;
			refGapIndex = match;
		}
		else {
			refGapScore = gapExtendScore;
			refGapPrev	= north;
			refGapIndex = refGap;
		}
	}
	if (refGapScore > maxScore) 
		refGapScore = infty;

	refGapScore = std::max(refGapScore, 2*minScore);
	return refGapScore;
}


ssize_t AlignVertex::ScoreQryGap(ssize_t gapOpenCost, ssize_t gapExtendCost, ssize_t maxScore, ssize_t minScore) {

	ssize_t gapOpenScore, gapExtendScore;

	GetWest();

	if (west == NULL) {
		qryGapScore = infty;
	}
	else {
		gapOpenScore = west->matchScore + gapOpenCost;
		gapExtendScore = west->qryGapScore + gapExtendCost;
		
		if (gapOpenScore < gapExtendScore) {
			qryGapScore = gapOpenScore;
			qryGapPrev	= west;
			qryGapIndex = match;
		}
		else {
			qryGapScore = gapExtendScore;
			qryGapPrev	= west;
			qryGapIndex = qryGap;
		}
	}
	if (qryGapScore > maxScore)
		qryGapScore = infty;
	
	qryGapScore = std::max(qryGapScore, 2*minScore);
	return qryGapScore;
}


ssize_t AlignVertex::ScoreMatch(DNASequence &refSeq, DNASequence &qrySeq,
															ssize_t refPos, ssize_t qryPos,
															ssize_t refStartPos, ssize_t qryStartPos, ssize_t extDir,
															IntMatrix &scoreMat, ssize_t maxScore, ssize_t minScore) {
	ssize_t diagScore;
	GetNorthWest();
	if (northwest == NULL || refPos < 0 || qryPos < 0)
		diagScore = infty;
	else {
		ssize_t rn, qn;
		rn = nuc_index[(unsigned char)refSeq.seq[refPos*extDir + refStartPos]];
		qn = nuc_index[(unsigned char)qrySeq.seq[qryPos*extDir + qryStartPos]];
		diagScore = northwest->matchScore + scoreMat[rn][qn];
		assert(refPos*extDir+refStartPos < refSeq.length && refPos*extDir+refStartPos >= 0);
		assert(qryPos*extDir+qryStartPos < qrySeq.length && qryPos*extDir+qryStartPos >= 0);
	}
	
	matchScore = std::min(diagScore, 
												std::min(refGapScore, qryGapScore));

	if (matchScore == diagScore) {
		matchPrev	 = northwest;
		matchIndex = match;
	}
	else if (matchScore == refGapScore) {
		matchPrev	 = this;
		matchIndex = refGap;
	}
	else if (matchScore == qryGapScore) {
		matchPrev	 = this;
		matchIndex = qryGap;
	}

	if (matchScore > maxScore) 
		matchScore = infty;

	matchScore = std::max(matchScore, 2*minScore);
	return matchScore;
}

ssize_t AlignVertex::Perfect(ssize_t maxScore) {
	return refGapScore < maxScore and qryGapScore < maxScore and matchScore < maxScore;
}

ssize_t AlignVertex::Capable(ssize_t maxScore) {
	return refGapScore < maxScore or qryGapScore < maxScore or matchScore < maxScore;
}


