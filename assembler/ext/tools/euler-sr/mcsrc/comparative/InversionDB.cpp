/***************************************************************************
 * Title:          InversionDB.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "InversionDB.h"
#include "blocks/dbOrthoPos.h"
#include "InversionUtils.h"

void BinInversionList(DBQuery &query,
		      InversionList &invList, 
		      std::string seqName,
		      std::vector<ssize_t> &startPos,
		      std::vector<ssize_t> &endPos,
		      BinMap &binnedInversions) {

  std::string tableName, tempTableName, chainTableName, netTableName;
  ssize_t refLen;
  std::string seqTableName = "sequences";
  InversionList::iterator invIt;
  ssize_t humanStartPos, humanEndPos;
  ssize_t humanStartStrand, humanEndStrand;
  std::set<std::string> tempTables;
  ssize_t startDeleted, endDeleted;
  for (invIt = invList.begin(); invIt != invList.end(); ++invIt) {
    // Map the inversion onto the human section.
    tableName = "human_" + (*invIt)->species + "_" + seqName;
    tempTableName = tableName + "_temp";
    chainTableName = tableName + "_chain";
    netTableName = tableName + "_net";
      
    if (!GetSequenceLength(query, seqTableName, (*invIt)->species, seqName, refLen)) {
      std::cout  << "could not get sequence length for " << (*invIt)->species << std::endl;
      exit(1);
    }
    ssize_t startSuccess, endSuccess;
    startSuccess = LookupOrthologousPosition(query, 
					     tableName,
					     0, (*invIt)->startPos,
					     refLen,
					     humanStartPos,
					     humanStartStrand,
					     tempTables, 100, startDeleted, 0, 0);
    // This may take a couple of tries.  
    endSuccess   = LookupOrthologousPosition(query, 
					     tableName, 
					     0,
					     (*invIt)->endPos,
					     refLen,
					     humanEndPos,
					     humanEndStrand,
					     tempTables, 100, endDeleted, 0,0);			      


    if (startSuccess and endSuccess and 
	humanStartPos == humanEndPos and (*invIt)->startPos != (*invIt)->endPos) {
      // Sometimes there are overlapping boundary conditions when 
      // looking up coordinates that are inverted.  Look just past the inversion
      // then map back
      std::cout << "start and end mapped to the same site, looking out of hte boundaries " << std::endl;
      startSuccess = LookupOrthologousPosition(query, 
					       tableName,
					       0, (*invIt)->startPos - 10,
					       refLen,
					       humanStartPos,
					       humanStartStrand,
					       tempTables, 100, startDeleted, 0, 0);
      humanStartPos += 10;
      endSuccess   = LookupOrthologousPosition(query, 
					       tableName, 
					       0, 
					       (*invIt)->endPos + 10,
					       refLen,
					       humanEndPos,
					       humanEndStrand,
					       tempTables, 100, endDeleted, 0, 0);			      
      humanEndPos -= 10;
    }
    
    if (startSuccess and endSuccess) { 
      ssize_t s;
      // Locate the position on human
      ssize_t startWin, endWin;
      ssize_t found = 0;
      double lenNew, lenOld;
      if (humanStartPos > humanEndPos) {
	ssize_t th = humanStartPos;
	humanStartPos = humanEndPos;
	humanEndPos = th;
      }
      std::cout << "posfor " << tableName << " " << (*invIt)->startPos << " " << (*invIt)->endPos 
		<< " " << humanStartPos << " " << humanEndPos << std::endl;
      for (s = 0; s < startPos.size() and !found; found || s++ ) {
	if ((humanStartPos >= startPos[s] and humanStartPos <= endPos[s]) or
	    (humanEndPos   >= startPos[s] and humanEndPos   <= endPos[s]) or
	    (humanStartPos <  startPos[s] and humanEndPos    >  endPos[s])) {
	  lenNew = humanEndPos - humanStartPos;
	  lenOld = endPos[s] - startPos[s];
	  std::cout << "matched " << humanStartPos << " " 
		    << humanEndPos << " " << startPos[s] << " " << endPos[s] << std::endl;
	  if ( (lenNew >= lenOld and lenNew / lenOld < 5.5) or 
		(lenOld >= lenNew and lenOld / lenNew < 5.5) ) {
	    // Found a match.
	    std::cout << "adding " << (*invIt)->species << "  " 
		      << humanStartPos << " (" << (*invIt)->startPos << ") "
		      << humanEndPos << " (" << (*invIt)->endPos << ") to " 
		      << startPos[s] << " " << endPos[s] << std::endl;
	    
	    if (humanStartPos < startPos[s])
	      startPos[s] = humanStartPos;
	    
	    if (humanEndPos > endPos[s])
	      endPos[s] = humanEndPos;
	    
	    found = 1;
	  }
	  else {
	    std::cout << "lengths didn't match " << lenNew << " " << lenOld 
		      << " " << lenNew / lenOld << " " << lenOld / lenNew << " "
		      << humanStartPos << " " << humanEndPos << " to " << startPos[s] << " " << endPos[s] 
		      << std::endl;
	  }
	}
      }
      if (!found) {
	std::cout << "starting " << (*invIt)->species << "  " 
		  << humanStartPos << " (" << (*invIt)->startPos << ") "
		  << humanEndPos << " (" << (*invIt)->endPos << ") " << std::endl;
	
	startPos.push_back(humanStartPos);
	endPos.push_back(humanEndPos);
      }
      binnedInversions[s].push_back(*invIt);
    }
  }
}

void VerifyInversion(DBQuery &query,
		     std::string refSpecies,
		     std::string qrySpecies,
		     std::string sequence,
		     StringSet &temporaryTables, ssize_t maxTempTables,
		     ssize_t startPos, ssize_t endPos,
		     ssize_t doLocal,
		     ssize_t &type,
		     std::map<std::string, DNASequence*> &seqMap,
		     FloatMatrix &scoreMat)
{
			     
  std::string qrySeqName, refSeqName;
  std::string lavTableName, chainTableName, netTableName, tempTableName;
  
  qrySeqName = qrySpecies + "." + sequence + ".fa";
  refSeqName = refSpecies + "." + sequence + ".fa";
  lavTableName = refSpecies + "_" + qrySpecies + "_" + sequence;
  chainTableName = lavTableName + "_chain";
  netTableName   = lavTableName + "_net";
  tempTableName  = lavTableName + "_temp";
  ssize_t verbosity= 0;
  //  std::cout << "verifyinversion creating temporary table: " << tempTableName << std::endl; ;
  // Store the lav alignments that are on the main net.
  CreateTemporaryTable(query,lavTableName, chainTableName, netTableName, tempTableName, 0,
		       temporaryTables, maxTempTables, verbosity);
  //  std::cout << " done " << std::endl;

  // 
  // Find the alignments that overlap with the inversion.
  // 
  query << "SELECT * from " << tempTableName << " where " 
	<< " (   tStart <= " << startPos 
	<< " AND tEnd   >= " << startPos << ") OR "
	<< " (   tStart <= " << endPos 
	<< " AND tEnd   >= " << endPos << ") OR "
	<< " (   tStart > " << startPos 
	<< " AND tEnd  < " << endPos << ")";

  query.execute();
  DNASequence *qrySeq, *refSeq;  
  //  std::cout << "getting " << refSeqName << std::endl;;
  if (seqMap.find(refSeqName) == seqMap.end()) {
    refSeq = new DNASequence;
    SeqReader::GetSeq(refSeqName, *refSeq, SeqReader::noConvert);
    //    std::cout << "got seq " << refSeqName << " of len: " << refSeq->length << std::endl;
    seqMap[refSeqName] = refSeq;
  }
  else {
    refSeq = seqMap[refSeqName];
  }
  // Read in the query sequence if need be 
  if (seqMap.find(qrySeqName) == seqMap.end()) {
    qrySeq = new DNASequence;
    SeqReader::GetSeq(qrySeqName, *qrySeq, SeqReader::noConvert);
    seqMap[qrySeqName] = qrySeq;
  }
  else {
    qrySeq = seqMap[qrySeqName];
  }
  MYSQL_RES *res;
  MYSQL_ROW row;
  res = mysql_store_result(query.conn);
  ssize_t numRows = 0;
  type = -1;
  if (res != NULL) 
    numRows = mysql_num_rows(res);
  ssize_t alignmentFound = 0;
  if (res != 0 and numRows > 0) {
    // Found an alignment overlapping with the inversion.  This
    // is probably not sufficient, and will need to have a
    // somewhat high score.

    // Continue to the next alignment.
    LAVRow lavRow;
    std::vector<LAVRow> alignment; 
    ssize_t overlapLength = 0;
    while (row = mysql_fetch_row(res)) {
      GetBlock(row, lavRow);
      alignment.push_back(lavRow);
    }
    double score;
    ssize_t lengthAligned;
    double percentIdentity;
    // Found an alignment in the forward strand, score it to make sure it makes sense.
    ssize_t numMatches, numMismatches;
    score = ScoreRows(*refSeq, *qrySeq, 
		      alignment, startPos, endPos, scoreMat, 
		      lengthAligned,  
		      numMatches, numMismatches, percentIdentity );
    
    /*
      std::cout << "score of overlap: " << score << " " 
      << lengthAligned << " " << percentIdentity 
      << " " << lengthAligned / (1.0 * endPos - startPos) << std::endl;
    */
    if ((endPos- startPos > 0 ) and 
				lengthAligned > (endPos - startPos)/2 ) {
      // The alignment covered at least 50% of the inversion, that's good enough.
      // Chances are the alignment just covers a few percent of the sequence.
      alignmentFound = 1;
      type = WITHOUT_INVERSION;
    }
  }
  if (!alignmentFound) {
    // No alignments were found in the area of the inversion.
    // Try to align the inversion to the unaligned sequence.
    ssize_t qStart, qEnd, tStart, tEnd;

    // Didn't locate any alignments that overlap with the inversion here.  
    // Try to locate the inversion
    // by sw alignment of the unaligned orthologous region.

    // Step 1.  Find the boundaries of the unaligend region where the inversion should be.

    ssize_t maxTEnd, minTStart;
    // Get the block before the gap
    qStart = -1;
    qEnd   = -1;
    LAVRow prevRow, nextRow;

    FetchSurroundingBlocks(query, tempTableName, 
			   startPos, endPos, -1, -1, prevRow, nextRow, verbosity);
	  
    tStart = prevRow.tEnd; tEnd = nextRow.tStart;
    qStart = prevRow.qEnd; qEnd = nextRow.qStart;
    
    if (qStart > qEnd) {
      ssize_t tempPos = qStart;
      qStart = qEnd;
      qEnd = tempPos;
    }

    // Check for deletion.
    ssize_t deletion;
    deletion = CheckForDeletion(tStart, tEnd, qStart, qEnd, 10.0);
    if (deletion) {
      //      foundWith = 0; foundWithout = 0; foundUnknown = 0; foundGapped = 1;
      //      std::cout << "got deletion in temp pos " << std::endl;
      type = DELETED;
    }
    else {
      if (!doLocal or qStart < 0 or qEnd < 0 or (qEnd - qStart < 0)) {
	//	foundWith = 0; foundWithout = 0; foundUnknown = 1;
	type = UNKNOWN;
      }
      else {
	// Store the coordinates of the unaligned region.
	      
	// Step 2. Extract the sequences corresponding to the
	// inversion and the unaligned sequence. 
	      
	// Step 2.1.  Open the query sequence.
	//	    qrySeq.Reset();
	if (doLocal) {
	  DoLocalAlign(*refSeq, *qrySeq,
		       startPos, endPos, qStart, qEnd,
		       type);
	  //	  std::cout << "did local align and got type: " << type << std::endl;
	}
      }
    }
  }
}
