/***************************************************************************
 * Title:          CommunicationLib.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _COMMUNICATION_LIB_H
#define _COMMUNICATION_LIB_H

#include "DNASequence.h"

  // Enumerate jobs for a task;
#define	doExit 1
#define	doIdle 2
#define	doUniqueMask 3
#define	doCommonMask 4
#define doEnumerate  5

class CommunicationLib {
public:
  // message ids, for keeping track of what
  // type of message is being sent/received
  static ssize_t msgCommonMaskCompute;
  static ssize_t msgCommonMaskResult;
  static ssize_t msgUniqueMaskCompute;
  static ssize_t msgUniqueMaskResult;
  static ssize_t msgEnumerateCompute;
  static ssize_t msgEnumerateResult;
  static ssize_t msgTaskType;
  static ssize_t msgUniqueMaskParams;
  static ssize_t msgCommonMaskParams;
  static ssize_t msgReport;
  static ssize_t msgAny;
  // task ids
  static ssize_t idAny;

  static ssize_t SendSequences(ssize_t myId,
		    ssize_t index1, ssize_t index2,
		    DNASequence &seq1,
		    DNASequence &seq2,
		    ssize_t destId,
		    ssize_t messageTag);
  
  static ssize_t ReceiveSequences(ssize_t sourceId,
		       ssize_t messageId,
		       ssize_t *index1, ssize_t *index2,
		       DNASequence &seq1,
		       DNASequence &seq2) ;

  static void SendMaskingInformation(ssize_t childId,
			     ssize_t dist,
			     ssize_t doIndel,
			     ssize_t hashLength,
			     ssize_t wordLength);

  static void SendCommonInformation(ssize_t childId,
				    ssize_t hashLength,
				    ssize_t wordLength);

  static void ReceiveCommonInformation(ssize_t parentId,
				       ssize_t *hashLength,
				       ssize_t *wordLength);

  static void ReceiveMaskingInformation(ssize_t parentId,
					ssize_t *dist,
					ssize_t *doIndel,
					ssize_t *hashLength,
					ssize_t *wordLength);

  static ssize_t ReceiveEnumeration(ssize_t childId,
				ssize_t *&enumerations,
				ssize_t *&locations,
				ssize_t *index0, ssize_t *index1);

  static void SendEnumeration(ssize_t parentId,
			      ssize_t myTid,
			      ssize_t length,
			      ssize_t *enumeration,
			      ssize_t *indices,
			      ssize_t index0,
			      ssize_t index1);


  static ssize_t SendPerformanceReport(ssize_t numLookups,
				    ssize_t numReferences,
				    ssize_t destId);

  static void ReceivePerformanceReport(ssize_t sourceId,
				       ssize_t *numLookups,
				       ssize_t *numReferences);


  static void InitiateTask(ssize_t childId,
			   ssize_t task);
  static void ReceiveTask(ssize_t parentId, ssize_t *task);
};


#endif
