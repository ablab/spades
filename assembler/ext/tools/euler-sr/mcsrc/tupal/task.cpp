/***************************************************************************
 * Title:          task.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include "pvm3.h"
#include "DNASequence.h"
#include "TupleLib.h"
#include "CommunicationLib.h"



int main(int argc, char** argv) {
  ssize_t myTid;                  /* my task id */
  
  ssize_t nhost, narch;
  struct pvmhostinfo *hostp;

  // Initialize PVM
  /* enroll in pvm */
  myTid = pvm_mytid();
  pvm_config( &nhost, &narch, &hostp );

  // Learn who will be giving the orders
  ssize_t parentTid = pvm_parent();
  ssize_t taskType;

  // Always work on two sequences
  DNASequence seq1, seq2;

  pvm_initsend( PvmDataDefault );
  char hostname[200];
  gethostname(hostname, 200);
  pvm_pkstr(hostname);
  pvm_send(parentTid, 1);
  // parameters for unique masking
  ssize_t dist;
  ssize_t doIndel;
  ssize_t hashLength, wordLength;
  ssize_t index1, index2;

  while (1) {
    CommunicationLib::ReceiveTask(parentTid, &taskType);
    switch (taskType) { 
    case doExit:
      // Done processing.
      pvm_exit();
      return 0;
      continue;
    case doUniqueMask:
      // Receive masking paramters.
      CommunicationLib::ReceiveMaskingInformation(parentTid,
						  &dist,
						  &doIndel,
						  &hashLength,
						  &wordLength);
      CommunicationLib::ReceiveSequences(parentTid, 
					 CommunicationLib::msgUniqueMaskCompute,
					 &index1, &index2,
					 seq1, seq2);
      MaskNotUnique(seq1, seq2, dist, doIndel, hashLength, wordLength);

      CommunicationLib::SendSequences(myTid, 
				      index1, index2,
				      seq1, seq2,
				      parentTid,
				      CommunicationLib::msgUniqueMaskResult);
      CommunicationLib::SendPerformanceReport(Node::totalLookups, 
					      Node::totalReferences, 
					      parentTid);
      continue;
    case doCommonMask:
      // Get information
      CommunicationLib::ReceiveCommonInformation(parentTid,
						 &hashLength,
						 &wordLength);

      CommunicationLib::ReceiveSequences(parentTid, 
					 CommunicationLib::msgCommonMaskCompute,
					 &index1, &index2,
					 seq1, seq2);

      // Do work
      UnmaskShared(seq1, seq2, hashLength, wordLength);

      // Send it back
      CommunicationLib::SendSequences(myTid, 
				      index1, index2,
				      seq1, seq2,
				      parentTid,
				      CommunicationLib::msgCommonMaskResult);
      continue;
    case doEnumerate:
      CommunicationLib::ReceiveCommonInformation(parentTid,
						 &hashLength,
						 &wordLength);
      
      CommunicationLib::ReceiveSequences(parentTid, 
					 CommunicationLib::msgEnumerateCompute,
					 &index1, &index2,
					 seq1, seq2);

      ssize_t i;
      ssize_t *locations, *enumerations;
      locations = new ssize_t[seq1.length];
      enumerations = new ssize_t[seq1.length];

      for (i = 0; i < seq1.length; i++) {
	locations[i] = -1;
	enumerations[i] = 0;
      }

      EnumerateUnique(seq1, seq2,
		      hashLength, wordLength, 
		      enumerations,
		      locations
		      );

      CommunicationLib::SendEnumeration(parentTid,
					myTid,
					seq1.length,
					enumerations,
					locations,
					index1, index2);

      delete locations;
      delete enumerations;
      continue;
    }
  }
  return 0;
  pvm_exit();
}
