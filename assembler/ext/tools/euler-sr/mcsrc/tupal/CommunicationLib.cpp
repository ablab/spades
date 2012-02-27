/***************************************************************************
 * Title:          CommunicationLib.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "CommunicationLib.h"
#include "pvm3.h"

ssize_t CommunicationLib::msgUniqueMaskCompute = 1;
ssize_t CommunicationLib::msgUniqueMaskResult  = 2;
ssize_t CommunicationLib::msgTaskType          = 3;
ssize_t CommunicationLib::msgUniqueMaskParams  = 4;
ssize_t CommunicationLib::msgCommonMaskCompute = 5;
ssize_t CommunicationLib::msgCommonMaskResult  = 6;
ssize_t CommunicationLib::msgEnumerateCompute = 7;
ssize_t CommunicationLib::msgEnumerateResult  = 8;
ssize_t CommunicationLib::msgReport            = 9;
ssize_t CommunicationLib::msgCommonMaskParams  = 4;
ssize_t CommunicationLib::msgAny = -1;
ssize_t CommunicationLib::idAny = -1;


ssize_t CommunicationLib::SendSequences(ssize_t myId,
				    ssize_t index1, ssize_t index2,
				    DNASequence &seq1,
				    DNASequence &seq2,
				    ssize_t destId,
				    ssize_t messageTag) {

  ssize_t info;
  ssize_t length1, length2;
  info = pvm_initsend( PvmDataDefault );
  // Pack who this is from since is possible the recipient
  // is listening for anything.
  info = pvm_pkint(&myId, 1,1);

  // pack the indices. These may be used to describe the 
  // sequences.
  info = pvm_pkint(&index1, 1,1);
  assert(info >= 0);
  info = pvm_pkint(&index2, 1,1);
  assert(info >= 0);

  // pack the lengths
  length1 = seq1.length;
  length2 = seq2.length;
  info = pvm_pkint(&length1, 1,1);
  assert(info >= 0);
  info = pvm_pkint(&length2, 1,1);
  assert(info >= 0);

  // pack the sequence
  info = pvm_pkbyte(seq1.seq, seq1.length, 1);
  assert(info >= 0);
  info = pvm_pkbyte(seq2.seq, seq2.length, 1);
  assert(info >= 0);


  // now send information about each sequence
  pvm_pkint(&seq1.startEnumeration, 1,1);
  pvm_pkint(&seq2.startEnumeration, 1,1);

  pvm_pkint(&seq1.startPosition, 1,1);
  pvm_pkint(&seq2.startPosition, 1,1);
  
  info = pvm_send(destId, messageTag);
  assert(info >= 0);

  return 1; // Handle errors better later on.
}


ssize_t CommunicationLib::ReceiveSequences(ssize_t sourceId,
				       ssize_t messageId,
				       ssize_t *index1, ssize_t *index2,
				       DNASequence &seq1,
				       DNASequence &seq2) 
{
  ssize_t bufid, info;
  // Receive the result of some computation
  ssize_t length1, length2;
  ssize_t childId;
  bufid = pvm_recv(sourceId, messageId);
  
  // Receive the index of the child (to give more work)
  info = pvm_upkint(&childId, 1, 1);
  // Receive what fragments were masked
  info = pvm_upkint(index1, 1, 1);
  info = pvm_upkint(index2, 1, 1);
  // Receive the lengths of the fragments
  info = pvm_upkint(&length1, 1,1);
  info = pvm_upkint(&length2, 1,1);
  assert(length1 >= 0);
  assert(length2 >= 0);
  // Receive the masked fragments
  // allocate new sequences
  seq1.Reset(length1);
  seq2.Reset(length2);
  info = pvm_upkbyte(seq1.seq, seq1.length, 1);
  info = pvm_upkbyte(seq2.seq, seq2.length, 1);


  pvm_upkint(&seq1.startEnumeration, 1,1);
  pvm_upkint(&seq2.startEnumeration, 1,1);

  pvm_upkint(&seq1.startPosition, 1,1);
  pvm_upkint(&seq2.startPosition, 1,1);

  return childId;
}

void CommunicationLib::SendMaskingInformation(ssize_t childId,
					      ssize_t dist,
					      ssize_t doIndel,
					      ssize_t hashLength,
					      ssize_t wordLength) {
  pvm_initsend( PvmDataDefault ) ;
  pvm_pkint(&dist, 1,1);
  pvm_pkint(&doIndel, 1,1);
  pvm_pkint(&hashLength, 1,1);
  pvm_pkint(&wordLength, 1,1);
  pvm_send(childId, msgUniqueMaskParams);
}



void CommunicationLib::ReceiveMaskingInformation(ssize_t parentId,
						 ssize_t *dist,
						 ssize_t *doIndel,
						 ssize_t *hashLength,
						 ssize_t *wordLength) {
  pvm_recv(parentId, CommunicationLib::msgUniqueMaskParams);
  pvm_upkint(dist, 1,1);
  pvm_upkint(doIndel, 1,1);
  pvm_upkint(hashLength, 1,1);
  pvm_upkint(wordLength, 1,1);
}


void CommunicationLib::SendCommonInformation(ssize_t childId,
					     ssize_t hashLength,
					     ssize_t wordLength) {
  pvm_initsend( PvmDataDefault ) ;
  pvm_pkint(&hashLength, 1,1);
  pvm_pkint(&wordLength, 1,1);
  pvm_send(childId, msgCommonMaskParams);
}

void CommunicationLib::ReceiveCommonInformation(ssize_t parentId,
						 ssize_t *hashLength,
						 ssize_t *wordLength) {
  pvm_recv(parentId, msgCommonMaskParams);
  pvm_upkint(hashLength, 1,1);
  pvm_upkint(wordLength, 1,1);
}


void CommunicationLib::SendEnumeration(ssize_t parentId,
				       ssize_t myId,
				       ssize_t length,
				       ssize_t *enumeration,
				       ssize_t *indices,
				       ssize_t index0,
				       ssize_t index1) {
  std::cout << "sending enumeration for indices " 
	    << index0 << " " << index1 << std::endl;
    
  pvm_initsend( PvmDataDefault );
  pvm_pkint(&myId, 1, 1);
  pvm_pkint(&index0, 1, 1);
  pvm_pkint(&index1, 1, 1);
  pvm_pkint(&length, 1, 1);
  pvm_pkint(enumeration, length, 1);
  pvm_pkint(indices, length, 1);
  pvm_send(parentId, msgEnumerateResult);
}


ssize_t CommunicationLib::ReceiveEnumeration(ssize_t childId,
					  ssize_t *&enumerations,
					  ssize_t *&locations,
					  ssize_t *index0, ssize_t *index1) {
  ssize_t length;
  pvm_recv(childId, msgEnumerateResult);
  pvm_upkint(&childId, 1, 1);
  pvm_upkint(index0, 1, 1);
  pvm_upkint(index1, 1, 1);
  pvm_upkint(&length, 1, 1);
  enumerations = new ssize_t[length];
  locations    = new ssize_t[length];
  pvm_upkint(enumerations, length, 1);
  pvm_upkint(locations, length, 1);
  return childId;
}
					  
void CommunicationLib::ReceiveTask(ssize_t parentId, ssize_t *task) {
  pvm_recv(parentId, msgTaskType);
  pvm_upkint(task, 1,1);
}

void CommunicationLib::InitiateTask(ssize_t childId,
				    ssize_t task) {
  pvm_initsend( PvmDataDefault );
  pvm_pkint(&task, 1,1);
  pvm_send(childId, msgTaskType);
}

ssize_t CommunicationLib::SendPerformanceReport(ssize_t numLookups,
					    ssize_t numReferences,
					    ssize_t destId) {
  pvm_initsend( PvmDataDefault );
  pvm_pkint(&numLookups, 1, 1);
  pvm_pkint(&numReferences, 1, 1);
  ssize_t info;
  info = pvm_send(destId, CommunicationLib::msgReport);
  return info;
}

void CommunicationLib::ReceivePerformanceReport(ssize_t sourceId,
						ssize_t *numLookups,
						ssize_t *numReferences) {

  pvm_recv( sourceId, CommunicationLib::msgReport );
  pvm_upkint( numLookups, 1, 1 );
  pvm_upkint( numReferences, 1 ,1 );

}
