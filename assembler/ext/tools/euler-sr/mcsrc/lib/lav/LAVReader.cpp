/***************************************************************************
 * Title:          LAVReader.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/01/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <assert.h>
#include "LAVReader.h"
#include "LAVSequence.h"
#include "LAVAlignedContig.h"
#include "utils.h"
// Description of the parsing of the lav file.
// 
// The first routine determines what type of block is going to be read
// based off of the initial character of the lav file.
// This may be:
//  d  -  description block
//  s  -  sequences
//  h  -  contig headers
//  a  -  alignments


ssize_t LAVReader::GetStanzaType(std::ifstream &lavFile, char &type) {
  // Should we continue?
  if (lavFile.eof())
    return 0;
  lavFile >> type;
  return 1;
}

ssize_t LAVReader::DiscardLine(std::ifstream &lavIn) {
  while (!lavIn.eof() && lavIn.peek() != '\n')
    lavIn.get();
  lavIn.get();
  return !lavIn.eof();
}

ssize_t LAVReader::SkipWS(std::ifstream &lavIn) {
  char c;
  while (!lavIn.eof() && ((c = lavIn.peek())  && ( c == ' ' || 
						 c == '\t' || 
						 c == '\n')))
    lavIn.get();
  return !lavIn.eof();
}

ssize_t LAVReader::ReadStanza(std::ifstream &lavIn, LAVFile &lavFile, LAVAlignedContig *&alignedContig) {
  char type = -1;
  if (lavIn.eof()) return 0; // early bail out;
  lavIn >> type;

  LAVBlock *lavBlock = NULL;
  switch(type) {
  case 'd': 
    if (!LAVReader::ReadDescription(lavIn, lavFile)) return 0;
    break;
  case 's':
    alignedContig = new LAVAlignedContig;
    if (!LAVReader::ReadSequenceNames(lavIn, *alignedContig )) return 0;
    lavFile.alignments.push_back(alignedContig);
    break;
  case 'h':
    if (!LAVReader::ReadContigNames(lavIn, *alignedContig)) return 0;
    break;
  case 'a':
    lavBlock = new LAVBlock;
    assert(alignedContig != NULL);
    if (!LAVReader::ReadAlignment(lavIn, *lavBlock)) return 0;
    lavBlock->strand = alignedContig->qryContig.strand;
    alignedContig->alignments.push_back(lavBlock);
    break;
  case 'm':
  case 'x':
    LAVReader::ReadMaskedRegion(lavIn);
    break;
  case '#':
    LAVReader::ReadComment(lavIn);
    break;
  case EOF:
    return 0;
  default:
    std::cout << "unrecognized char in lav file: " << type << std::endl;
    exit(1);
  }

  return 1;
}

ssize_t LAVReader::ReadComment(std::ifstream &lavIn) {
  LAVReader::DiscardLine(lavIn);
  return !lavIn.eof();
}

ssize_t LAVReader::ReadBlockStart(std::ifstream &lavIn) {
  while (!lavIn.eof() && lavIn.peek() != '{') {
    lavIn.get();
  }
  lavIn.get();
  return !lavIn.eof();
}

ssize_t LAVReader::ReadBlockEnd(std::ifstream &lavIn) {
  while (!lavIn.eof() && lavIn.peek() != '}') {
    lavIn.get();
  }
  lavIn.get();
  return !lavIn.eof();
}


ssize_t LAVReader::ReadDescription(std::ifstream &lavIn, LAVFile &lavFile) {
  if (!LAVReader::ReadBlockStart(lavIn)) return 0;
  if (!LAVReader::GetString(lavIn, lavFile.blastzOpts)) return 0;
  if (!LAVReader::ReadBlockEnd(lavIn)) return 0;
  return !lavIn.eof();
}

ssize_t LAVReader::GetString(std::ifstream &lavIn, std::string &str) {
  str = "";
  // advance to the first '"'
  char c = -1;
  while (!lavIn.eof() && (c = lavIn.peek()) != '"' && c != -1  && c != 0 && c != '}') {
    lavIn.get();
  }
  if (c == '}') {
    std::cout << "malformatted string " << std::endl;
  }
  if (c == -1) return -1;
  lavIn.get(); // get rid of the first '"'

  // Get the string
  while (!lavIn.eof() && lavIn.peek() != '"') {
    str += lavIn.get();
  }
  lavIn.get(); // get rid of the second '"'
  return !lavIn.eof();
}

ssize_t LAVReader::ReadSequence(std::ifstream &lavIn, LAVSequence &sequence) {
  GetString(lavIn, sequence.sequenceName);
  lavIn >> sequence.start >> sequence.end >> sequence.strand >> sequence.contig;
  return !lavIn.eof();
}
  
ssize_t LAVReader::ReadSequenceNames(std::ifstream &lavIn, LAVAlignedContig &alignedContig) {
  if (!LAVReader::ReadBlockStart(lavIn)) return 0;
  // Read two sequences
  if (!LAVReader::ReadSequence(lavIn, alignedContig.refContig)) return 0;
  if (!LAVReader::ReadSequence(lavIn, alignedContig.qryContig)) return 0;
  if (!LAVReader::ReadBlockEnd(lavIn)) return 0;
  return !lavIn.eof();
}

ssize_t LAVReader::ReadContigNames(std::ifstream &lavIn, LAVAlignedContig &alignedContig) {
  if (!LAVReader::ReadBlockStart(lavIn)) return 0;
  if (!LAVReader::GetString(lavIn, alignedContig.refContigName)) return 0 ;
  if (!LAVReader::GetString(lavIn, alignedContig.qryContigName)) return 0;
  if (!LAVReader::ReadBlockEnd(lavIn)) return 0;
  return !lavIn.eof();
}

ssize_t LAVReader::BlockEnded(std::ifstream &lavIn) {
  if (lavIn.eof()) return 1;
  return lavIn.peek() == '}';
}

ssize_t LAVReader::ReadAlignment(std::ifstream &lavIn,
			      LAVBlock &lavBlock) {
  char lineType;
  ssize_t refAlBegin, refAlEnd, qryAlBegin, qryAlEnd, alIdent;
  if (!LAVReader::ReadBlockStart(lavIn)) return 0;
  while (!LAVReader::BlockEnded(lavIn)) {
    lavIn >> lineType;
    switch(lineType){ 
    case 's':
      lavIn >> lavBlock.score;
      break;
    case 'b':
      lavIn >> lavBlock.refBegin >> lavBlock.qryBegin;
      break;
    case 'e':
      lavIn >> lavBlock.refEnd >> lavBlock.qryEnd;
      break;
    case 'l':
      lavIn >> refAlBegin >> qryAlBegin >> refAlEnd >> qryAlEnd >> alIdent;
      lavBlock.refALBegin.push_back(refAlBegin);
      lavBlock.refALEnd.push_back(refAlEnd);
      lavBlock.qryALBegin.push_back(qryAlBegin);
      lavBlock.qryALEnd.push_back(qryAlEnd);
      lavBlock.alIdentity.push_back(alIdent);
      LAVReader::SkipWS(lavIn);
      break;
    }
  }
  if (!LAVReader::ReadBlockEnd(lavIn)) return 0;
  return !lavIn.eof();
}


void LAVReader::ReadLAVFile(std::string lavFileName, LAVFile &lavFile) {
  std::ifstream lavIn;

  openck(lavFileName, lavIn, std::ios::in, std::cout, 0);
  // read all stanzas
  LAVAlignedContig *alignedContig;
  alignedContig = NULL;
  while (LAVReader::ReadStanza(lavIn, lavFile, alignedContig)) {

  }
  lavIn.close();
}


  
ssize_t LAVReader::ReadMaskedRegion(std::ifstream &lavIn) {
  if (!LAVReader::ReadBlockStart(lavIn)) return 0;
  if (!LAVReader::ReadBlockEnd(lavIn)) return 0;
  return !lavIn.eof();
}
