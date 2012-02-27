/***************************************************************************
 * Title:          LAVReader.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _LAVREADER_
#define _LAVREADER_ 

#include <iostream>
#include <fstream>

#include "LAVFile.h"
#include "LAVBlock.h"
#include "LAVSequence.h"

class LAVReader {
  static ssize_t GetStanzaType(std::ifstream &lavFile, char &type);
  static ssize_t DiscardLine(std::ifstream &lavIn);
  static ssize_t SkipWS(std::ifstream &lavIn);
  static ssize_t ReadStanza(std::ifstream &lavIn, LAVFile &lavFile, LAVAlignedContig *&lavAlignedContig);
  static ssize_t ReadBlockStart(std::ifstream &lavIn);
  static ssize_t ReadBlockEnd(std::ifstream &lavIn);
  static ssize_t ReadDescription(std::ifstream &lavIn, LAVFile &lavFile) ;
  static ssize_t GetString(std::ifstream &lavIn, std::string &str);
  static ssize_t ReadSequence(std::ifstream &lavIn, LAVSequence &sequence) ;
  static ssize_t ReadSequenceNames(std::ifstream &lavIn, LAVAlignedContig &alignedContig);
  static ssize_t ReadContigNames(std::ifstream &lavIn, LAVAlignedContig &alignedContig);
  static ssize_t BlockEnded(std::ifstream &lavIn);
  static ssize_t ReadAlignment(std::ifstream &lavFile, LAVBlock &lavBlock);
  static ssize_t ReadComment(std::ifstream &lavIn);
  static ssize_t ReadMaskedRegion(std::ifstream &lavIn);
public:
  static void ReadLAVFile(std::string lavFileName, LAVFile &lavFile);
};


#endif
