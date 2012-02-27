/***************************************************************************
 * Title:          EmbossAlign.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/25/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <unistd.h>

#include "compatibility.h"
#include "EmbossAlign.h"
#include "utils.h"


ssize_t EmbossAlign(std::string &commandLine, 
		DNASequence &refSeq, DNASequence &qrySeq, 
		EmbossAlignment &alignment) {

  pid_t pid = getpid();
  char pidstr[20];
  sprintf(pidstr, PRI_PID, pid);
  
  std::string refFileName, qryFileName, alignFileName; 
  MakeTempName(refFileName, "embAlignRef.fasta");
  MakeTempName(qryFileName, "embAlignQuery.fasta");
  MakeTempName(alignFileName, "embAlignOut");

  std::ofstream refOut, qryOut;
  openck(refFileName, refOut, std::ios::out);
  refSeq.PrintSeq(refOut);
  refOut << std::endl;
  refOut.close();
  openck(qryFileName, qryOut, std::ios::out);
  qrySeq.PrintSeq(qryOut);
  qryOut << std::endl;  qryOut.close(); 
  
  commandLine += refFileName + " " + qryFileName + " -outfile=" + alignFileName + " >& /dev/null ";
  system((const char*) commandLine.c_str());
  
  ssize_t retval;
  retval = SRSPairParser::ParseSRSFile(alignFileName, alignment);
  
  std::string rmQry, rmRef, rmAlign;
  rmQry   = "rm " + qryFileName;
  rmRef   = "rm " + refFileName;
  rmAlign = "rm " + alignFileName;
  
  system((const char*) rmQry.c_str());
  system((const char*) rmRef.c_str());
  system((const char*) rmAlign.c_str());

  return retval;
}
