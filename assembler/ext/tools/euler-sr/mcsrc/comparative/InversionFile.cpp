/***************************************************************************
 * Title:          InversionFile.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <ctype.h>

#include "ValidatedInversion.h"
#include "InversionFile.h"

#include "utils.h"
#include "mctypes.h"

ssize_t ReadInversionFile(std::string &invFileName, 
		      InversionMap &invMap,
		      StringSet &species) {
  std::ifstream inInv;
  ssize_t line = 0;
  std::string title;
  openck(invFileName, inInv);
  std::string refSpec, qrySpec, seq;
  while (inInv) {
    ++line;
    if (inInv.peek() == '\n')
      return 1;
    if (!isdigit(inInv.peek())) {
      inInv >> refSpec >> qrySpec >> seq;

      inInv.get(); // disccard newline 
      // parse a title line
      // an example of a line is: bzalign/armadillo.ENm001.fa.baboon.ENm001.fa.lav
      /*
      ssize_t curPos = 0;
      ssize_t nextPos;
      ssize_t tempPos;

      // Parse a title.
      
      if ( (tempPos = title.find( "/", curPos)) != title.npos)
	curPos = tempPos + 1;

      if ( (nextPos = title.find(".", curPos )) != title.npos ) {
	refSpec = title.substr(curPos, nextPos - curPos );
	curPos = nextPos + 1;
	if ((nextPos = title.find(".", curPos)) != title.npos ) {
	  seq = title.substr(curPos, nextPos - curPos );
	  curPos = title.find(".", curPos+1)+1; // skip file extension.
	  //	  curPos = nextPos + 4; // skip ".fa", bad way to do this
	  if ((nextPos = title.find(".", curPos)) != title.npos) {
	    qrySpec = title.substr(curPos, nextPos - curPos );
	  }
	}

      */
      //      std::cout << "read: " << refSpec << " qry: " << qrySpec << std::endl;
      if (species.find(refSpec) == species.end())
	species.insert(refSpec);
      if (species.find(qrySpec) == species.end())
	species.insert(qrySpec);

    }
    else {
      // parse a location line
      ssize_t tStart, tEnd, qStart, qEnd, strand;
      if (! (inInv >> tStart >> tEnd >> qStart >> qEnd >> strand)) {
	std::cout << "error reading invFile on line: " << line << std::endl;
	exit(1);
      }
      inInv.get();
      Inversion* inversion;
      inversion = new Inversion;
      inversion->species = qrySpec;
      inversion->sequence = seq;
      inversion->tStart = tStart;
      inversion->tEnd   = tEnd;
      inversion->qStart = qStart;
      inversion->qEnd   = qEnd;
      invMap[refSpec][qrySpec].push_back(inversion);
    }
  }
  return 1;
}
