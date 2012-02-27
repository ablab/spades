/***************************************************************************
 * Title:          LAVPrinter.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "LAVPrinter.h"

#include "LAVAlignedContig.h"

std::ostream&  LAVPrinter::PrintLAVFile(LAVFile &lavFile, std::ostream  &out) {
  out << "#:lav" << std::endl;
  out << "d {" << std::endl;
  if (lavFile.blastzOpts == "") 
    out << "\"No alignment parameters\"" << std::endl;
  else
    out << "\"" << lavFile.blastzOpts << "\"" << std::endl;
  out << "}" << std::endl;
  ssize_t i;
  for (i = 0; i < lavFile.alignments.size(); i++) {
    out << "#:lav" << std::endl;
    lavFile.alignments[i]->PrintAlignment(out);
  }
  out << "#:eof" << std::endl;
  return out;
}
