/***************************************************************************
 * Title:          LAVPrinter.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _LAV_PRINTER
#define _LAV_PRINTER

#include <iostream>
#include <fstream>
#include <ostream>
#include <string>

#include "LAVFile.h"
#include "utils.h"
class LAVPrinter {
public:
  static void PrintLAVFile(LAVFile &lavFile, std::string &outName) {
    std::ofstream out;
    openck(outName, out, std::ios::out);
    PrintLAVFile(lavFile, out);
  }
  static std::ostream& PrintLAVFile(LAVFile &lavFile, std::ostream  &out);
};

#endif
