/***************************************************************************
 * Title:          LAVContig.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _LAV_CONTIG
#define _LAV_CONTIG

class LAVContig {
public:
  std::string contigName;
  ssize_t start, end, strand, contig;
};

#endif
