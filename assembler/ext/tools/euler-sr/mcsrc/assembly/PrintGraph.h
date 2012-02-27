/***************************************************************************
 * Title:          PrintGraph.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "DeBruijnGraph.h"

#include <set>
#include <iostream>

void PrintGraphSubset(IntervalGraph &g, 
											std::set<ssize_t> &vertexSet,
											ssize_t minEdgeLength,
											std::ostream &graphOut);
