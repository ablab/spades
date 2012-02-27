/***************************************************************************
 * Title:          GraphReader.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef GRAPH_READER_H_
#define GRAPH_READER_H_
#include "IntervalGraph.h"
#include "utils.h"
#include <iostream>

ssize_t GetDimensions(std::string &graphFileName, ssize_t &numVertices, ssize_t &numEdges);
ssize_t ReadGraph(std::string &graphFileName, IntervalGraph &graph);

#endif
