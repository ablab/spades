/***************************************************************************
 * Title:          FasterSearchForCycle.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
ssize_t IntervalGraph::FindShortContainedDAG(ssize_t sourceVertex, ssize_t maxEdgeLength, 
																				 std::vector<ssize_t> &optimalPath) {

