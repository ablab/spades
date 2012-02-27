/***************************************************************************
 * Title:          cut.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
										ssize_t tnextedge, tnextindex;
										ssize_t tcuredge, tcurindex;
										tnextedge = paths[readIndex][nextPathPos].edge;
										tnextindex = paths[readIndex][nextPathPos].index;
										tcuredge  = paths[readIndex][curPathPos].edge;
										tcurindex = paths[readIndex][curPathPos].index;
										std::cout << "repalced for length: " 
															<< edges[tnextedge].intervals[tnextindex].readPos - 
											( edges[tcuredge].intervals[tcurindex].readPos +
												edges[tcuredge].intervals[tcurindex].length - 
												vertices[edges[tcuredge].dest].vertexSize );
										std::cout << " replacing forward " << curPathPos + 1 << " ... " << nextPathPos-1
															<< " with " << altPathEdges.size() << std::endl;
										ssize_t tcl = 0;
										ssize_t t;
										for (t = 0; t < altPathEdges.size(); t++ ) {
											tcl += edges[altPathEdges[t]].length -
												vertices[edges[altPathEdges[t]].dest].vertexSize;
										}
										std::cout << " replaced length: " << tcl << std::endl;

and for the compliment path:
											ssize_t cnextedge, cnextindex;
											ssize_t ccuredge, ccurindex;
											cnextedge = paths[compReadIndex][compPathNext].edge;
											cnextindex = paths[compReadIndex][compPathNext].index;
											ccuredge  = paths[compReadIndex][compPathPos].edge;
											ccurindex = paths[compReadIndex][compPathPos].index;
											std::cout << " replaced rev length: "  
																<< edges[cnextedge].intervals[cnextindex].readPos -
												( edges[ccuredge].intervals[ccurindex].readPos + 
													edges[ccuredge].intervals[ccurindex].length -
													vertices[edges[ccuredge].dest].vertexSize );
											std::cout << "replacing comp " << compPathPos + 1 << " ... " << compPathNext-1
																<< " with " << compPathEdges.size() << std::endl;
											ssize_t crl = 0;
											ssize_t c;
											for (c = 0; c < compPathEdges.size(); c++) {
												crl += edges[compPathEdges[c]].length - 
													vertices[edges[compPathEdges[c]].dest].vertexSize;
											}
											std::cout << "the length of the replacement path is: " << crl << std::endl;

	
	// lastNewEdge is the last new edge that was spliced in
	// The path should not have any gaps in it, so that 
	// intervalPos_i + intervalLength_i - vertexSize = intervalPos_(i+1)
	// The new path may pass through edges that make this not true
	// so we simply modify the length of the last interval so that 
	// it is.
	// Since the two edges overlap by 'vertexSize', this is subtracted.
	//
	// If the sum on the left is too short, we simply increase the length
	// of the interval (cause an insertion).
	// If the sum is too large, rather than deleting the extra, we shift the path
	// positions for the rest of the path.
			ssize_t nextEdge;
			ssize_t nextEdgeInterval;
			ssize_t edgeIndex;
			ssize_t endPos;
			ssize_t newEdge;
			newEdge = newEdges[e];
			// Let the new edge index is 'edgeIndex', and the path in
			// that edge passes through 'intvIndex'.
			// Also let the next edge in the path is 'nextEdgeIndex'
			// The path paths[readIndex][e + pathStart] ends at position 'endPos' within read 'readIndex.
			// The next edge should begin (readPos) at 'endPos' - size(edges[edgeIndex].dest).
			// If not, adjust the length of the interval in 'edgeIndex' to be so.
			
			edgeIndex = newEdges[e];
			endPos  = edges[edgeIndex].intervals[lastInterval].readPos + 
				edges[edgeIndex].intervals[lastInterval].length - 
				vertices[edges[edgeIndex].dest].vertexSize;

			if (edges[nextEdge].intervals[nextEdgeInterval].readPos != endPos) {
				if (endPos - edges[nextEdge].intervals[nextEdgeInterval].readPos > 0) {
					// the inserted edge is a deletion, just add that to the endo of the interval
					edges[newEdge].intervals[lastInterval].length += 
						endPos - edges[nextEdge].intervals[nextEdgeInterval].readPos;
				}
				else {
					ssize_t i;
					ssize_t bp;
					ssize_t diff = edges[nextEdge].intervals[nextEdgeInterval].readPos - endPos;
					if (readIndex > (pathLengths.size() / 2)) 
						bp = readIndex - (pathLengths.size() / 2) ;
					else
						bp = (pathLengths.size() / 2) + readIndex;

					for (i = pathEnd + newSpace + 1; 
							 i < pathLengths[readIndex] + newSpace && 
								 paths[readIndex][i].edge != -1 &&
								 paths[readIndex][i].index != -1; i++) {
						nextEdge = paths[readIndex][i].edge;
						nextEdgeInterval = paths[readIndex][i].index;
						edges[nextEdge].intervals[nextEdgeInterval].readPos += diff;
					}
				}
				assert(edges[newEdge].intervals[lastInterval].length >= 0);
			}
		}
	}


From RerouteSimplePathIntervals

	
  /*
    while (inInterval < numInIntervals and
    outInterval < numOutIntervals) {

    while ( inInterval < numInIntervals and
    outInterval < numOutIntervals and
    (edges[inEdge].intervals[inInterval].read !=
    edges[outEdge].intervals[outInterval].read)) {
    if (edges[inEdge].intervals[inInterval].read <
    edges[outEdge].intervals[outInterval].read)
    inInterval++;
    else
    outInterval++;
    }
    
    // If an entire list has been processed, just continue
    if (inInterval >= numInIntervals or 
    outInterval >= numOutIntervals) {
    continue;
    }
    // Otherwise, have reached intervals that are the same
    // do essentially the same processing
    ssize_t intervalRemoved;
    while (inInterval < numInIntervals  and
    outInterval < numOutIntervals and
    (edges[inEdge].intervals[inInterval].read ==
    edges[outEdge].intervals[outInterval].read)) {
    // Look to see if the interval 'inInterval' is followed
    // immediately by 'outInterval'
    if (IsIntervalMarkedForRemoval(outEdge, outInterval)) {
    outInterval++;
    continue;
    }

    if (!IsIntervalMarkedForRemoval(inEdge, inInterval) and
    edges[inEdge].intervals[inInterval].readPos + 
    edges[inEdge].intervals[inInterval].length - 
    vertices[edges[inEdge].dest].vertexSize ==
    edges[outEdge].intervals[outInterval].readPos) {

    // If so, merge 'outInterval' into 'inInterval'
    // Since both intervals take the vertex into account
    // the vertex size is incorporated into the outInterval
    // length, so it should be subtracted off to avoid double
    // counting.
    edges[inEdge].intervals[inInterval].length += 
    edges[outEdge].intervals[outInterval].length - vertexSize;
    assert(edges[inEdge].intervals[inInterval].length >= 0);
    // This out interval has been merged into the in interval.
    // We could remove it from the list, but that's too much
    // memory management.  Just simply mark the start pos of the 
    // out interval as 0
    ssize_t read, pathPos;

    // Flag this interval as being merged
    edges[outEdge].intervals[outInterval].readPos = -1;

    read = edges[outEdge].intervals[outInterval].read;
    pathPos = edges[outEdge].intervals[outInterval].pathPos;
    MarkIntervalForRemoval(outEdge, outInterval);
    DeleteReadInterval(read, pathPos);
													 
    // Advance to the next read interval in both the in edge
    // and out edge

    ++inInterval;
    ++outInterval;
    // We will have to regrow the intervals on inEdge
    // by the number of intervals that start in 
    // outEdge.  
    ++numMerged;
    }
    else {
    // both 'inInterval' and 'outInterval' are for the same
    // read, but the read passes through the edges multiple 
    // times, so 'outInterval' does not immediately succeed 'inInterval'.
	
    // Move forward to the earliest starting interval
    if (edges[inEdge].intervals[inInterval].readPos < 
    edges[outEdge].intervals[outInterval].readPos) {
    inInterval++;
    }
    else {
    outInterval++;
    }
    }
    }
    }
  */


This was cut from SplitVertices.cpp as a way to unlink edges from split vertices.
				// Fix the edge from the src vertex to the dest by removing it
				ssize_t outEdge, outEdgeIndex;
				ssize_t outUnlinked, inUnlinked;
				outUnlinked = inUnlinked = 0;
				for (outEdgeIndex = (*splitVertices[vSrc][vertexSize - 1].value).FirstOut();
						 outEdgeIndex != (*splitVertices[vSrc][vertexSize - 1].value).EndOut();
				outEdgeIndex = (*splitVertices[vSrc][vertexSize - 1].value).NextOut(outEdgeIndex)) {
					outEdge = (*splitVertices[vSrc][vertexSize-1].value).out[outEdgeIndex];
					// This points to the 
					if (newEdges[outEdge].dest == vDest) {
						// 
						newEdges[outEdge].dest = -1;
						newEdges[outEdge].src  = -1;
						(*splitVertices[vSrc][vertexSize-1].value).out[outEdgeIndex] = -1;
						ssize_t inEdge, inEdgeIndex;
						for (inEdgeIndex = (*splitVertices[vDest][0].value).FirstIn();
								 inEdgeIndex != (*splitVertices[vDest][0].value).EndIn();
								 inEdgeIndex = (*splitVertices[vDest][0].value).NextIn(inEdgeIndex)){
							if ((*splitVertices[vDest][0].value).in[inEdgeIndex] == outEdge) {
								(*splitVertices[vDest][0].value).in[inEdgeIndex] = -1;
								inUnlinked = 1;
								break;
							}
						}
						// unlinked the out, job done.
						outUnlinked = 1;
						break;
					}
				}

				assert(outUnlinked);
				assert(inUnlinked);
