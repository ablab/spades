// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: David Alexander

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "Coverage.hpp"

using std::vector;

namespace ConsensusCore {

    void CoverageInWindow(int tStartDim,
                          uint32_t *tStart,
                          int tEndDim,
                          uint32_t *tEnd,
                          uint32_t winStart,
                          int winLen,
                          uint32_t* coverage)
    {
        using std::max;
        using std::min;

        assert (tStartDim == tEndDim);

        int nReads = tStartDim;
        uint32_t winEnd = winStart + winLen;
        std::fill_n(coverage, winLen, 0);
        for (int read=0; read < nReads; read++)
        {
            uint32_t tStart_ = tStart[read];
            uint32_t tEnd_   = tEnd[read];
            for (uint32_t pos=max(tStart_, winStart); pos < min(tEnd_, winEnd); pos++)
            {
                coverage[pos-winStart] += 1;
            }
        }
    }


    #define CHUNK_SIZE 10000

    vector<Interval> CoveredIntervals(uint32_t minCoverage,
                                      int tStartDim,
                                      uint32_t *tStart,
                                      int tEndDim,
                                      uint32_t *tEnd,
                                      uint32_t winStart,
                                      int winLen)
    {
        using std::make_pair;

        assert (tStartDim == tEndDim);
        // assert(isSorted(tStart));  // find out how to get this ... it's C++11

        // Approach: divide into chunks, find coverage in each chunk,
        // then scan for covered windows ... careful to anneal windows
        // spanning chunk boundaries.  We also rely on the sortedness of the
        // tStart to restrict our attention to

        uint32_t winEnd = winStart + winLen;
        uint32_t coverage[CHUNK_SIZE];
        int currentIntervalStart = -1;
        vector<Interval> intervals;

        uint32_t startRowInChunk = 0;
        for (uint32_t chunkStart = winStart;
             chunkStart < winEnd;
             chunkStart += CHUNK_SIZE)
        {
            uint32_t chunkEnd = std::min(chunkStart + CHUNK_SIZE, winEnd);

            // We compute a conservative guess of the rows that are involved in
            // this chunk.  Not every row in the range [startRowInChunk, endRowInChunk)
            // actually overlaps the chunk, but no rows not in that range intersect the chunk.
            // startRowInChunk is computed by scanning from where it was for the last chunk.
            // This is the best we can do within pulling in the nBackRead stuff.
            uint32_t endRowInChunk = std::lower_bound(tStart, tStart + tStartDim, chunkEnd) - tStart;
            for (; ((startRowInChunk < endRowInChunk) & (tEnd[startRowInChunk] < chunkStart));
                 startRowInChunk++);

            CoverageInWindow((endRowInChunk-startRowInChunk), tStart+startRowInChunk,
                             (endRowInChunk-startRowInChunk), tEnd+startRowInChunk,
                             chunkStart, CHUNK_SIZE, coverage);
            int j = 0;
            while (j < chunkEnd - chunkStart)
            {
                if (coverage[j] >= minCoverage)
                {
                    if (currentIntervalStart == -1)
                    {
                        currentIntervalStart = chunkStart + j;
                    }
                }
                else
                {
                    if (currentIntervalStart != -1)
                    {
                        intervals.push_back(make_pair((uint32_t)currentIntervalStart,
                                                      chunkStart + j));
                        currentIntervalStart = -1;
                    }
                }
                j++;
            }
        }
        if (currentIntervalStart != -1)
        {
            intervals.push_back(make_pair((uint32_t)currentIntervalStart, winEnd));
        }
        return intervals;
    }
}
