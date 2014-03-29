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

#pragma once

#include <xmmintrin.h>
#include <utility>
#include <vector>

#include "LFloat.hpp"
#include "Types.hpp"
#include "Utils.hpp"

namespace ConsensusCore {

    class SparseVector
    {
    public:  // Constructor, destructor
        SparseVector(int logicalLength, int beginRow, int endRow);
        SparseVector(const SparseVector& other);
        ~SparseVector();

        // Ensures there is enough allocated storage to
        // hold entries for at least [beginRow, endRow) (plus padding);
        // clears existing entries.
        void ResetForRange(int beginRow, int endRow);

    public:
        const float& operator()(int i) const;
        float Get(int i) const;
        void Set(int i, float v);
        __m128 Get4(int i) const;
        void Set4(int i, __m128 v);
        void Clear();

    public:
        int AllocatedEntries() const;
        void CheckInvariants() const;

    private:
        // Expand the range of rows for which we have backing storage,
        // while preserving contents.  The arguments will become the
        // new allocated bounds, so caller should add padding if desired
        // before calling.
        void ExpandAllocated(int newAllocatedBegin, int newAllocatedEnd);

    private:
        std::vector<float>* storage_;

        // the "logical" length of the vector, of which only
        // a subset of entries are actually allocated
        int logicalLength_;

        // row numbers in the abstraction we are presenting
        int allocatedBeginRow_;
        int allocatedEndRow_;

        // analytics
        int nReallocs_;
    };
}

#include "Matrix/SparseVector-inl.hpp"
