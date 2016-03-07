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

#include <algorithm>
#include <cassert>
#include <utility>

#include "Matrix/SparseMatrix.hpp"
#include "LFloat.hpp"

using std::min;
using std::max;

namespace ConsensusCore {
    //
    // Nullability
    //
    inline const SparseMatrix&
    SparseMatrix::Null()
    {
        static SparseMatrix* nullObj = new SparseMatrix(0, 0);
        return *nullObj;
    }

    inline bool
    SparseMatrix::IsNull() const
    {
        return (Rows() == 0 && Columns() == 0);
    }

    //
    // Size information
    //
    inline const int
    SparseMatrix::Rows() const
    {
        return nRows_;
    }

    inline const  int
    SparseMatrix::Columns() const
    {
        return nCols_;
    }

    //
    // Entry range queries per column
    //
    inline void
    SparseMatrix::StartEditingColumn(int j, int hintBegin, int hintEnd)
    {
        assert(columnBeingEdited_ == -1);
        columnBeingEdited_ = j;
        if (columns_[j] != NULL)
        {
            columns_[j]->ResetForRange(hintBegin, hintEnd);
        } else {
            columns_[j] = new SparseVector(Rows(), hintBegin, hintEnd);
        }
    }

    inline void
    SparseMatrix::FinishEditingColumn(int j, int usedRowsBegin, int usedRowsEnd)
    {
        assert(columnBeingEdited_ == j);
        usedRanges_[j] = std::make_pair(usedRowsBegin, usedRowsEnd);
        DEBUG_ONLY(CheckInvariants(columnBeingEdited_));
        columnBeingEdited_ = -1;
    }

    inline std::pair<int, int>
    SparseMatrix::UsedRowRange(int j) const
    {
        return usedRanges_[j];
    }

    inline bool
    SparseMatrix::IsColumnEmpty(int j) const
    {
        return (usedRanges_[j].first >= usedRanges_[j].second);
    }

    //
    // Accessors
    //
    inline const float&
    SparseMatrix::operator() (int i, int j) const
    {
        static const float emptyCell = Zero<lfloat>();
        if (columns_[j] == NULL)
        {
            return emptyCell;
        }
        else
        {
            return (*columns_[j])(i);
        }
    }

    inline float
    SparseMatrix::Get(int i, int j) const
    {
        return (*this)(i, j);
    }

    inline void
    SparseMatrix::Set(int i, int j, float v)
    {
        columns_[j]->Set(i, v);
    }

    inline void
    SparseMatrix::ClearColumn(int j)
    {
        usedRanges_[j] = std::make_pair(0, 0);
        columns_[j]->Clear();
        DEBUG_ONLY(CheckInvariants(j);)
    }

    //
    // SSE
    //
    inline __m128
    SparseMatrix::Get4(int i, int j) const
    {
        if (columns_[j] == NULL)
        {
            return Zero4<lfloat>();
        }
        else
        {
            return columns_[j]->Get4(i);
        }
    }

    inline void
    SparseMatrix::Set4(int i, int j, __m128 v4)
    {
        columns_[j]->Set4(i, v4);
    }
}
