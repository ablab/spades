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

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>
#include <cassert>
#include <utility>

#include "Matrix/DenseMatrix.hpp"
#include "Utils.hpp"


using boost::numeric::ublas::matrix_column;

namespace ConsensusCore {
    //
    // Nullability
    //
    inline const DenseMatrix&
    DenseMatrix::Null()
    {
        static DenseMatrix* nullObj = new DenseMatrix(0, 0);
        return *nullObj;
    }

    inline bool
    DenseMatrix::IsNull() const
    {
        return (Rows() == 0 && Columns() == 0);
    }

    //
    // Size information
    //
    inline const int
    DenseMatrix::Rows() const
    {
        return size1();
    }

    inline const  int
    DenseMatrix::Columns() const
    {
        return size2();
    }

    //
    // Entry range queries per column
    //
    inline void
    DenseMatrix::StartEditingColumn(int j, int hintBegin, int hintEnd)
    {
        assert(columnBeingEdited_ == -1);
        columnBeingEdited_ = j;
        ClearColumn(j);
    }

    inline void
    DenseMatrix::FinishEditingColumn(int j, int usedRowsBegin, int usedRowsEnd)
    {
        assert(columnBeingEdited_ == j);
        usedRanges_[j] = std::make_pair(usedRowsBegin, usedRowsEnd);
        DEBUG_ONLY(CheckInvariants(columnBeingEdited_));
        columnBeingEdited_ = -1;
    }

    inline std::pair<int, int>
    DenseMatrix::UsedRowRange(int j) const
    {
        return usedRanges_[j];
    }

    inline bool
    DenseMatrix::IsColumnEmpty(int j) const
    {
        return (usedRanges_[j].first >= usedRanges_[j].second);
    }

    //
    // Accessors
    //
    inline void
    DenseMatrix::Set(int i, int j, float v)
    {
        assert(columnBeingEdited_ == j);
        boost_dense_matrix::operator()(i, j) = v;
    }

    inline float
    DenseMatrix::Get(int i, int j) const
    {
        return (*this)(i, j);
    }

    inline const float&
    DenseMatrix::operator() (int i, int j) const
    {
        return boost_dense_matrix::operator()(i, j);
    }

    inline void
    DenseMatrix::ClearColumn(int j)
    {
        DEBUG_ONLY(CheckInvariants(j);)
        // (Rely on the fact that the underlying memory is stored
        // contiguously)
        int begin, end;
        boost::tie(begin, end) = usedRanges_[j];
        std::fill_n((float*)&boost_dense_matrix::operator()(begin, j),  // NOLINT
                    end - begin,
                    value_type());
        usedRanges_[j] = std::make_pair(0, 0);
        DEBUG_ONLY(CheckInvariants(j);)
    }

    //
    // SSE
    //
    inline __m128
    DenseMatrix::Get4(int i, int j) const
    {
        assert(0 <= i && i <= Rows() - 4);
        return _mm_loadu_ps(&boost_dense_matrix::operator()(i, j).value);
    }

    inline void
    DenseMatrix::Set4(int i, int j, __m128 v4)
    {
        assert(columnBeingEdited_ == j);
        assert(0 <= i && i <= Rows() - 4);
        _mm_storeu_ps(&boost_dense_matrix::operator()(i, j).value, v4);
    }
}

