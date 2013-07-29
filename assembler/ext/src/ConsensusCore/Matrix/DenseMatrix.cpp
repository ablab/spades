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


#include "Matrix/DenseMatrix.hpp"

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <cassert>

#include "LFloat.hpp"

using boost::numeric::ublas::matrix;
using boost::numeric::ublas::row_major;

namespace ConsensusCore {

    // Performance insensitive routines are not inlined

    DenseMatrix::DenseMatrix(int rows, int cols)
        : boost_dense_matrix(rows, cols),
          usedRanges_(cols, std::make_pair(0, 0)),
          columnBeingEdited_(-1)
    {
        for (int j = 0; j < cols; j++)
        {
            CheckInvariants(j);
        }
    }

    DenseMatrix::~DenseMatrix()
    {}

    int
    DenseMatrix::UsedEntries() const
    {
        // use column ranges
        int filledEntries = 0;
        for (int col = 0; col < Columns(); ++col)
        {
            int start, end;
            boost::tie(start, end) = UsedRowRange(col);
            filledEntries += (end - start);
        }
        return filledEntries;
    }

    int
    DenseMatrix::AllocatedEntries() const
    {
        return Rows() * Columns();
    }

    void
    DenseMatrix::ToHostMatrix(float** mat, int* rows, int* cols) const
    {
        // TODO(dalexander): make sure SWIG client deallocates this memory -- use %newobject flag
        matrix<lfloat, row_major> rowMajorPeer(*this);
        *mat = new float[Rows() * Columns()];
        std::copy(rowMajorPeer.data().begin(), rowMajorPeer.data().end(), *mat);
        *rows = Rows();
        *cols = Columns();
    }

    void
    DenseMatrix::CheckInvariants(int column) const
    {
        // make sure no used entries are outside of the bands
        int start, end;
        boost::tie(start, end) = UsedRowRange(column);
        assert(0 <= start && start <= end && end <= Rows());
        for (int i = 0; i < Rows(); i++)
        {
            if (!(start <= i && i < end))
            {
                assert ((*this)(i, column) == value_type());
            }
        }
    }
}
