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
#include <string>
#include <utility>

#include <stdint.h>


//
// Forward declarations
//
namespace ConsensusCore {
    struct BandingOptions;
    class DenseMatrix;
    struct QuiverConfig;
    class PairwiseAlignment;
    class AnnotatedPairwiseAlignment;
    struct PoaConfig;
    struct QvModelParams;
    struct EdnaModelParams;
    struct ChannelSequenceFeatures;
    struct QvSequenceFeatures;
    struct SequenceFeatures;
    class SparseMatrix;
    class Mutation;
}

namespace ConsensusCore {
namespace detail {
    class ViterbiCombiner;
    class SumProductCombiner;
}}


//
// Exception types
//

namespace ConsensusCore {
    /// \brief Abstract base class for "error"-type exceptions.  Do
    ///        not catch these.
    class ErrorBase
    {
    public:
        virtual std::string Message() const throw() = 0;
        virtual ~ErrorBase() {}
    };


    /// \brief Abstract base class for exceptions, which user code
    ///        may safely catch.
    class ExceptionBase
    {
    public:
        virtual std::string Message() const throw() = 0;
        virtual ~ExceptionBase() {}
    };


    /// \brief An exception signaling an error in ConsensusCore's internal logic.
    class InternalError : public ErrorBase
    {
    public:
        InternalError()
            : msg_("Internal error encountered!")
        {}

        explicit InternalError(const std::string& msg)
            : msg_(msg)
        {}

        std::string Message() const throw()
        {
            return msg_;
        }

    private:
        std::string msg_;
    };

    class InvalidInputError : public ErrorBase
    {
    public:
        InvalidInputError()
            : msg_("Invalid input!")
        {}

        explicit InvalidInputError(const std::string& msg)
            : msg_(msg)
        {}

        std::string Message() const throw()
        {
            return msg_;
        }

    private:
        std::string msg_;
    };

    class NotYetImplementedException : public ErrorBase
    {
    public:
        std::string Message() const throw()
        {
            return "Feature not yet implemented";
        }
    };
}


typedef std::pair<uint32_t, uint32_t> Interval;
