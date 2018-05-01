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

#include <vector>
#include <string>
#include <utility>
#include <tuple>

#include "ConsensusCore/Types.hpp"
#include "ConsensusCore/Mutation.hpp"

namespace ConsensusCore
{
    /// \brief An object representing a Poa (partial-order alignment) graph
    class PoaGraph
    {
        class Impl;
        Impl* impl;

    public:  // Flags enums for specifying GraphViz output features
        enum {
            COLOR_NODES    = 0x1,
            VERBOSE_NODES = 0x2
        };

    public:
        void AddSequence(const std::string& sequence, const PoaConfig& config);

        // TODO(dalexander): move this method to PoaConsensus so we don't have to use a tuple
        // interface here (which was done to avoid a circular dep on PoaConsensus).
#ifndef SWIG
        std::tuple<std::string, float, std::vector< std::pair<Mutation*, float> >* >
        FindConsensus(const PoaConfig& config) const;
#endif  // !SWIG

        int NumSequences() const;
        std::string ToGraphViz(int flags = 0) const;
        void WriteGraphVizFile(std::string filename, int flags = 0) const;
        PoaGraph();
        ~PoaGraph();
    };
}

