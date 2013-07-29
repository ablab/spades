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

#include <boost/utility.hpp>
#include <string>
#include <utility>
#include <vector>

#include "Poa/PoaConfig.hpp"
#include "Poa/PoaGraph.hpp"

namespace ConsensusCore
{
    using boost::noncopyable;

    /// \brief A multi-sequence consensus obtained from a partial-order alignment
    class PoaConsensus : private noncopyable
    {
        std::string consensusSequence_;
        PoaConfig config_;
        std::vector< std::pair<Mutation*, float> >* variants_;
        PoaGraph* poaGraph_;
        float score_;

    public:
        explicit PoaConsensus(const PoaConfig& config);
        ~PoaConsensus();

    public:
        static const PoaConsensus* FindConsensus(const std::vector<std::string>& reads,
                                                 const PoaConfig& config);
        static const PoaConsensus* FindConsensus(const std::vector<std::string>& reads,
                                                 bool global);
        static const PoaConsensus* FindConsensus(const std::vector<std::string>& reads);

    public:
        const PoaGraph* Graph() const;
        float Score() const;

        // Consensus interface methods
        std::string Sequence() const;
        std::string ToString() const;

        // Mutations interface
        const std::vector< std::pair<Mutation*, float> >* Mutations() const;
    };
}
