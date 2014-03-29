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

// Author: Patrick Marks, David Alexander

#pragma once

#include <string>
#include <vector>

#include "Types.hpp"

namespace ConsensusCore
{
    enum MutationType
    {
        INSERTION, DELETION, SUBSTITUTION
    };

    /// \brief Single mutation to a template sequence.
    class Mutation
    {
    private:
        MutationType type_;
        int start_;
        int end_;
        std::string newBases_;

        bool CheckInvariants() const;

    public:
        Mutation(MutationType type, int start, int end, std::string newBases);
        Mutation(MutationType type, int position, char base);

        MutationType Type() const;

        bool IsSubstitution() const;
        bool IsInsertion() const;
        bool IsDeletion() const;

        /// \brief Template positions of the mutation.
        /// Convention: the bases of the template changed by the mutation are [ Start, End ).
        //  For a substitution, tpl[Start..End) are mutated; for a deletion, tpl[Start..End) are removed.
        /// In the case of an insertion, Start=End= template position before the new bases are inserted.
        int Start() const;
        int End() const;

        std::string NewBases() const;
        int LengthDiff() const;
        std::string ToString() const;

    public:
        bool operator==(const Mutation& other) const;
        bool operator<(const Mutation& other) const;
    };


    std::string ApplyMutation(const Mutation& mut, const std::string& tpl);
    std::string ApplyMutations(const std::vector<Mutation*>& muts, const std::string& tpl);

    std::string MutationsToTranscript(const std::vector<Mutation*>& muts,
                                      const std::string& tpl);

    std::vector<int> TargetToQueryPositions(const std::vector<Mutation*>& mutations,
                                            const std::string& tpl);
}

#include "Mutation-inl.hpp"
