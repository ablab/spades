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

#include "Mutation.hpp"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <string>
#include <vector>

#include <iostream>

#include "Types.hpp"
#include "PairwiseAlignment.hpp"
#include "Utils.hpp"

using std::max;

namespace ConsensusCore
{
    std::string
    Mutation::ToString() const
    {
        using boost::str;
        using boost::format;

        switch (Type())
        {
            case INSERTION:    return str(format("Insertion (%s) @%d") % newBases_ % start_);
            case DELETION:     return str(format("Deletion @%d:%d") % start_ % end_);
            case SUBSTITUTION: return str(format("Substitution (%s) @%d:%d") % newBases_ % start_ % end_);
            default: ShouldNotReachHere();
        }
    }

    static void
    _ApplyMutationInPlace(const Mutation& mut, int start, std::string* tpl)
    {
        if (mut.IsSubstitution())
        {
            (*tpl).replace(start, mut.End() - mut.Start(), mut.NewBases());
        }
        else if (mut.IsDeletion())
        {
            (*tpl).erase(start, mut.End() - mut.Start());
        }
        else if (mut.IsInsertion())
        {
            (*tpl).insert(start, mut.NewBases());
        }
    }

    std::string
    ApplyMutation(const Mutation& mut, const std::string& tpl)
    {
        std::string tplCopy(tpl);
        _ApplyMutationInPlace(mut, mut.Start(), &tplCopy);
        return tplCopy;
    }

    struct compareMutationPointers {
        bool operator() (Mutation* lhs, Mutation* rhs) { return *lhs < *rhs; }
    };

    std::string
    ApplyMutations(const std::vector<Mutation*>& muts, const std::string& tpl)
    {
        std::string tplCopy(tpl);
        std::vector<Mutation*> sortedMuts(muts);
        std::sort(sortedMuts.begin(), sortedMuts.end(), compareMutationPointers());
        int runningLengthDiff = 0;
        foreach (const Mutation* mut, sortedMuts)
        {
            _ApplyMutationInPlace(*mut, mut->Start() + runningLengthDiff, &tplCopy);
            runningLengthDiff += mut->LengthDiff();
        }
        return tplCopy;
    }

    std::string MutationsToTranscript(const std::vector<Mutation*>& mutations,
                                      const std::string& tpl)
    {
        std::vector<Mutation*> sortedMuts(mutations);
        std::sort(sortedMuts.begin(), sortedMuts.end(), compareMutationPointers());

        // Build an alignnment transcript corresponding to these mutations.
        int tpos = 0;
        std::string transcript = "";
        foreach (const Mutation* m, sortedMuts)
        {
            for (; tpos < m->Start(); ++tpos)
            {
                transcript.push_back('M');
            }

            if (m->IsInsertion())
            {
                transcript += std::string(m->LengthDiff(), 'I');
            }
            else if (m->IsDeletion())
            {
                transcript += std::string(-m->LengthDiff(), 'D');
                tpos += -m->LengthDiff();
            }
            else if (m->IsSubstitution())
            {
                int len = m->End() - m->Start();
                transcript += std::string(len, 'R');
                tpos += len;
            }
            else
            {
                ShouldNotReachHere();
            }
        }
        for (; tpos < (int)tpl.length(); ++tpos)
        {
            transcript.push_back('M');
        }
        return transcript;
    }

    // MutatedTemplatePositions:
    //  * Returns a vector of length (tpl.length()+1), which, roughly speaking,
    //    indicates the positions in the mutated template tpl' of the characters
    //    in tpl.
    //  * More precisely, for any slice [s, e) of tpl, letting:
    //      - t[s, e) denote the induced substring of the template;
    //      - m[s, e) denote the subvector of mutations with Position
    //        in [s, e);
    //      - t' denote the mutated template; and
    //      - t[s, e)' denote the result of applying mutation m[s, e) to t[s, e),
    //    the resultant vector mtp satisfies t'[mtp[s], mtp[e]) == t[s,e)'.
    //  * Example:
    //               01234567                           0123456
    //              "GATTACA" -> (Del T@2, Ins C@5) -> "GATACCA";
    //    here mtp = 01223567, which makes sense, because for instance:
    //      - t[0,3)=="GAT" has become t'[0,2)=="GA";
    //      - t[0,2)=="GA"  remains "GA"==t'[0,2);
    //      - t[4,7)=="ACA" has become t[3,7)=="ACCA",
    //      - t[5,7)=="CA"  remains "CA"==t'[5,7).
    //
    std::vector<int> TargetToQueryPositions(const std::vector<Mutation*>& mutations,
                                            const std::string& tpl)
    {
        return TargetToQueryPositions(MutationsToTranscript(mutations, tpl));
    }
}
