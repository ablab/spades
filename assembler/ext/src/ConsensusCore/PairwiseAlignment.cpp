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

#include "PairwiseAlignment.hpp"

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <string>
#include <vector>

#include "Types.hpp"
#include "Sequence.hpp"
#include "Utils.hpp"

namespace ConsensusCore {

    std::string PairwiseAlignment::Target() const
    {
        return target_;
    }

    std::string PairwiseAlignment::Query() const
    {
        return query_;
    }

    float PairwiseAlignment::Accuracy() const
    {
        return ((float)(Matches())) / Length();
    }

    std::string PairwiseAlignment::Transcript() const
    {
        return transcript_;
    }

    int PairwiseAlignment::Matches() const
    {
        return std::count(transcript_.begin(), transcript_.end(), 'M');
    }

    int PairwiseAlignment::Errors() const
    {
        return Length() - Matches();
    }

    int PairwiseAlignment::Mismatches() const
    {
        return std::count(transcript_.begin(), transcript_.end(), 'R');
    }

    int PairwiseAlignment::Insertions() const
    {
        return std::count(transcript_.begin(), transcript_.end(), 'I');
    }

    int PairwiseAlignment::Deletions() const
    {
        return std::count(transcript_.begin(), transcript_.end(), 'D');
    }

    int PairwiseAlignment::Length() const
    {
        return target_.length();
    }

    PairwiseAlignment::PairwiseAlignment(const std::string& target, const std::string& query)
        :target_(target),
         query_(query),
         transcript_(target_.length(), 'Z')
    {
        if (target_.length() != query_.length())
        {
            throw InvalidInputError();
        }
        for (unsigned int i = 0; i < target_.length(); i++) {
            char t = target_[i];
            char q = query_[i];
            char tr;

            if (t == '-' && q == '-') { throw InvalidInputError(); }
            else if (t == q)          { tr = 'M'; }
            else if (t == '-')        { tr = 'I'; }
            else if (q == '-')        { tr = 'D'; }
            else                      { tr = 'R'; } // NOLINT

            transcript_[i] = tr;
        }
    }

    NeedlemanWunschParams::NeedlemanWunschParams(float matchScore,
                                                 float mismatchScore,
                                                 float insertScore,
                                                 float deleteScore)
        : MatchScore(matchScore),
          MismatchScore(mismatchScore),
          InsertScore(insertScore),
          DeleteScore(deleteScore)
    {}

    NeedlemanWunschParams DefaultNeedlemanWunschParams()
    {
        return NeedlemanWunschParams(0.0, -1.0, -1.0, -1.0);
    }


    static inline float MAX3(float a, float b, float c)
    {
        return std::max((a), std::max((b), (c)));
    }

    static inline int ARGMAX3(float a, float b, float c)
    {
        if      (a >= b && a >= c) return 0;
        else if (b >= c)           return 1;
        else                       return 2;
    }

    PairwiseAlignment*
    Align(const std::string& target,
          const std::string& query,
          NeedlemanWunschParams params)
    {
        using boost::numeric::ublas::matrix;

        int I = query.length();
        int J = target.length();
        matrix<float> Score(I + 1, J + 1);

        Score(0, 0) = 0;
        for (int i = 1; i <= I; i++) { Score(i, 0) = i * params.InsertScore; }
        for (int j = 1; j <= J; j++) { Score(0, j) = j * params.DeleteScore; }
        for (int i = 1; i <= I; i++)
        {
            for (int j = 1; j <= J; j++)
            {
                bool isMatch = (query[i - 1] == target[j - 1]);
                Score(i, j) = MAX3(Score(i - 1, j - 1) + (isMatch ? params.MatchScore :
                                                                    params.MismatchScore),
                                   Score(i - 1, j)     + params.InsertScore,
                                   Score(i,     j - 1) + params.DeleteScore);
            }
        }

        // Traceback, build up reversed aligned query, aligned target
        std::string raQuery, raTarget;
        int i = I, j = J;
        while (i > 0 || j > 0)
        {
            int move;
            if (i == 0) {
                move = 2;  // only deletion is possible
            } else if (j == 0) {
                move = 1;  // only insertion is possible
            } else {
                bool isMatch = (query[i - 1] == target[j - 1]);
                move = ARGMAX3(Score(i - 1, j - 1) + (isMatch ? params.MatchScore :
                                                                params.MismatchScore),
                               Score(i - 1, j)     + params.InsertScore,
                               Score(i,     j - 1) + params.DeleteScore);
            }
            // Incorporate:
            if (move == 0)
            {
                i--;
                j--;
                raQuery.push_back(query[i]);
                raTarget.push_back(target[j]);
            }
            // Insert:
            else if (move == 1)
            {
                i--;
                raQuery.push_back(query[i]);
                raTarget.push_back('-');
            }
            // Delete:
            else if (move == 2)
            {
                j--;
                raQuery.push_back('-');
                raTarget.push_back(target[j]);
            }
        }

        return new PairwiseAlignment(Reverse(raTarget), Reverse(raQuery));
    }


    //
    //  Code for lifting target coordinates into query coordinates.
    //


    static bool addsToTarget(char transcriptChar)
    {
        return (transcriptChar == 'M' ||
                transcriptChar == 'R' ||
                transcriptChar == 'D');
    }

    static int targetLength(const std::string& alignmentTranscript)
    {
        return std::count_if(alignmentTranscript.begin(), alignmentTranscript.end(), addsToTarget);
    }

#ifndef NDEBUG
    static bool addsToQuery(char transcriptChar)
    {
        return (transcriptChar == 'M' ||
                transcriptChar == 'R' ||
                transcriptChar == 'I');
    }

    static int queryLength(const std::string& alignmentTranscript)
    {
        return std::count_if(alignmentTranscript.begin(), alignmentTranscript.end(), addsToQuery);
    }
#endif  // !NDEBUG


    // TargetPositionsInQuery:
    // * Returns a vector of targetLength(transcript) + 1, which,
    //   roughly speaking, indicates the positions in the query string of the
    //   the characters in the target, as induced by an alignment with the
    //   given transcript string.
    // * More precisely, given an alignment (T, Q, X)  (x=transcript),
    //   letting T[s, e) denote any slice of T,
    //    - [s',e') denote the subslice of indices of Q aligned to T[s, e),
    //    - ntp = NewTargetPositions(X)
    //   we have
    //      [s', e') = [ntp(s), ntp(e))
    //
    // * Ex:
    //     MMM -> 0123
    //     DMM -> 0012,  MMD -> 0122, MDM -> 0112
    //     IMM -> 123,   MMI -> 013,  MIM -> 023
    //     MRM, MIDM, MDIM -> 0123
    std::vector<int> TargetToQueryPositions(const std::string& transcript)
    {
        std::vector<int> ntp;
        ntp.reserve(targetLength(transcript) + 1);

        int targetPos = 0;
        int queryPos = 0;
        foreach (char c, transcript)
        {
            if (c == 'M' || c == 'R')
            {
                ntp.push_back(queryPos);
                targetPos++;
                queryPos++;
            }
            else if (c == 'D')
            {
                ntp.push_back(queryPos);
                targetPos++;
            }
            else if (c == 'I')
            {
                queryPos++;
            }
            else
            {
                ShouldNotReachHere();
            }
        }
        ntp.push_back(queryPos);

        assert((int)ntp.size() == targetLength(transcript) + 1);
        assert(ntp[targetLength(transcript)] == queryLength(transcript));
        return ntp;
    }

    std::vector<int> TargetToQueryPositions(const PairwiseAlignment& aln)
    {
        return TargetToQueryPositions(aln.Transcript());
    }
}
