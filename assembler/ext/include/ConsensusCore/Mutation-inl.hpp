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

#include "Mutation.hpp"

#include <string>

#include "Utils.hpp"

namespace ConsensusCore {

    inline
    Mutation::Mutation(MutationType type, int start, int end, std::string newBases)
        : type_(type),
          start_(start),
          end_(end),
          newBases_(newBases)
    {
        if (!CheckInvariants()) throw InvalidInputError();
    }

    inline
    Mutation::Mutation(MutationType type, int position, char base)
        : type_(type),
          start_(position)
    {
        if (type == INSERTION) {
            end_ = position;
        } else {
            end_ = position + 1;
        }
        newBases_ = (type == DELETION ? "" : std::string(1, base));
        if (!CheckInvariants()) throw InvalidInputError();
    }

    inline bool
    Mutation::CheckInvariants() const
    {
        if (!((type_ == INSERTION && (start_ == end_) && newBases_.length() > 0)  ||
              (type_ == DELETION  && (start_ < end_)  && newBases_.length() == 0) ||
              (type_ == SUBSTITUTION && (start_ < end_) && ((int)(newBases_.length()) == end_ - start_))))
        {
            return false;
        }
        foreach (char base, newBases_)
        {
            if (!(base == 'A' ||
                  base == 'C' ||
                  base == 'G' ||
                  base == 'T'))
                return false;
        }
        return true;
    }

    inline bool
    Mutation::IsSubstitution() const
    {
        return (type_ == SUBSTITUTION);
    }

    inline bool
    Mutation::IsInsertion() const
    {
        return (type_ == INSERTION);
    }

    inline bool
    Mutation::IsDeletion() const
    {
        return (type_ == DELETION);
    }

    inline int
    Mutation::Start() const
    {
        return start_;
    }

    inline int
    Mutation::End() const
    {
        return end_;
    }

    inline std::string
    Mutation::NewBases() const
    {
        return newBases_;
    }

    inline MutationType
    Mutation::Type() const
    {
        return type_;
    }

    inline int
    Mutation::LengthDiff() const
    {
        if (IsInsertion())
            return newBases_.length();
        else if (IsDeletion())
            return start_ - end_;
        else
            return 0;
    }

    inline bool
    Mutation::operator==(const Mutation& other) const
    {
        return (Start()    == other.Start() &&
                End()      == other.End()   &&
                Type()     == other.Type()  &&
                NewBases() == other.NewBases());
    }

    inline bool
    Mutation::operator<(const Mutation& other) const
    {
        if (Start() != other.Start()) { return Start() < other.Start(); }
        if (End()   != other.End())   { return End()   < other.End();   }
        if (Type()  != other.Type())  { return Type()  < other.Type();  }
        return NewBases() < other.NewBases();
    }
}
