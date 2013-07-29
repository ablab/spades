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
#include <boost/range.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>
#include <cassert>
#include <string>
#include <vector>

#include "Types.hpp"

namespace ConsensusCore
{
    // Feature/Features object usage caveats:
    //  - Feature and Features objects _must_ be stored by value, not reference
    //  - The underlying array must be allocated using new[]
    template <typename T>
    class Feature : private boost::shared_array<T>
    {
    public:
        // \brief Allocate a new feature object, copying content from ptr.
        Feature(const T* ptr, int length)
            : boost::shared_array<T>(new T[length]),
              length_(length)
        {
            assert(length >= 0);
            std::copy(ptr, ptr + length, get());
        }

        // \brief Allocate and zero-fill a new feature object of given length.
        explicit Feature(int length)
            : boost::shared_array<T>(new T[length]()),
              length_(length)
        {
            assert(length >= 0);
        }

        int Length() const
        {
            return length_;
        }

        const T& operator[](int i) const
        {
            return this->boost::shared_array<T>::operator[](i);
        }

        T& operator[](int i)
        {
            return this->boost::shared_array<T>::operator[](i);
        }

        T ElementAt(int i) const
        {
            return (*this)[i];
        }

    private:
        int length_;

#ifndef SWIG
    public:
        T* get()
        {
            return this->boost::shared_array<T>::get();
        }

        const T* get() const
        {
            return this->boost::shared_array<T>::get();
        }

        operator std::string() const;
#endif  // !SWIG
    };


#ifndef SWIG
    //
    // Support for boost::foreach
    //
    template<typename T>
    inline const T* range_begin(const Feature<T>& f)
    {
        return f.get();
    }

    template<typename T>
    inline const T* range_end(const Feature<T>& f)
    {
        return f.get() + f.Length();
    }

    template<typename T>
    inline T* range_begin(Feature<T>& f) // NOLINT
    {
        return f.get();
    }

    template<typename T>
    inline T* range_end(Feature<T>& f)  // NOLINT
    {
        return f.get() + f.Length();
    }
#endif  // !SWIG

    typedef Feature<float> FloatFeature;
    typedef Feature<char> CharFeature;
    typedef Feature<int> IntFeature;
}


#ifndef SWIG
namespace boost
{
    template<typename T>
    struct range_const_iterator<ConsensusCore::Feature<T> >
    {
        typedef const T* type;
    };

    template<typename T>
    struct range_mutable_iterator<ConsensusCore::Feature<T> >
    {
        typedef T* type;
    };
}
#endif  // !SWIG
