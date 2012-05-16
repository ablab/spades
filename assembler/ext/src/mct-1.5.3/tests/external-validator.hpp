// This file is part of Miscellaneous Container Templates.
//
//             https://launchpad.net/libmct/
//
// Copyright (c) 2012 Paul Pogonyshev.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#ifndef MCT_TESTS_EXTERNAL_VALIDATOR_HPP
#define MCT_TESTS_EXTERNAL_VALIDATOR_HPP


#include "common.hpp"


#if HAVE_BOOST_INTERPROCESS


typedef  pair <interp::offset_ptr <const void>, int>  external_validation_data;

// Can't use strings in interprocess context because they contain pointer(s).
class external_validation_error
{
  char  _message[0x1000];

public:

  external_validation_error (const string& message)
  {
    const size_t  length = min (message.length (), sizeof (_message) - 1);
    std::copy (message.data (), message.data () + length, _message);
    _message[length] = 0;
  }

  string
  message () const
  {  return _message;  }
};


#endif  // HAVE_BOOST_INTERPROCESS


#endif  // Multi-inclusion guard.


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
