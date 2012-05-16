// This file is part of Miscellaneous Container Templates.
//
//             https://launchpad.net/libmct/
//
// Copyright (c) 2009, 2010, 2011, 2012 Paul Pogonyshev.
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


#ifndef MCT_IMPL_OPTIMIZATION_HPP
#define MCT_IMPL_OPTIMIZATION_HPP


#ifdef __GNUC__


# define MCT_OPTIMIZATION_LIKELY(condition)             \
  __builtin_expect ((condition) ? true : false, true)

# define MCT_OPTIMIZATION_UNLIKELY(condition)           \
  __builtin_expect ((condition) ? true : false, false)

# define MCT_OPTIMIZATION_EXPECT(expression, value)     \
  __builtin_expect ((expression), (value))


#else  // not defined __GNUC__


# define MCT_OPTIMIZATION_LIKELY(condition)             (condition)
# define MCT_OPTIMIZATION_UNLIKELY(condition)           (condition)
# define MCT_OPTIMIZATION_EXPECT(expression, value)     (expression)


#endif


#endif  // Multi-inclusion guard.


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
