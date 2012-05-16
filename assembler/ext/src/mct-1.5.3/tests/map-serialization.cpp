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


#include "tests/config.hpp"


#if HAVE_BOOST_SERIALIZATION


// Even though we perform our own integrity validation, and even comparisons to a known
// good implementation in our tests, we still let the tables validate themselves with
// MCT_SELF_VALIDATION.  If nothing else, this ensures that MCT_SELF_VALIDATION behaves
// correctly on valid tables.
#undef  MCT_SELF_VALIDATION
#define MCT_SELF_VALIDATION     1

#include "tests/map-common.hpp"
#include "tests/serialization-common.hpp"

#include <mct/hash-map-serialization.hpp>


BOOST_AUTO_TEST_SUITE (map_serialization)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_serialization_1,
                               implementation, test_map_implementations)
{
  typedef  typename implementation::template resolve <int, int>::map_type  map_type;

  for (size_t num_elements = 0; num_elements != 100; ++num_elements)
    {
      vector <pair <int, int> >  initializers
        (create_vector (Data <pair <int, int> >::generate, num_elements));

      map_type  map (initializers.begin (), initializers.end ());

      test_serialization (map);
    }
}


BOOST_AUTO_TEST_SUITE_END ()


#endif  // HAVE_BOOST_SERIALIZATION


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
