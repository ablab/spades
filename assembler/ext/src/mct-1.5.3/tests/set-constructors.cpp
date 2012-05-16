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


// Even though we perform our own integrity validation, and even comparisons to a known
// good implementation in our tests, we still let the tables validate themselves with
// MCT_SELF_VALIDATION.  If nothing else, this ensures that MCT_SELF_VALIDATION behaves
// correctly on valid tables.
#undef  MCT_SELF_VALIDATION
#define MCT_SELF_VALIDATION     1

#include "tests/set-common.hpp"


BOOST_AUTO_TEST_SUITE (set_constructors)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_range_constructor_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));
  if (implementation::stable_order)
    assert_identical_order (set, data::values1 ());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_set_range_constructor_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;
  set_type  set (data::values1 ().begin (), data::values1 ().begin ());
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_initializer_constructor_1,
                               implementation, test_set_implementations)
{
  typedef  typename implementation::template resolve <int>::set_type  set_type;

  set_type  set1 { 1, 2, 3, 4 };
  set_type  set2 { 1, 3, 1, 3, 2, 2, 2 };

  if (implementation::stable_order)
    {
      assert_identical_order (set1, vector <int> { 1, 2, 3, 4 });
      assert_identical_order (set2, vector <int> { 1, 3, 1, 3, 2, 2, 2 });
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_set_initializer_constructor_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;
  set_type  set ({ });
}

#endif  // MCT_CXX0X_SUPPORTED


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_copy_constructor_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set1 (RANGE (data::values1 ()));
  set_type  set2 (RANGE (data::values2 ()));
  set_type  set3 (set1);
  set_type  set4 (set2);

  BOOST_CHECK_EQUAL (set3, set1);
  BOOST_CHECK_EQUAL (set4, set2);

  if (implementation::stable_order)
    {
      assert_identical_order (set3, set1);
      assert_identical_order (set4, set2);
    }
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_move_constructor_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (std::move (set_type (RANGE (data::values1 ()))));
  BOOST_CHECK_EQUAL (set, set_type (RANGE (data::values1 ())));
}

#endif  // MCT_CXX0X_SUPPORTED


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
