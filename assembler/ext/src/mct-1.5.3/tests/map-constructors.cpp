// This file is part of Miscellaneous Container Templates.
//
//             https://launchpad.net/libmct/
//
// Copyright (c) 2011, 2012 Paul Pogonyshev.
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

#include "tests/map-common.hpp"


BOOST_AUTO_TEST_SUITE (map_constructors)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_range_constructor_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));
  if (implementation::stable_order)
    assert_identical_order (map, value_data::keys1 ());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_map_range_constructor_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;
  map_type  map (value_data::values1 ().begin (), value_data::values1 ().begin ());
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_initializer_constructor_1,
                               implementation, test_map_implementations)
{
  typedef  typename implementation::template resolve <int, string>::map_type  map_type;

  map_type  map1 { { 1, "foo" }, { 2, "bar" }, { 3, "baz" }, { 4, "ham" } };
  map_type  map2 { { 1, "foo" }, { 3, "bar" }, { 1, "baz" }, { 3, "ham" }, { 2, "spam" },
                   { 2, "egg" }, { 2, "bop" } };

  map_type  expected { { 1, "foo" }, { 3, "bar" }, { 2, "spam" } };
  BOOST_CHECK_EQUAL (map2, expected);

  if (implementation::stable_order)
    {
      assert_identical_order (map1, vector <int> { 1, 2, 3, 4 });
      assert_identical_order (map2, vector <int> { 1, 3, 1, 3, 2, 2, 2 });
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_map_initializer_constructor_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;
  map_type  map ({ });
}

#endif  // MCT_CXX0X_SUPPORTED


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_copy_constructor_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map1 (RANGE (value_data::values1 ()));
  map_type  map2 (RANGE (value_data::values2 ()));
  map_type  map3 (map1);
  map_type  map4 (map2);

  BOOST_CHECK_EQUAL (map3, map1);
  BOOST_CHECK_EQUAL (map4, map2);

  if (implementation::stable_order)
    {
      assert_identical_order (map3, keys_of (map1));
      assert_identical_order (map4, keys_of (map2));
    }
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_move_constructor_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (std::move (map_type (RANGE (value_data::values1 ()))));
  BOOST_CHECK_EQUAL (map, map_type (RANGE (value_data::values1 ())));
}

#endif  // MCT_CXX0X_SUPPORTED


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
