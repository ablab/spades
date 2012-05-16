// This file is part of Miscellaneous Container Templates.
//
//             https://launchpad.net/libmct/
//
// Copyright (c) 2010, 2011, 2012 Paul Pogonyshev.
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


BOOST_AUTO_TEST_SUITE (map_insert)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_insert_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;
  for (typename vector <key_type>::const_iterator key_scan = key_data::values1 ().begin ();
       key_scan != key_data::values1 ().end (); ++key_scan)
    map.insert (make_pair (*key_scan, mapped_type ()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_insert_range_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;
  map.insert (RANGE (value_data::values1 ()));
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_insert_moving_1,
                               implementation, test_map_implementations)
{
  // Note that we use 'pair <const ..., ...>' in this test.  While inserting a pair
  // without constified first member works fine, this involves some hidden conversions and
  // those call copy or move constructors, thus making the test fail.  We only want to
  // test whether the containers themselves don't copy/move when they shouldn't, so the
  // usage is a little bit cumbersome.

  // FIXME: This part of the test is currently disabled, because it is not expected to
  //        pass.  The problem is that 'value_type' in a map is 'pair <const key_type,
  //        mapped_type>' and constness of the first field prevents move constructor from
  //        kicking in.  I guess it can be solved and eventually should be, but this is
  //        not very important.
  if (false)
  {
    typedef  typename implementation::template resolve <must_not_be_copied, int>::map_type
             map_type;

    // Make it large, so it never rehashes and thus copies elements.
    map_type  map1 (1000);

    for (int k = 0; k < 100; ++k)
      {
        map1.insert (pair <const must_not_be_copied, int> (must_not_be_copied (k), 0));
        map1.validate_integrity ();
      }

    map_type  map2 (1000);

    for (int k = 0; k < 100; ++k)
      {
        map2[must_not_be_copied (k)] = 0;
        map2.validate_integrity ();
      }
  }

  {
    typedef  typename implementation::template resolve <int, must_not_be_copied>::map_type
             map_type;

    map_type  map1 (1000);

    for (int k = 0; k < 100; ++k)
      {
        map1.insert (pair <const int, must_not_be_copied> (k, must_not_be_copied (k)));
        map1.validate_integrity ();
      }
  }

  {
    typedef  typename implementation::template resolve <must_not_be_moved, int>::map_type
             map_type;

    map_type  map1 (1000);

    for (int k = 0; k < 100; ++k)
      {
        must_not_be_moved                    key (k);
        pair <const must_not_be_moved, int>  value (key, 0);

        map1.insert (value);
        map1.validate_integrity ();
      }

    map_type  map2 (1000);

    for (int k = 0; k < 100; ++k)
      {
        must_not_be_moved  key (k);

        map2[key] = 0;
        map2.validate_integrity ();
      }
  }

  {
    typedef  typename implementation::template resolve <int, must_not_be_moved>::map_type
             map_type;

    map_type  map1 (1000);

    for (int k = 0; k < 100; ++k)
      {
        must_not_be_moved                    mapped (k);
        pair <const int, must_not_be_moved>  value (k, mapped);

        map1.insert (value);
        map1.validate_integrity ();
      }
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_insert_initializers_1,
                               implementation, test_map_implementations)
{
  typedef  typename implementation::template resolve <int, string>::map_type  map_type;
  typedef  typename robust_for <map_type>::type                               robust_type;

  map_type  map { { 0, "foo" }, { 1, "bar" }, { 2, "baz" } };
  map.insert ({ { 3, "ham" }, { 4, "spam" }, { 5, "egg" } });

  vector <pair <int, string> >  expected_vector { { 0, "foo" }, { 1, "bar"  }, { 2, "baz" },
                                                  { 3, "ham" }, { 4, "spam" }, { 5, "egg" } };
  robust_type                   expected (expected_vector.begin (), expected_vector.end ());

  map.validate_integrity ();
  BOOST_CHECK_EQUAL (map, expected);

  if (implementation::stable_order)
    assert_identical_order (map, keys_of (expected));
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_emplace_1,
                               implementation, test_map_implementations)
{
  typedef  typename implementation::template resolve <int_wrapper, int>::map_type  map_type;

  map_type  map;
  for (int k = 0; k < 50; ++k)
    {
      map.emplace (make_pair (k, k));
      map.validate_integrity ();
    }

  BOOST_CHECK_EQUAL (map.size (), 50u);
  if (implementation::stable_order)
    assert_identical_order (map, create_vector (int_wrapper::wrap, map.size ()));
}

#endif


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_insert_erase_1,
                               parameters, test_map_parameters_normal)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;

  for (typename vector <assignable_value_type>::const_iterator value_scan
         = value_data::values1 ().begin ();
       value_scan != value_data::values1 ().end (); ++value_scan)
    {
      // Simply insert and immediately erase one element.
      map.insert (*value_scan);
      BOOST_CHECK (map.erase (value_scan->first));
    }
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
