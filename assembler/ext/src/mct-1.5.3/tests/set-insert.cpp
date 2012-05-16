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


BOOST_AUTO_TEST_SUITE (set_insert)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_insert_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set;

  for (typename vector <type>::const_iterator scan = data::values1 ().begin ();
       scan != data::values1 ().end (); ++scan)
    set.insert (*scan);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_insert_range_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set;
  set.insert (RANGE (data::values1 ()));
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_insert_moving_1,
                               implementation, test_set_implementations)
{
  {
    typedef  typename implementation::template resolve <must_not_be_copied>::set_type  set_type;

    // Make it large, so it never rehashes and thus copies elements.
    set_type  set (1000);

    for (int k = 0; k < 100; ++k)
      {
        set.insert (must_not_be_copied (k));
        set.validate_integrity ();
      }
  }

  {
    typedef  typename implementation::template resolve <must_not_be_moved>::set_type  set_type;

    set_type  set (1000);

    for (int k = 0; k < 100; ++k)
      {
        must_not_be_moved  value (k);

        set.insert (value);
        set.validate_integrity ();
      }
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_insert_initializers_1,
                               implementation, test_set_implementations)
{
  typedef  typename implementation::template resolve <int>::set_type  set_type;
  typedef  typename robust_for <set_type>::type                       robust_type;

  set_type  set { 0, 1, 2 };
  set.insert ({ 3, 4, 5 });

  vector <int>  expected_vector { 0, 1, 2, 3, 4, 5 };
  robust_type   expected (expected_vector.begin (), expected_vector.end ());

  set.validate_integrity ();
  BOOST_CHECK_EQUAL (set, expected);

  if (implementation::stable_order)
    assert_identical_order (set, expected);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_emplace_1,
                               implementation, test_set_implementations)
{
  typedef  typename implementation::template resolve <int_wrapper>::set_type  set_type;

  set_type  set;
  for (int k = 0; k < 50; ++k)
    {
      set.emplace (k);
      set.validate_integrity ();
    }

  BOOST_CHECK_EQUAL (set.size (), 50u);
  if (implementation::stable_order)
    assert_identical_order (set, create_vector (int_wrapper::wrap, set.size ()));
}

#endif  // MCT_CXX0X_SUPPORTED


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_insert_erase_1,
                               parameters, test_set_parameters_normal)
{
  COMMON_SET_TEST_SETUP;

  set_type  set;

  for (typename vector <type>::const_iterator scan = data::values1 ().begin ();
       scan != data::values1 ().end (); ++scan)
    {
      // Simply insert and immediately erase one element.
      set.insert (*scan);
      BOOST_CHECK (set.erase (*scan));
    }
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
