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

#include "tests/set-common.hpp"


BOOST_AUTO_TEST_SUITE (set_linked)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_front_1,
                               parameters, test_set_parameters_std_equal_stable_order)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));
  BOOST_CHECK_EQUAL (set.front (), data::values1 () [0]);

  set.pop_front ();
  BOOST_CHECK_EQUAL (set.front (), data::values1 () [1]);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_back_1,
                               parameters, test_set_parameters_std_equal_stable_order)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));
  BOOST_CHECK_EQUAL (set.back (), *(range_end (data::values1 ()) - 1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_back_2,
                               parameters, test_set_parameters_std_equal_linked)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));
  BOOST_CHECK_EQUAL (set.back (), *(range_end (data::values1 ()) - 1));

  set.pop_back ();
  BOOST_CHECK_EQUAL (set.back (), *(range_end (data::values1 ()) - 2));
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_erase_after_1,
                               parameters, test_set_parameters_forward)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  set.erase_after (set.before_begin ());

  // Erasing the last item.
  for (const_iterator scan = set.begin (); ;)
    {
      const_iterator  next = scan;
      ++next;

      if (next != set.before_end ())
        scan = next;
      else
        {
          set.erase_after (scan);
          break;
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_erase_after_range_1,
                               parameters, test_set_parameters_forward)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  // Erase everything but the first and the last items.
  set.erase_after (set.begin (), set.before_end ());

  BOOST_CHECK_EQUAL (set.size  (), 2u);
  BOOST_CHECK_EQUAL (set.front (), data::values1 () [0]);
  BOOST_CHECK_EQUAL (set.back  (), data::values1 () [range_size (data::values1 ()) - 1u]);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_relink_1,
                               parameters, test_set_parameters_linked)
{
  COMMON_SET_TEST_SETUP;

  vector <type>  initializers (create_vector (data::generate, 50));
  set_type       set (initializers.begin (), initializers.end ());

  assert_identical_order (set, initializers);

  set.relink (set.find (data::generate (40)), set.find (data::generate (20)));

  initializers.erase  (find (initializers.begin (), initializers.end (), data::generate (20)));
  initializers.insert (find (initializers.begin (), initializers.end (), data::generate (40)),
                       data::generate (20));

  assert_identical_order (set, initializers);

  // The following "moves" have no effect.

  set.relink (set.find (data::generate (5)), set.find (data::generate (5)));

  assert_identical_order (set, initializers);

  set.relink (set.find (data::generate (15)), set.find (data::generate (14)));

  assert_identical_order (set, initializers);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_relink_after_1,
                               parameters, test_set_parameters_forward)
{
  COMMON_SET_TEST_SETUP;

  vector <type>  initializers (create_vector (data::generate, 50));
  set_type       set (initializers.begin (), initializers.end ());

  assert_identical_order (set, initializers);

  set.relink_after (set.before_end (), set.before_begin ());

  type  first (initializers.front ());
  initializers.erase     (initializers.begin ());
  initializers.push_back (first);

  assert_identical_order (set, initializers);

  // The following "moves" have no effect.

  set.relink_after (set.find (data::generate (5)), set.find (data::generate (5)));

  assert_identical_order (set, initializers);

  set.relink_after (set.find (data::generate (15)), set.find (data::generate (14)));

  assert_identical_order (set, initializers);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_reverse_1,
                               parameters, test_set_parameters_stable_order)
{
  COMMON_SET_TEST_SETUP;

  vector <type>  initializers (create_vector (data::generate, 50));
  set_type       set (initializers.begin (), initializers.end ());

  assert_identical_order (set, initializers);

  set.reverse ();
  reverse (initializers.begin (), initializers.end ());

  assert_identical_order (set, initializers);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_sort_1,
                               implementation, test_set_implementations_stable_order)
{
  typedef  typename implementation::template resolve <int>::set_type  set_type;

  const vector <int>&  values = Data <int>::values3 ();

  for (size_t k = 0; k < values.size (); ++k)
    {
      vector <int>  first_values (values.begin (), values.begin () + k);
      set_type      set          (RANGE (first_values));
      assert_identical_order (set, first_values);

      std::sort (RANGE (first_values));
      set.sort ();
      assert_identical_order (set, first_values);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_sort_stability_1,
                               implementation, test_set_implementations_stable_order)
{
  typedef  typename implementation::template resolve <int>::set_type  set_type;

  const vector <int>&  values = Data <int>::values3 ();

  for (size_t k = 0; k < values.size (); ++k)
    {
      vector <int>  first_values (values.begin (), values.begin () + k);
      set_type      set          (RANGE (first_values));
      assert_identical_order (set, first_values);

      std::stable_sort (RANGE (first_values), even_first ());
      set.sort (even_first ());
      assert_identical_order (set, first_values);
    }
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
