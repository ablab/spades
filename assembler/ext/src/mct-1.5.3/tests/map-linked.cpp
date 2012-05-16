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


BOOST_AUTO_TEST_SUITE (map_linked)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_front_1,
                               parameters, test_map_parameters_std_equal_stable_order)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));
  BOOST_CHECK_EQUAL (assignable_value_type (map.front ()), value_data::values1 () [0]);

  map.pop_front ();
  BOOST_CHECK_EQUAL (assignable_value_type (map.front ()), value_data::values1 () [1]);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_back_1,
                               parameters, test_map_parameters_std_equal_stable_order)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));
  BOOST_CHECK_EQUAL (assignable_value_type (map.back ()),
                     *(range_end (value_data::values1 ()) - 1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_back_2,
                               parameters, test_map_parameters_std_equal_linked)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));
  BOOST_CHECK_EQUAL (assignable_value_type (map.back ()),
                     *(range_end (value_data::values1 ()) - 1));

  map.pop_back ();
  BOOST_CHECK_EQUAL (assignable_value_type (map.back ()),
                     *(range_end (value_data::values1 ()) - 2));
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_erase_after_1,
                               parameters, test_map_parameters_forward)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  map.erase_after (map.before_begin ());

  // Erasing the last item.
  for (const_iterator scan = map.begin (); ;)
    {
      const_iterator  next = scan;
      ++next;

      if (next != map.before_end ())
        scan = next;
      else
        {
          map.erase_after (scan);
          break;
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_erase_after_range_1,
                               parameters, test_map_parameters_forward)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  // Erase everything but the first and the last items.
  map.erase_after (map.begin (), map.before_end ());

  BOOST_CHECK_EQUAL (map.size  (), 2u);
  BOOST_CHECK_EQUAL (assignable_value_type (map.front ()), value_data::values1 () [0]);
  BOOST_CHECK_EQUAL (assignable_value_type (map.back  ()),
                     value_data::values1 () [range_size (value_data::values1 ()) - 1u]);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_relink_1,
                               parameters, test_map_parameters_linked)
{
  COMMON_MAP_TEST_SETUP;

  vector <assignable_value_type>  initializers (create_vector (value_data::generate, 50));
  map_type                        map (initializers.begin (), initializers.end ());

  assert_identical_order (map, keys_of (initializers));

  map.relink (map.find (key_data::generate (40)), map.find (key_data::generate (20)));

  initializers.erase  (find (initializers.begin (), initializers.end (),
                             value_data::generate (20)));
  initializers.insert (find (initializers.begin (), initializers.end (),
                             value_data::generate (40)),
                       value_data::generate (20));

  assert_identical_order (map, keys_of (initializers));

  // The following "moves" have no effect.

  map.relink (map.find (key_data::generate (5)), map.find (key_data::generate (5)));

  assert_identical_order (map, keys_of (initializers));

  map.relink (map.find (key_data::generate (15)), map.find (key_data::generate (14)));

  assert_identical_order (map, keys_of (initializers));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_relink_after_1,
                               parameters, test_map_parameters_forward)
{
  COMMON_MAP_TEST_SETUP;

  vector <assignable_value_type>  initializers (create_vector (value_data::generate, 50));
  map_type                        map (initializers.begin (), initializers.end ());

  assert_identical_order (map, keys_of (initializers));

  map.relink_after (map.before_end (), map.before_begin ());

  value_type  first (initializers.front ());
  initializers.erase     (initializers.begin ());
  initializers.push_back (first);

  assert_identical_order (map, keys_of (initializers));

  // The following "moves" have no effect.

  map.relink_after (map.find (key_data::generate (5)), map.find (key_data::generate (5)));

  assert_identical_order (map, keys_of (initializers));

  map.relink_after (map.find (key_data::generate (15)), map.find (key_data::generate (14)));

  assert_identical_order (map, keys_of (initializers));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_reverse_1,
                               parameters, test_map_parameters_stable_order)
{
  COMMON_MAP_TEST_SETUP;

  vector <assignable_value_type>  initializers (create_vector (value_data::generate, 50));
  map_type                        map (initializers.begin (), initializers.end ());

  assert_identical_order (map, keys_of (initializers));

  map.reverse ();
  reverse (initializers.begin (), initializers.end ());

  assert_identical_order (map, keys_of (initializers));
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_sort_1,
                               implementation, test_map_implementations_stable_order)
{
  typedef  typename implementation::template resolve <int, int>::map_type  map_type;

  const vector <pair <int, int> >&  values = Data <pair <int, int> >::values3 ();

  for (size_t k = 0; k < values.size (); ++k)
    {
      vector <pair <int, int> >  first_values (values.begin (), values.begin () + k);
      map_type                   map          (RANGE (first_values));

      assert_identical_order (map, keys_of (first_values));

      std::sort (RANGE (first_values));
      map.sort ();
      assert_identical_order (map, keys_of (first_values));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_sort_stability_1,
                               implementation, test_map_implementations_stable_order)
{
  typedef  typename implementation::template resolve <int, int>::map_type  map_type;

  const vector <pair <int, int> >&  values = Data <pair <int, int> >::values3 ();

  for (size_t k = 0; k < values.size (); ++k)
    {
      vector <pair <int, int> >  first_values (values.begin (), values.begin () + k);
      map_type                   map          (RANGE (first_values));
      assert_identical_order (map, keys_of (first_values));

      std::stable_sort (RANGE (first_values), even_first ());
      map.sort (even_first ());
      assert_identical_order (map, keys_of (first_values));
    }
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
