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


BOOST_AUTO_TEST_SUITE (map_iterators)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_iterators_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  size_t  num_elements = 0;
  for (const_iterator iter = map.begin (); iter != map.end (); ++iter)
    ++num_elements;

  BOOST_CHECK_EQUAL (num_elements, map.size ());
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_iterators_2,
                               parameters, test_map_parameters_normal)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  size_t  num_elements = 0;
  for (const_iterator iter = map.end (); iter != map.begin (); --iter)
    ++num_elements;

  BOOST_CHECK_EQUAL (num_elements, map.size ());
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_iterator_composite_1,
                               parameters, test_map_parameters_normal)
{
  COMMON_MAP_TEST_SETUP;

  map_type        map (RANGE (value_data::values1 ()));
  size_t          size = map.size ();
  const_iterator  scan = map.begin ();

  for (size_t k = 0; k != size; ++k, ++scan)
    {
      BOOST_CHECK (contains (map, scan->first));

      if (k != 0)
        {
          const_iterator  copy (scan);
          --copy;
          ++copy;
          BOOST_CHECK_EQUAL (copy, scan);

          const_iterator  copy2 (scan);
          copy2--;
          copy2++;
          BOOST_CHECK_EQUAL (copy2, scan);
        }

      const_iterator  copy (scan);
      ++copy;
      --copy;
      BOOST_CHECK_EQUAL (copy, scan);

      const_iterator  copy2 (scan);
      copy2++;
      copy2--;
      BOOST_CHECK_EQUAL (copy2, scan);
    }

  BOOST_CHECK_EQUAL (scan, map.end ());
}


// Following tests assume MCT_CHECK_PRECONDITIONS mode.

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_invalid_iterators_1,
                               parameters, test_map_parameters_robust_iterator_validation)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  iterator  end = map.end ();
  BOOST_CHECK_THROW (++end, logic_error);

  key_type  key;
  BOOST_CHECK_THROW (key = end->first, logic_error);
  MCT_UNUSED (key);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_invalid_iterators_2,
                               parameters, test_map_parameters_normal_robust_iterator_validation)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  iterator  begin = map.begin ();
  BOOST_CHECK_THROW (--begin, logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_erase_iterator_precondition_1,
                               parameters, test_map_parameters_normal)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));
  iterator  it = map.begin ();

  // This validation doesn't require fully robust checking to be available.
  map.erase (it);
  BOOST_CHECK_THROW (map.erase (it), logic_error);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_modification_through_iterator_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;

  map.insert (make_pair (key_data::values1 () [0], mapped_type ()));
  map.assign_through_iterator (map.find (key_data::values1 () [0]), mapped_data::values1 () [0]);

  BOOST_CHECK_EQUAL (map.at (key_data::values1 () [0]), mapped_data::values1 () [0]);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_before_begin_before_end_iterators_1,
                               parameters, test_map_parameters_forward)
{
  COMMON_MAP_TEST_SETUP;

  // While incrementing the before-* iterators is not particularly useful, it should still
  // work.  We test both empty and non-empty maps.

  map_type  map1;

  BOOST_CHECK_EQUAL (map1.before_begin (),   map1.before_end ());
  BOOST_CHECK_NE    (map1.before_begin (),   map1.end ());
  BOOST_CHECK_EQUAL (++map1.before_begin (), map1.begin ());
  BOOST_CHECK_EQUAL (++map1.before_end (),   map1.end ());
  BOOST_CHECK_EQUAL (++map1.before_begin (), ++map1.before_end ());

  map_type  map2 (RANGE (value_data::values1 ()));

  BOOST_CHECK_NE    (map2.before_begin (),   map2.before_end ());
  BOOST_CHECK_NE    (map2.before_begin (),   map2.end ());
  BOOST_CHECK_EQUAL (++map2.before_begin (), map2.begin ());
  BOOST_CHECK_EQUAL (++map2.before_end (),   map2.end ());
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
