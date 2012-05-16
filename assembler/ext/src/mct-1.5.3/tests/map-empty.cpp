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


BOOST_AUTO_TEST_SUITE (map_empty)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_map_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;

  BOOST_CHECK       (map.empty ());
  BOOST_CHECK_EQUAL (map.size (), 0u);
  BOOST_CHECK       (!map.uses_dynamic_memory ());
  BOOST_CHECK_EQUAL (map.used_memory (), sizeof (map_base_type));
  BOOST_CHECK_EQUAL (map.begin  (), map.end  ());
  BOOST_CHECK_EQUAL (map.cbegin (), map.cend ());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_map_explicit_size_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (1000);

  BOOST_CHECK (!map.uses_dynamic_memory ());

  map.insert (make_pair (key_data::generate (0), mapped_data::generate (0)));

  BOOST_CHECK    (map.uses_dynamic_memory ());
  BOOST_CHECK_GE (map.bucket_count (), 1000);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_map_copy_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;
  map_type  map2 (map);

  BOOST_CHECK       (map2.empty ());
  BOOST_CHECK_EQUAL (map2.size (), 0u);
  BOOST_CHECK       (!map2.uses_dynamic_memory ());
  BOOST_CHECK_EQUAL (map2.begin  (), map2.end  ());
  BOOST_CHECK_EQUAL (map2.cbegin (), map2.cend ());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_map_assignment_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;
  map_type  map2;

  map = map2;

  BOOST_CHECK       (map.empty ());
  BOOST_CHECK_EQUAL (map.size (), 0u);
  BOOST_CHECK       (!map.uses_dynamic_memory ());
  BOOST_CHECK_EQUAL (map.begin  (), map.end  ());
  BOOST_CHECK_EQUAL (map.cbegin (), map.cend ());

  map.insert (RANGE (value_data::values1 ()));
  map = map2;

  // Check the same invariants if assigning over a non-empty map.
  BOOST_CHECK       (map.empty ());
  BOOST_CHECK_EQUAL (map.size (), 0u);
  BOOST_CHECK       (!map.uses_dynamic_memory ());
  BOOST_CHECK_EQUAL (map.begin  (), map.end  ());
  BOOST_CHECK_EQUAL (map.cbegin (), map.cend ());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_map_erase_1,
                               parameters, test_map_parameters_normal)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;

  for (typename vector <key_type>::const_iterator scan = key_data::values1 ().begin ();
       scan != key_data::values1 ().end (); ++scan)
    BOOST_CHECK_EQUAL (map.erase (*scan), 0u);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_map_rehash_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;

  BOOST_CHECK (!map.uses_dynamic_memory ());

  map.rehash (1000);

  BOOST_CHECK    (!map.uses_dynamic_memory ());
  BOOST_CHECK_GE (map.bucket_count (), 1000);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_map_rehash_cleared_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));
  map.clear ();

  BOOST_CHECK (map.uses_dynamic_memory ());

  map.rehash (2 * map.bucket_count ());

  BOOST_CHECK (!map.uses_dynamic_memory ());
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
