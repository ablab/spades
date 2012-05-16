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


BOOST_AUTO_TEST_SUITE (map_other)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_find_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  for (typename vector <key_type>::const_iterator key_scan = value_data::keys1 ().begin ();
       key_scan != value_data::keys1 ().end (); ++key_scan)
    BOOST_CHECK (contains (map, *key_scan));

  for (typename vector <key_type>::const_iterator key_scan = key_data::values2 ().begin ();
       key_scan != key_data::values2 ().end (); ++key_scan)
    {
      BOOST_CHECK_EQUAL (contains (map, *key_scan),
                         (find (RANGE (value_data::keys1 ()), *key_scan)
                          != range_end (value_data::keys1 ())));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_count_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  for (typename vector <key_type>::const_iterator key_scan = value_data::keys1 ().begin ();
       key_scan != value_data::keys1 ().end (); ++key_scan)
    BOOST_CHECK_EQUAL (map.count (*key_scan), 1u);

  for (typename vector <key_type>::const_iterator key_scan = key_data::values2 ().begin ();
       key_scan != key_data::values2 ().end (); ++key_scan)
    {
      BOOST_CHECK_EQUAL (map.count (*key_scan),
                         ((find (RANGE (value_data::keys1 ()), *key_scan)
                           != range_end (value_data::keys1 ()))
                          ? 1u : 0u));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_equal_range_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  for (typename vector <key_type>::const_iterator key_scan = value_data::keys1 ().begin ();
       key_scan != value_data::keys1 ().end (); ++key_scan)
    {
      pair <const_iterator, const_iterator>  range (map.equal_range (*key_scan));
      BOOST_CHECK       (map.key_eq () (range.first->first, *key_scan));
      BOOST_CHECK_EQUAL (++range.first, range.second);
    }

  for (typename vector <key_type>::const_iterator key_scan = key_data::values2 ().begin ();
       key_scan != key_data::values2 ().end (); ++key_scan)
    {
      pair <const_iterator, const_iterator>  range (map.equal_range (*key_scan));

      if (find (RANGE (value_data::keys1 ()), *key_scan) == range_end (value_data::keys1 ()))
        BOOST_CHECK_EQUAL (range.first, range.second);
      else
        {
          BOOST_CHECK       (map.key_eq () (range.first->first, *key_scan));
          BOOST_CHECK_EQUAL (++range.first, range.second);
        }
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_assignment_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map1 (RANGE (value_data::values1 ()));
  map_type  map2;

  BOOST_CHECK_EQUAL (map2.size (), 0u);

  map2 = map1;
  BOOST_CHECK_EQUAL (map2, map1);

  if (implementation::stable_order)
    assert_identical_order (map1, keys_of (map2));
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_self_assignment_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map1 (RANGE (value_data::values1 ()));
  map_type  map2 (RANGE (value_data::values1 ()));

  BOOST_CHECK_EQUAL (map1, map2);
  if (implementation::stable_order)
    assert_identical_order (map1, keys_of (map2));

  map1 = map1;

  BOOST_CHECK_EQUAL (map1, map2);
  if (implementation::stable_order)
    assert_identical_order (map1, keys_of (map2));
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_move_assignment_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;

  map = std::move (map_type (RANGE (value_data::values1 ())));
  BOOST_CHECK_EQUAL (map, map_type (RANGE (value_data::values1 ())));
}

// See comment in 'hash_table_base::operator= (hash_table_base&&)' for rationale.
BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_self_move_assignment_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map1 (RANGE (value_data::values1 ()));
  map_type  map2 (RANGE (value_data::values1 ()));

  BOOST_CHECK_EQUAL (map1, map2);
  if (implementation::stable_order)
    assert_identical_order (map1, keys_of (map2));

  map1 = std::move (map1);

  BOOST_CHECK_EQUAL (map1, map2);
  if (implementation::stable_order)
    assert_identical_order (map1, keys_of (map2));
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_initializer_assignment_1,
                               implementation, test_map_implementations)
{
  typedef  typename implementation::template resolve <int, string>::map_type  map_type;

  map_type  map { { 1, "foo" }, { 2, "bar" }, { 3, "baz" }, { 4, "bop" } };
  map = { { 4, "ham" }, { 5, "spam" }, { 6, "egg" } };

  if (implementation::stable_order)
    assert_identical_order (map, vector <int> { 4, 5, 6 });
}

#endif  // MCT_CXX0X_SUPPORTED


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_swap_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map1 (RANGE (value_data::values1 ()));
  map_type  map2 (RANGE (value_data::values2 ()));
  map_type  map3 (RANGE (value_data::values1 ()));
  map_type  map4 (RANGE (value_data::values2 ()));

  BOOST_CHECK_EQUAL (map1, map3);
  BOOST_CHECK_EQUAL (map2, map4);

  if (implementation::stable_order)
    {
      assert_identical_order (map1, keys_of (map3));
      assert_identical_order (map2, keys_of (map4));
    }

  swap (map3, map4);

  BOOST_CHECK_EQUAL (map1, map4);
  BOOST_CHECK_EQUAL (map2, map3);

  if (implementation::stable_order)
    {
      assert_identical_order (map1, keys_of (map4));
      assert_identical_order (map2, keys_of (map3));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_clear_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));
  map.clear ();
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_erase_value_1,
                               parameters, test_map_parameters_normal)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map1 (RANGE (value_data::values1 ()));
  map_type  map2 (range_begin (value_data::values1 ()) + 1, range_end (value_data::values1 ()));

  // Make sure funny comparison functions don't ruin this test.
  if (map1.size () != map2.size ())
    {
      BOOST_CHECK       (map1.erase (value_data::values1 () [0].first));
      BOOST_CHECK_EQUAL (map1, map2);

      if (implementation::stable_order)
        assert_identical_order (map1, keys_of (map2));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_erase_position_1,
                               parameters, test_map_parameters_normal)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  const size_t  num_passes_expected = map.size ();
  size_t        num_passes          = 0;

  for (iterator scan = map.begin (), end = map.end (); scan != end; ++num_passes)
    scan = map.erase (scan);

  BOOST_CHECK       (map.empty ());
  BOOST_CHECK_EQUAL (num_passes, num_passes_expected);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_erase_range_1,
                               parameters, test_map_parameters_normal)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map (RANGE (value_data::values1 ()));

  map.erase (map.begin (), map.end ());
  BOOST_CHECK (map.empty ());
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_subscription_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map;
  int       k = 0;
  for (typename vector <key_type>::const_iterator key_scan = key_data::values1 ().begin ();
       key_scan != key_data::values1 ().end (); ++key_scan)
    map.assign_at_key (*key_scan, mapped_data::generate (k++));
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_comparison_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  map_type  map1 (RANGE (value_data::values1 ()));
  map_type  map2 (RANGE (value_data::values1 ()));
  map_type  map3 (RANGE (value_data::values2 ()));
  map_type  map4 (RANGE (value_data::values2 ()));

  BOOST_CHECK_EQUAL (map1, map2);
  BOOST_CHECK_EQUAL (map3, map4);
  BOOST_CHECK_NE    (map1, map3);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_growth_1,
                               parameters, test_map_parameters_std_equal)
{
  COMMON_MAP_TEST_SETUP;

  // Note: fairly small, because testing after each step is slow.
  map_type  map;
  for (int k = 0; k < 200; ++k)
    map.insert (make_pair (key_data::generate (k), mapped_data::generate (k)));

  BOOST_CHECK_EQUAL (map.size (), 200u);
  if (implementation::stable_order)
    assert_identical_order (map, create_vector (key_data::generate, map.size ()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_shrink_1,
                               parameters, test_map_parameters_std_equal_normal)
{
  COMMON_MAP_TEST_SETUP;

  const size_t                    num_items = 200;
  vector <assignable_value_type>  items (create_vector (value_data::generate, num_items));
  map_type                        map   (items.begin (), items.end ());

  BOOST_CHECK_EQUAL (map.size (), num_items);

  const_iterator  almost_end = map.end ();
  for (int k = 0; k < 10; ++k)
    --almost_end;

  map.erase (map.begin (), almost_end);

  BOOST_CHECK_EQUAL (map.size (), 10u);
  if (implementation::stable_order)
    assert_identical_order (map, create_vector (key_data::generate, num_items, num_items - 10));

  map.rehash (32);

  BOOST_CHECK_EQUAL (map.bucket_count (), 32u);
  if (implementation::stable_order)
    assert_identical_order (map, create_vector (key_data::generate, num_items, num_items - 10));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_tiny_map_1,
                               parameters, test_map_parameters)
{
  COMMON_MAP_TEST_SETUP;

  static  const float  load_factors[] = { 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99 };

  for (const float* scan = range_begin (load_factors); scan != range_end (load_factors); ++scan)
    {
      BOOST_TEST_CHECKPOINT ("testing with maximum load factor " << *scan);

      map_type  map (1);
      map.max_load_factor (*scan);

      for (int k = 0; k < 50; ++k)
        map.insert (value_data::generate (k));

      if (implementation::stable_order)
        assert_identical_order (map, create_vector (key_data::generate, 50));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_max_load_factor_decrease_1,
                               implementation, test_map_implementations)
{
  typedef  test_table <typename implementation::template resolve <int, int>::map_type>  map_type;

  map_type  map;

  // Make sure 'map' is full for the default maximum load factor.
  for (int k = 0; ; k++)
    {
      map_type  copy (map);
      copy.insert (make_pair (k, k));

      if (copy.bucket_count () == map.bucket_count ())
        map.swap (copy);
      else
        break;
    }

  // Just something very small; implementation will clamp to an acceptable value.
  map.max_load_factor (0.0f);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_intrusive_bucket_size_1,
                               implementation, test_map_implementations)
{
  {
    typedef  test_table <typename implementation::template resolve <ham, int>::map_type>  map_type;

    map_type  map;

    BOOST_CHECK       (!map_type::implementation_type::KEEPS_HASHES);
    BOOST_CHECK_EQUAL (sizeof (typename map_type::iterator), sizeof (pair <ham, int>*));
  }

  {
    typedef  test_table <typename implementation::template resolve <int, ham>::map_type>  map_type;

    map_type  map;

    BOOST_CHECK       (!map_type::implementation_type::KEEPS_HASHES);
    BOOST_CHECK_EQUAL (sizeof (typename map_type::iterator), sizeof (pair <int, ham>*));
  }

  {
    typedef  test_table <typename implementation::template resolve <ham, ham>::map_type>  map_type;

    map_type  map;

    BOOST_CHECK       (!map_type::implementation_type::KEEPS_HASHES);
    BOOST_CHECK_EQUAL (sizeof (typename map_type::iterator), sizeof (pair <ham, ham>*));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_cross_class_comparison_1,
                               implementation, test_map_implementations)
{
  typedef  test_table <typename implementation::template resolve <int, int>::map_type>  map_type;

  vector <pair <int, int> >  values1 (create_vector (Data <pair <int, int> >::generate, 10));
  vector <pair <int, int> >  values2 (values1);

  values2[3].second = 100;

  map_type  map1 (RANGE (values1));
  map_type  map2 (RANGE (values2));

  // Damn macros, can't have commas in arguments.
  typedef  closed_hash_map <int, int>  closed_hash_map_int;
  BOOST_CHECK_EQUAL (map1, closed_hash_map_int (RANGE (values1)));
  BOOST_CHECK_NE    (map1, closed_hash_map_int (RANGE (values2)));

  typedef  linked_hash_map <int, int>  linked_hash_map_int;
  BOOST_CHECK_EQUAL (map1, linked_hash_map_int (RANGE (values1)));
  BOOST_CHECK_NE    (map1, linked_hash_map_int (RANGE (values2)));

  typedef  huge_linked_hash_map <int, int>  huge_linked_hash_map_int;
  BOOST_CHECK_EQUAL (map1, huge_linked_hash_map_int (RANGE (values1)));
  BOOST_CHECK_NE    (map1, huge_linked_hash_map_int (RANGE (values2)));

  typedef  forward_hash_map <int, int>  forward_hash_map_int;
  BOOST_CHECK_EQUAL (map1, forward_hash_map_int (RANGE (values1)));
  BOOST_CHECK_NE    (map1, forward_hash_map_int (RANGE (values2)));

  typedef  huge_forward_hash_map <int, int>  huge_forward_hash_map_int;
  BOOST_CHECK_EQUAL (map1, huge_forward_hash_map_int (RANGE (values1)));
  BOOST_CHECK_NE    (map1, huge_forward_hash_map_int (RANGE (values2)));
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
