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


BOOST_AUTO_TEST_SUITE (set_other)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_find_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  for (typename vector <type>::const_iterator scan = data::values1 ().begin ();
       scan != data::values1 ().end (); ++scan)
    BOOST_CHECK (contains (set, *scan));

  for (typename vector <type>::const_iterator scan = data::values2 ().begin ();
       scan != data::values2 ().end (); ++scan)
    {
      BOOST_CHECK_EQUAL (contains (set, *scan),
                         find (RANGE (data::values1 ()), *scan) != range_end (data::values1 ()));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_count_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  for (typename vector <type>::const_iterator scan = data::values1 ().begin ();
       scan != data::values1 ().end (); ++scan)
    BOOST_CHECK_EQUAL (set.count (*scan), 1u);

  for (typename vector <type>::const_iterator scan = data::values2 ().begin ();
       scan != data::values2 ().end (); ++scan)
    {
      BOOST_CHECK_EQUAL (set.count (*scan),
                         (find (RANGE (data::values1 ()), *scan) != range_end (data::values1 ())
                          ? 1u : 0u));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_equal_range_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  for (typename vector <type>::const_iterator scan = data::values1 ().begin ();
       scan != data::values1 ().end (); ++scan)
    {
      pair <const_iterator, const_iterator>  range (set.equal_range (*scan));
      BOOST_CHECK       (set.key_eq () (*range.first, *scan));
      BOOST_CHECK_EQUAL (++range.first, range.second);
    }

  for (typename vector <type>::const_iterator scan = data::values2 ().begin ();
       scan != data::values2 ().end (); ++scan)
    {
      pair <const_iterator, const_iterator>  range (set.equal_range (*scan));

      if (find (RANGE (data::values1 ()), *scan) == range_end (data::values1 ()))
        BOOST_CHECK_EQUAL (range.first, range.second);
      else
        {
          BOOST_CHECK       (set.key_eq () (*range.first, *scan));
          BOOST_CHECK_EQUAL (++range.first, range.second);
        }
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_comparison_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set1 (RANGE (data::values1 ()));
  set_type  set2 (RANGE (data::values1 ()));
  set_type  set3 (RANGE (data::values2 ()));
  set_type  set4 (RANGE (data::values2 ()));

  BOOST_CHECK_EQUAL (set1, set2);
  BOOST_CHECK_EQUAL (set3, set4);
  BOOST_CHECK_NE    (set1, set3);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_assignment_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set1 (RANGE (data::values1 ()));
  set_type  set2;

  BOOST_CHECK_EQUAL (set2.size (), 0u);

  set2 = set1;
  BOOST_CHECK_EQUAL (set2, set1);

  if (implementation::stable_order)
    assert_identical_order (set1, set2);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_self_assignment_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set1 (RANGE (data::values1 ()));
  set_type  set2 (RANGE (data::values1 ()));

  BOOST_CHECK_EQUAL (set1, set2);
  if (implementation::stable_order)
    assert_identical_order (set1, set2);

  set1 = set1;

  BOOST_CHECK_EQUAL (set1, set2);
  if (implementation::stable_order)
    assert_identical_order (set1, set2);
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_move_assignment_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set;

  set = std::move (set_type (RANGE (data::values1 ())));
  BOOST_CHECK_EQUAL (set, set_type (RANGE (data::values1 ())));
}

// See comment in 'hash_table_base::operator= (hash_table_base&&)' for rationale.
BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_self_move_assignment_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set1 (RANGE (data::values1 ()));
  set_type  set2 (RANGE (data::values1 ()));

  BOOST_CHECK_EQUAL (set1, set2);
  if (implementation::stable_order)
    assert_identical_order (set1, set2);

  set1 = std::move (set1);

  BOOST_CHECK_EQUAL (set1, set2);
  if (implementation::stable_order)
    assert_identical_order (set1, set2);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_initializer_assignment_1,
                               implementation, test_set_implementations)
{
  typedef  typename implementation::template resolve <int>::set_type  set_type;

  set_type  set { 1, 2, 3, 4 };
  set = { 4, 5, 6 };

  if (implementation::stable_order)
    assert_identical_order (set, vector <int> { 4, 5, 6 });
}

#endif  // MCT_CXX0X_SUPPORTED


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_swap_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set1 (RANGE (data::values1 ()));
  set_type  set2 (RANGE (data::values2 ()));
  set_type  set3 (RANGE (data::values1 ()));
  set_type  set4 (RANGE (data::values2 ()));

  BOOST_CHECK_EQUAL (set1, set3);
  BOOST_CHECK_EQUAL (set2, set4);

  if (implementation::stable_order)
    {
      assert_identical_order (set1, set3);
      assert_identical_order (set2, set4);
    }

  swap (set3, set4);

  BOOST_CHECK_EQUAL (set1, set4);
  BOOST_CHECK_EQUAL (set2, set3);

  if (implementation::stable_order)
    {
      assert_identical_order (set1, set4);
      assert_identical_order (set2, set3);
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_clear_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));
  set.clear ();
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_erase_value_1,
                               parameters, test_set_parameters_normal)
{
  COMMON_SET_TEST_SETUP;

  set_type  set1 (RANGE (data::values1 ()));
  set_type  set2 (range_begin (data::values1 ()) + 1, range_end (data::values1 ()));

  // Make sure funny comparison functions don't ruin this test.
  if (set1.size () != set2.size ())
    {
      BOOST_CHECK       (set1.erase (data::values1 () [0]));
      BOOST_CHECK_EQUAL (set1, set2);

      if (implementation::stable_order)
        assert_identical_order (set1, set2);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_erase_position_1,
                               parameters, test_set_parameters_normal)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  const size_t  num_passes_expected = set.size ();
  size_t        num_passes          = 0;

  for (iterator scan = set.begin (), end = set.end (); scan != end; ++num_passes)
    scan = set.erase (scan);

  BOOST_CHECK       (set.empty ());
  BOOST_CHECK_EQUAL (num_passes, num_passes_expected);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_erase_range_1,
                               parameters, test_set_parameters_normal)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  set.erase (set.begin (), set.end ());
  BOOST_CHECK (set.empty ());
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_growth_1,
                               parameters, test_set_parameters_std_equal)
{
  COMMON_SET_TEST_SETUP;

  // Note: fairly small, because testing after each step is slow.
  set_type  set;
  for (int k = 0; k < 200; ++k)
    set.insert (data::generate (k));

  BOOST_CHECK_EQUAL (set.size (), 200u);
  if (implementation::stable_order)
    assert_identical_order (set, create_vector (data::generate, set.size ()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_shrink_1,
                               parameters, test_set_parameters_std_equal_normal)
{
  COMMON_SET_TEST_SETUP;

  const size_t   num_items = 200;
  vector <type>  items (create_vector (data::generate, num_items));
  set_type       set   (items.begin (), items.end ());

  BOOST_CHECK_EQUAL (set.size (), num_items);

  const_iterator  almost_end = set.end ();
  for (int k = 0; k < 10; ++k)
    --almost_end;

  set.erase (set.begin (), almost_end);

  BOOST_CHECK_EQUAL (set.size (), 10u);
  if (implementation::stable_order)
    assert_identical_order (set, create_vector (data::generate, num_items, num_items - 10));

  set.rehash (32);

  BOOST_CHECK_EQUAL (set.bucket_count (), 32u);
  if (implementation::stable_order)
    assert_identical_order (set, create_vector (data::generate, num_items, num_items - 10));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_tiny_set_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  static  const float  load_factors[] = { 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99 };

  for (const float* scan = range_begin (load_factors); scan != range_end (load_factors); ++scan)
    {
      BOOST_TEST_CHECKPOINT ("testing with maximum load factor " << *scan);

      set_type  set (1);
      set.max_load_factor (*scan);

      for (int k = 0; k < 50; ++k)
        set.insert (data::generate (k));

      if (implementation::stable_order)
        assert_identical_order (set, create_vector (data::generate, 50));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_max_load_factor_decrease_1,
                               implementation, test_set_implementations)
{
  typedef  test_table <typename implementation::template resolve <int>::set_type>  set_type;

  set_type  set;

  // Make sure 'set' is full for the default maximum load factor.
  for (int k = 0; ; k++)
    {
      set_type  copy (set);
      copy.insert (k);

      if (copy.bucket_count () == set.bucket_count ())
        set.swap (copy);
      else
        break;
    }

  // Just something very small; implementation will clamp to an acceptable value.
  set.max_load_factor (0.0f);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_intrusive_bucket_size_1,
                               implementation, test_set_implementations)
{
  typedef  test_table <typename implementation::template resolve <ham>::set_type>  set_type;

  set_type  set;

  BOOST_CHECK       (!set_type::implementation_type::KEEPS_HASHES);
  BOOST_CHECK_EQUAL (sizeof (typename set_type::iterator), sizeof (ham*));
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_cross_class_comparison_1,
                               implementation, test_set_implementations)
{
  typedef  test_table <typename implementation::template resolve <int>::set_type>  set_type;

  vector <int>  values1 (create_vector (Data <int>::generate, 10));
  vector <int>  values2 (values1);

  values2[3] = 100;

  set_type  set1 (RANGE (values1));
  set_type  set2 (RANGE (values2));

  BOOST_CHECK_EQUAL (set1, closed_hash_set <int> (RANGE (values1)));
  BOOST_CHECK_NE    (set1, closed_hash_set <int> (RANGE (values2)));

  BOOST_CHECK_EQUAL (set1, linked_hash_set <int> (RANGE (values1)));
  BOOST_CHECK_NE    (set1, linked_hash_set <int> (RANGE (values2)));

  BOOST_CHECK_EQUAL (set1, huge_linked_hash_set <int> (RANGE (values1)));
  BOOST_CHECK_NE    (set1, huge_linked_hash_set <int> (RANGE (values2)));

  BOOST_CHECK_EQUAL (set1, forward_hash_set <int> (RANGE (values1)));
  BOOST_CHECK_NE    (set1, forward_hash_set <int> (RANGE (values2)));

  BOOST_CHECK_EQUAL (set1, huge_forward_hash_set <int> (RANGE (values1)));
  BOOST_CHECK_NE    (set1, huge_forward_hash_set <int> (RANGE (values2)));
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
