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


// MCT_SELF_VALIDATION will break tests with throwing hasher or comparator.
#undef MCT_SELF_VALIDATION


#include "tests/map-common.hpp"

#include <boost/mpl/end.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/vector.hpp>


namespace
{

  template <typename Implementation, typename Key, typename Mapped, bool _keep_hashes = false>
  struct MapExceptionSafetyParameters
  {
    typedef  Implementation  implementation;
    typedef  Key             key_type;
    typedef  Mapped          mapped_type;

    static  const bool  keep_hashes  = _keep_hashes;
    static  const bool  stable_order = Implementation::stable_order;
    static  const bool  normal       = Implementation::normal;

    // Not 'test_map' here as we can't control where exception are actually thrown.
    typedef  typename Implementation::template resolve <key_type, mapped_type,
                                                        MCT_HASH_NAMESPACE::hash <key_type>,
                                                        equal_to <key_type>,
                                                        allocator <pair <const key_type,
                                                                         mapped_type> >,
                                                        keep_hashes>::map_type
             tested_type;
  };

  typedef  mpl::vector
             <MapExceptionSafetyParameters <PlainMap,       throws_on_copy,       int, false>,
              MapExceptionSafetyParameters <PlainMap,       throws_on_copy,       int, true>,
              MapExceptionSafetyParameters <PlainMap,       throws_in_comparator, int, false>,
              MapExceptionSafetyParameters <PlainMap,       throws_in_comparator, int, true>,
              MapExceptionSafetyParameters <PlainMap,       throws_in_hasher,     int, false>,
              MapExceptionSafetyParameters <PlainMap,       throws_in_hasher,     int, true>,
              MapExceptionSafetyParameters <LinkedMap,      throws_on_copy,       int, false>,
              MapExceptionSafetyParameters <LinkedMap,      throws_on_copy,       int, true>,
              MapExceptionSafetyParameters <LinkedMap,      throws_in_comparator, int, false>,
              MapExceptionSafetyParameters <LinkedMap,      throws_in_comparator, int, true>,
              MapExceptionSafetyParameters <LinkedMap,      throws_in_hasher,     int, false>,
              MapExceptionSafetyParameters <LinkedMap,      throws_in_hasher,     int, true>,
              MapExceptionSafetyParameters <HugeLinkedMap,  throws_on_copy,       int, false>,
              MapExceptionSafetyParameters <HugeLinkedMap,  throws_on_copy,       int, true>,
              MapExceptionSafetyParameters <HugeLinkedMap,  throws_in_comparator, int, false>,
              MapExceptionSafetyParameters <HugeLinkedMap,  throws_in_comparator, int, true>,
              MapExceptionSafetyParameters <HugeLinkedMap,  throws_in_hasher,     int, false>,
              MapExceptionSafetyParameters <HugeLinkedMap,  throws_in_hasher,     int, true> >
           exception_safety_test_map_parameters_normal;

  typedef  mpl::vector
             <MapExceptionSafetyParameters <ForwardMap,     throws_on_copy,       int, false>,
              MapExceptionSafetyParameters <ForwardMap,     throws_on_copy,       int, true>,
              MapExceptionSafetyParameters <ForwardMap,     throws_in_comparator, int, false>,
              MapExceptionSafetyParameters <ForwardMap,     throws_in_comparator, int, true>,
              MapExceptionSafetyParameters <ForwardMap,     throws_in_hasher,     int, false>,
              MapExceptionSafetyParameters <ForwardMap,     throws_in_hasher,     int, true>,
              MapExceptionSafetyParameters <HugeForwardMap, throws_on_copy,       int, false>,
              MapExceptionSafetyParameters <HugeForwardMap, throws_on_copy,       int, true>,
              MapExceptionSafetyParameters <HugeForwardMap, throws_in_comparator, int, false>,
              MapExceptionSafetyParameters <HugeForwardMap, throws_in_comparator, int, true>,
              MapExceptionSafetyParameters <HugeForwardMap, throws_in_hasher,     int, false>,
              MapExceptionSafetyParameters <HugeForwardMap, throws_in_hasher,     int, true> >
           exception_safety_test_map_parameters_forward;

  typedef  mpl::insert_range <exception_safety_test_map_parameters_normal,
                              mpl::end <exception_safety_test_map_parameters_normal>::type,
                              exception_safety_test_map_parameters_forward>::type
           exception_safety_test_map_parameters;


#define COMMON_EXCEPTION_SAFETY_MAP_TEST_SETUP                          \
  typedef  typename parameters::implementation  implementation;         \
  typedef  typename parameters::tested_type     map_type;               \
  typedef  typename map_type::iterator          iterator;               \
  typedef  typename map_type::const_iterator    const_iterator;         \
  typedef  typename map_type::key_type          key_type;               \
  typedef  typename map_type::mapped_type       mapped_type;            \
  typedef  typename map_type::value_type        value_type;             \
  typedef  pair <key_type, mapped_type>         assignable_value_type;  \
  typedef  typename key_type::enabler           exception_enabler;

}


BOOST_AUTO_TEST_SUITE (map_exception_safety)


// For constructors exception safety is practically limited to 'no leaks'.
BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_range_constructor_exception_safety_1,
                               parameters, exception_safety_test_map_parameters)
{
  COMMON_EXCEPTION_SAFETY_MAP_TEST_SETUP;

  assignable_value_type  initializers[50];
  for (size_t k = 0; k < range_size (initializers); ++k)
    initializers[k] = make_pair (key_type (k), mapped_type (k));

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      {
        // The exception will probably never get thrown for large 'fail_after' values.
        exception_enabler  _(fail_after);

        try
          {
            map_type (RANGE (initializers));
          }
        catch (expected_exception&)
          { }
      }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_copy_constructor_exception_safety_1,
                               parameters, exception_safety_test_map_parameters)
{
  COMMON_EXCEPTION_SAFETY_MAP_TEST_SETUP;

  assignable_value_type  initializers[50];
  for (size_t k = 0; k < range_size (initializers); ++k)
    initializers[k] = make_pair (key_type (k), mapped_type (k));

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      // Recreate base map each iteration, because our allocator only checks for "no
      // leaks" when all containers using it are destroyed.
      map_type  map (RANGE (initializers));

      {
        // The exception will probably never get thrown for large 'fail_after' values.
        exception_enabler  _(fail_after);

        try
          {
            map_type  copy (map);
          }
        catch (expected_exception&)
          { }
      }
    }
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_initializer_constructor_exception_safety_1,
                               parameters, exception_safety_test_map_parameters)
{
  COMMON_EXCEPTION_SAFETY_MAP_TEST_SETUP;

  initializer_list <value_type>  initializers (MAP_INITIALIZERS_0_50 (key_type, mapped_type));

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      exception_enabler  _(fail_after);

      try
        {
          map_type  map (initializers);
        }
      catch (expected_exception&)
        { }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_initializer_assignment_exception_safety_1,
                               parameters, exception_safety_test_map_parameters)
{
  COMMON_EXCEPTION_SAFETY_MAP_TEST_SETUP;

  initializer_list <value_type>  initializers (MAP_INITIALIZERS_0_50  (key_type, mapped_type));
  map_type                       start_with   (MAP_INITIALIZERS_40_60 (key_type, mapped_type));
  map_type                       after_assignment (initializers);
  vector <key_type>              assigned_order   (INITIALIZERS_0_50 (key_type));

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      map_type  map (start_with);
      bool      thrown = false;

      {
        exception_enabler  _(fail_after);

        try
          {
            map = initializers;
          }
        catch (expected_exception&)
          {
            thrown = true;
          }
      }

      map.validate_integrity ();

      if (!thrown)
        {
          BOOST_REQUIRE_EQUAL (map, after_assignment);
          if (implementation::stable_order)
            assert_identical_order (map, assigned_order);
        }
      else
        {
          BOOST_REQUIRE_EQUAL (map, start_with);
          if (implementation::stable_order)
            assert_identical_order (map, keys_of (start_with));
        }
    }
}

#endif  // MCT_CXX0X_SUPPORTED


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_insert_exception_safety_1,
                               parameters, exception_safety_test_map_parameters)
{
  COMMON_EXCEPTION_SAFETY_MAP_TEST_SETUP;

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      map_type  map;
      size_t    num_inserted = 0;

      {
        exception_enabler  _(fail_after);

        try
          {
            // Limit map size since not all types will eventually throw.
            for (; num_inserted < 1000; ++num_inserted)
              map.insert (make_pair (key_type (num_inserted), mapped_type (num_inserted)));
          }
        catch (expected_exception&)
          { }
      }

      map.validate_integrity ();
      BOOST_CHECK_EQUAL (map.size (), num_inserted);
    }
}

// Here we test "reinserts", which basically results in inserting over debris buckets.
// Forward maps are excluded since we need the erase() function.
BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_insert_exception_safety_2,
                               parameters, exception_safety_test_map_parameters_normal)
{
  COMMON_EXCEPTION_SAFETY_MAP_TEST_SETUP;

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      map_type  map;
      for (int k = 0; k < 100; ++k)
        map.insert (make_pair (key_type (k), mapped_type (k)));

      bool  thrown = false;

      {
        exception_enabler  _(fail_after);

        try
          {
            for (const_iterator scan = map.begin (); scan != map.end ();)
              {
                // Don't copy the type itself, as this may throw.
                const int  value = scan->first.value;
                scan = map.erase (scan);
                implementation::insert_before (map, scan, key_type (value), mapped_type (value));

                // Properly functioning map will not rehash itself in this loop, so
                // iterators must not be invalidated midway.
                BOOST_CHECK (scan == map.end () || map.valid_iterator (scan));
              }
          }
        catch (expected_exception&)
          {
            thrown = true;
          }
      }

      map.validate_integrity ();
      BOOST_CHECK_EQUAL (map.size (), 100u - (thrown ? 1 : 0));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_insert_range_exception_safety_1,
                               parameters, exception_safety_test_map_parameters)
{
  COMMON_EXCEPTION_SAFETY_MAP_TEST_SETUP;

  assignable_value_type  initializers[50];
  for (size_t k = 0; k < range_size (initializers); ++k)
    initializers[k] = make_pair (key_type (k), mapped_type (k));

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      // Some elements are there from the start, some will be inserted.
      size_t    num_preexisted_elements = (range_size (initializers) / 2);
      map_type  map (initializers, initializers + num_preexisted_elements);
      bool      thrown = false;

      {
        exception_enabler  _(fail_after);

        try
          {
            map.insert (range_begin (initializers) + num_preexisted_elements,
                        range_end   (initializers));
          }
        catch (expected_exception&)
          {
            thrown = true;
          }
      }

      map.validate_integrity ();

      if (!thrown)
        BOOST_CHECK_EQUAL (map.size (), range_size (initializers));
      else if (!implementation::stable_order || !implementation::normal)
        {
          // Only basic exception safety, so the map needn't contain any specific number
          // of elements.
          BOOST_CHECK_GE (map.size (), num_preexisted_elements);
          BOOST_CHECK_LE (map.size (), range_size (initializers));
        }
      else
        {
          // Linked tables guarantee strong exception safety for this method.
          BOOST_CHECK_EQUAL (map.size (), num_preexisted_elements);
        }

      for (size_t k = 0; k < map.size (); ++k)
        BOOST_CHECK (contains (map, key_type (k)));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_map_sort_exception_safety_1,
                               implementation, test_map_implementations_stable_order)
{
  typedef  typename implementation::template resolve <int, int>::map_type  map_type;

  const vector <pair <int, int> >&  values        = Data <pair <int, int> >::values3 ();
  vector <pair <int, int> >         sorted_values = values;

  std::sort (RANGE (sorted_values));

  // The loop keeps going until 'fail_after' is large enough that the sorting succeeds.
  for (int fail_after = 0; ; ++fail_after)
    {
      map_type  map (RANGE (values));
      assert_identical_order (map, keys_of (values));

      throwing_int_strict_weak_ordering::enabler  _(fail_after);

      try
        {
          map.sort (throwing_int_strict_weak_ordering ());
          assert_identical_order (map, keys_of (sorted_values));
          break;
        }
      catch (expected_exception&)
        {
          map.validate_integrity ();
        }
    }
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
