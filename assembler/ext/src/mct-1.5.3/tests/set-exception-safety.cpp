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


#include "tests/set-common.hpp"

#include <boost/mpl/end.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/vector.hpp>


namespace
{

  template <typename Implementation, typename Type, bool _keep_hashes = false>
  struct SetExceptionSafetyParameters
  {
    typedef  Implementation  implementation;
    typedef  Type            type;

    static  const bool  keep_hashes  = _keep_hashes;
    static  const bool  stable_order = Implementation::stable_order;
    static  const bool  normal       = Implementation::normal;

    // Not 'test_set' here as we can't control where exception are actually thrown.
    typedef  typename Implementation::template resolve <Type,
                                                        MCT_HASH_NAMESPACE::hash <type>,
                                                        equal_to <type>,
                                                        allocator <type>,
                                                        keep_hashes>::set_type
             tested_type;
  };


  typedef  mpl::vector
             <SetExceptionSafetyParameters <PlainSet,       throws_on_copy,       false>,
              SetExceptionSafetyParameters <PlainSet,       throws_on_copy,       true>,
              SetExceptionSafetyParameters <PlainSet,       throws_in_comparator, false>,
              SetExceptionSafetyParameters <PlainSet,       throws_in_comparator, true>,
              SetExceptionSafetyParameters <PlainSet,       throws_in_hasher,     false>,
              SetExceptionSafetyParameters <PlainSet,       throws_in_hasher,     true>,
              SetExceptionSafetyParameters <LinkedSet,      throws_on_copy,       false>,
              SetExceptionSafetyParameters <LinkedSet,      throws_on_copy,       true>,
              SetExceptionSafetyParameters <LinkedSet,      throws_in_comparator, false>,
              SetExceptionSafetyParameters <LinkedSet,      throws_in_comparator, true>,
              SetExceptionSafetyParameters <LinkedSet,      throws_in_hasher,     false>,
              SetExceptionSafetyParameters <LinkedSet,      throws_in_hasher,     true>,
              SetExceptionSafetyParameters <HugeLinkedSet,  throws_on_copy,       false>,
              SetExceptionSafetyParameters <HugeLinkedSet,  throws_on_copy,       true>,
              SetExceptionSafetyParameters <HugeLinkedSet,  throws_in_comparator, false>,
              SetExceptionSafetyParameters <HugeLinkedSet,  throws_in_comparator, true>,
              SetExceptionSafetyParameters <HugeLinkedSet,  throws_in_hasher,     false>,
              SetExceptionSafetyParameters <HugeLinkedSet,  throws_in_hasher,     true> >
           exception_safety_test_set_parameters_normal;

  typedef  mpl::vector
             <SetExceptionSafetyParameters <ForwardSet,     throws_on_copy,       false>,
              SetExceptionSafetyParameters <ForwardSet,     throws_on_copy,       true>,
              SetExceptionSafetyParameters <ForwardSet,     throws_in_comparator, false>,
              SetExceptionSafetyParameters <ForwardSet,     throws_in_comparator, true>,
              SetExceptionSafetyParameters <ForwardSet,     throws_in_hasher,     false>,
              SetExceptionSafetyParameters <ForwardSet,     throws_in_hasher,     true>,
              SetExceptionSafetyParameters <HugeForwardSet, throws_on_copy,       false>,
              SetExceptionSafetyParameters <HugeForwardSet, throws_on_copy,       true>,
              SetExceptionSafetyParameters <HugeForwardSet, throws_in_comparator, false>,
              SetExceptionSafetyParameters <HugeForwardSet, throws_in_comparator, true>,
              SetExceptionSafetyParameters <HugeForwardSet, throws_in_hasher,     false>,
              SetExceptionSafetyParameters <HugeForwardSet, throws_in_hasher,     true> >
           exception_safety_test_set_parameters_forward;

  typedef  mpl::insert_range <exception_safety_test_set_parameters_normal,
                              mpl::end <exception_safety_test_set_parameters_normal>::type,
                              exception_safety_test_set_parameters_forward>::type
           exception_safety_test_set_parameters;


#define COMMON_EXCEPTION_SAFETY_SET_TEST_SETUP                          \
  typedef  typename parameters::implementation  implementation;         \
  typedef  typename parameters::tested_type     set_type;               \
  typedef  typename set_type::iterator          iterator;               \
  typedef  typename set_type::const_iterator    const_iterator;         \
  typedef  typename set_type::value_type        type;                   \
  typedef  typename type::enabler               exception_enabler;

}


BOOST_AUTO_TEST_SUITE (set_exception_safety)


// For constructors exception safety is practically limited to 'no leaks'.
BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_range_constructor_exception_safety_1,
                               parameters, exception_safety_test_set_parameters)
{
  COMMON_EXCEPTION_SAFETY_SET_TEST_SETUP;

  type  initializers[50];
  for (size_t k = 0; k < range_size (initializers); ++k)
    initializers[k] = type (k);

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      {
        // The exception will probably never get thrown for large 'fail_after' values.
        exception_enabler  _(fail_after);

        try
          {
            set_type (RANGE (initializers));
          }
        catch (expected_exception&)
          { }
      }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_copy_constructor_exception_safety_1,
                               parameters, exception_safety_test_set_parameters)
{
  COMMON_EXCEPTION_SAFETY_SET_TEST_SETUP;

  type  initializers[50];
  for (size_t k = 0; k < range_size (initializers); ++k)
    initializers[k] = type (k);

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      // Recreate base set each iteration, because our allocator only checks for "no
      // leaks" when all containers using it are destroyed.
      set_type  set (RANGE (initializers));

      {
        // The exception will probably never get thrown for large 'fail_after' values.
        exception_enabler  _(fail_after);

        try
          {
            set_type  copy (set);
          }
        catch (expected_exception&)
          { }
      }
    }
}


#if MCT_CXX0X_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_initializer_constructor_exception_safety_1,
                               parameters, exception_safety_test_set_parameters)
{
  COMMON_EXCEPTION_SAFETY_SET_TEST_SETUP;

  initializer_list <type>  initializers (INITIALIZERS_0_50 (type));

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      exception_enabler  _(fail_after);

      try
        {
          set_type  set (initializers);
        }
      catch (expected_exception&)
        { }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_initializer_assignment_exception_safety_1,
                               parameters, exception_safety_test_set_parameters)
{
  COMMON_EXCEPTION_SAFETY_SET_TEST_SETUP;

  initializer_list <type>  initializers (INITIALIZERS_0_50  (type));
  set_type                 start_with   (INITIALIZERS_40_60 (type));
  set_type                 after_assignment (initializers);
  vector <type>            assigned_order   (initializers);

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      set_type  set (start_with);
      bool      thrown = false;

      {
        exception_enabler  _(fail_after);

        try
          {
            set = initializers;
          }
        catch (expected_exception&)
          {
            thrown = true;
          }
      }

      set.validate_integrity ();

      if (!thrown)
        {
          BOOST_REQUIRE_EQUAL (set, after_assignment);
          if (implementation::stable_order)
            assert_identical_order (set, assigned_order);
        }
      else
        {
          BOOST_REQUIRE_EQUAL (set, start_with);
          if (implementation::stable_order)
            assert_identical_order (set, start_with);
        }
    }
}

#endif  // MCT_CXX0X_SUPPORTED


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_insert_exception_safety_1,
                               parameters, exception_safety_test_set_parameters)
{
  COMMON_EXCEPTION_SAFETY_SET_TEST_SETUP;

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      set_type  set;
      size_t    num_inserted = 0;

      {
        exception_enabler  _(fail_after);

        try
          {
            // Limit set size since not all types will eventually throw.
            for (; num_inserted < 1000; ++num_inserted)
              set.insert (type (num_inserted));
          }
        catch (expected_exception&)
          { }
      }

      set.validate_integrity ();
      BOOST_CHECK_EQUAL (set.size (), num_inserted);
    }
}

// Here we test "reinserts", which basically results in inserting over debris buckets.
// Forward sets are excluded since we need the erase() function.
BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_insert_exception_safety_2,
                               parameters, exception_safety_test_set_parameters_normal)
{
  COMMON_EXCEPTION_SAFETY_SET_TEST_SETUP;

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      set_type  set;
      for (int k = 0; k < 100; ++k)
        set.insert (type (k));

      bool  thrown = false;

      {
        exception_enabler  _(fail_after);

        try
          {
            for (const_iterator scan = set.begin (); scan != set.end ();)
              {
                // Don't copy the type itself, as this may throw.
                const int  value = scan->value;
                scan = set.erase (scan);
                implementation::insert_before (set, scan, type (value));

                // Properly functioning set will not rehash itself in this loop, so
                // iterators must not be invalidated midway.
                BOOST_CHECK (scan == set.end () || set.valid_iterator (scan));
              }
          }
        catch (expected_exception&)
          {
            thrown = true;
          }
      }

      set.validate_integrity ();
      BOOST_CHECK_EQUAL (set.size (), 100u - (thrown ? 1 : 0));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_insert_range_exception_safety_1,
                               parameters, exception_safety_test_set_parameters)
{
  COMMON_EXCEPTION_SAFETY_SET_TEST_SETUP;

  type  initializers[100];
  for (size_t k = 0; k < range_size (initializers); ++k)
    initializers[k] = type (k);

  for (int fail_after = 0; fail_after < 100; ++fail_after)
    {
      BOOST_TEST_CHECKPOINT ("fail after " << fail_after << " operations");

      // Some elements are there from the start, some will be inserted.
      size_t    num_preexisted_elements = (range_size (initializers) / 2);
      set_type  set (initializers, initializers + num_preexisted_elements);
      bool      thrown = false;

      {
        exception_enabler  _(fail_after);

        try
          {
            set.insert (range_begin (initializers) + num_preexisted_elements,
                        range_end   (initializers));
          }
        catch (expected_exception&)
          {
            thrown = true;
          }
      }

      set.validate_integrity ();

      if (!thrown)
        BOOST_CHECK_EQUAL (set.size (), range_size (initializers));
      else if (!implementation::stable_order || !implementation::normal)
        {
          // Only basic exception safety, so the set needn't contain any specific number
          // of elements.
          BOOST_CHECK_GE (set.size (), num_preexisted_elements);
          BOOST_CHECK_LE (set.size (), range_size (initializers));
        }
      else
        {
          // Linked tables guarantee strong exception safety for this method.
          BOOST_CHECK_EQUAL (set.size (), num_preexisted_elements);
        }

      for (size_t k = 0; k < set.size (); ++k)
        BOOST_CHECK (contains (set, type (k)));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_sort_exception_safety_1,
                               implementation, test_set_implementations_stable_order)
{
  typedef  typename implementation::template resolve <int>::set_type  set_type;

  const vector <int>&  values        = Data <int>::values3 ();
  vector <int>         sorted_values = values;

  std::sort (RANGE (sorted_values));

  // The loop keeps going until 'fail_after' is large enough that the sorting succeeds.
  for (int fail_after = 0; ; ++fail_after)
    {
      set_type  set (RANGE (values));
      assert_identical_order (set, values);

      throwing_int_strict_weak_ordering::enabler  _(fail_after);

      try
        {
          set.sort (throwing_int_strict_weak_ordering ());
          assert_identical_order (set, sorted_values);
          break;
        }
      catch (expected_exception&)
        {
          set.validate_integrity ();
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
