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


BOOST_AUTO_TEST_SUITE (set_iterators)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_iterators_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  size_t  num_elements = 0;
  for (const_iterator iter = set.begin (); iter != set.end (); ++iter)
    ++num_elements;

  BOOST_CHECK_EQUAL (num_elements, set.size ());
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_iterators_2,
                               parameters, test_set_parameters_normal)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  size_t  num_elements = 0;
  for (const_iterator iter = set.end (); iter != set.begin (); --iter)
    ++num_elements;

  BOOST_CHECK_EQUAL (num_elements, set.size ());
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_iterator_composite_1,
                               parameters, test_set_parameters_normal)
{
  COMMON_SET_TEST_SETUP;

  set_type        set (RANGE (data::values1 ()));
  size_t          size = set.size ();
  const_iterator  scan = set.begin ();

  for (size_t k = 0; k != size; ++k, ++scan)
    {
      BOOST_CHECK (contains (set, *scan));

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

  BOOST_CHECK_EQUAL (scan, set.end ());
}


// Following tests assume MCT_CHECK_PRECONDITIONS mode.

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_invalid_iterators_1,
                               parameters, test_set_parameters_robust_iterator_validation)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  iterator  end = set.end ();
  BOOST_CHECK_THROW (++end, logic_error);

  type  value;
  BOOST_CHECK_THROW (value = *end, logic_error);
  MCT_UNUSED (value);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_invalid_iterators_2,
                               parameters, test_set_parameters_normal_robust_iterator_validation)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));

  iterator  begin = set.begin ();
  BOOST_CHECK_THROW (--begin, logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_erase_iterator_precondition_1,
                               parameters, test_set_parameters_normal)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));
  iterator  it = set.begin ();

  // This validation doesn't require fully robust checking to be available.
  set.erase (it);
  BOOST_CHECK_THROW (set.erase (it), logic_error);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (test_set_before_begin_before_end_iterators_1,
                               parameters, test_set_parameters_forward)
{
  COMMON_SET_TEST_SETUP;

  // While incrementing the before-* iterators is not particularly useful, it should still
  // work.  We test both empty and non-empty sets.

  set_type  set1;

  BOOST_CHECK_EQUAL (set1.before_begin (),   set1.before_end ());
  BOOST_CHECK_NE    (set1.before_begin (),   set1.end ());
  BOOST_CHECK_EQUAL (++set1.before_begin (), set1.begin ());
  BOOST_CHECK_EQUAL (++set1.before_end (),   set1.end ());
  BOOST_CHECK_EQUAL (++set1.before_begin (), ++set1.before_end ());

  set_type  set2 (RANGE (data::values1 ()));

  BOOST_CHECK_NE    (set2.before_begin (),   set2.before_end ());
  BOOST_CHECK_NE    (set2.before_begin (),   set2.end ());
  BOOST_CHECK_EQUAL (++set2.before_begin (), set2.begin ());
  BOOST_CHECK_EQUAL (++set2.before_end (),   set2.end ());
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
