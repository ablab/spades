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

#include "tests/set-common.hpp"


BOOST_AUTO_TEST_SUITE (set_empty)


BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_set_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set;

  BOOST_CHECK       (set.empty ());
  BOOST_CHECK_EQUAL (set.size (), 0u);
  BOOST_CHECK       (!set.uses_dynamic_memory ());
  BOOST_CHECK_EQUAL (set.used_memory (), sizeof (set_base_type));
  BOOST_CHECK_EQUAL (set.begin  (), set.end  ());
  BOOST_CHECK_EQUAL (set.cbegin (), set.cend ());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_set_explicit_size_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (1000);

  BOOST_CHECK (!set.uses_dynamic_memory ());

  set.insert (data::generate (0));

  BOOST_CHECK    (set.uses_dynamic_memory ());
  BOOST_CHECK_GE (set.bucket_count (), 1000);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_set_copy_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set;
  set_type  set2 (set);

  BOOST_CHECK       (set2.empty ());
  BOOST_CHECK_EQUAL (set2.size (), 0u);
  BOOST_CHECK       (!set2.uses_dynamic_memory ());
  BOOST_CHECK_EQUAL (set2.begin  (), set2.end  ());
  BOOST_CHECK_EQUAL (set2.cbegin (), set2.cend ());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_set_assignment_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set;
  set_type  set2;

  set = set2;

  BOOST_CHECK       (set.empty ());
  BOOST_CHECK_EQUAL (set.size (), 0u);
  BOOST_CHECK       (!set.uses_dynamic_memory ());
  BOOST_CHECK_EQUAL (set.begin  (), set.end  ());
  BOOST_CHECK_EQUAL (set.cbegin (), set.cend ());

  set.insert (RANGE (data::values1 ()));
  set = set2;

  // Check the same invariants if assigning over a non-empty set.
  BOOST_CHECK       (set.empty ());
  BOOST_CHECK_EQUAL (set.size (), 0u);
  BOOST_CHECK       (!set.uses_dynamic_memory ());
  BOOST_CHECK_EQUAL (set.begin  (), set.end  ());
  BOOST_CHECK_EQUAL (set.cbegin (), set.cend ());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_set_erase_1,
                               parameters, test_set_parameters_normal)
{
  COMMON_SET_TEST_SETUP;

  set_type  set;

  for (typename vector <type>::const_iterator scan = data::values1 ().begin ();
       scan != data::values1 ().end (); ++scan)
    BOOST_CHECK_EQUAL (set.erase (*scan), 0u);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_set_rehash_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set;

  BOOST_CHECK (!set.uses_dynamic_memory ());

  set.rehash (1000);

  BOOST_CHECK    (!set.uses_dynamic_memory ());
  BOOST_CHECK_GE (set.bucket_count (), 1000);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_empty_set_rehash_cleared_1,
                               parameters, test_set_parameters)
{
  COMMON_SET_TEST_SETUP;

  set_type  set (RANGE (data::values1 ()));
  set.clear ();

  // Hash tables mustn't just drop memory once they become empty, as that would lead to
  // terrible performance in cases where becoming empty and then being filled again is the
  // norm.
  BOOST_CHECK (set.uses_dynamic_memory ());

  set.rehash (2 * set.bucket_count ());

  // However, when asked to resize internal table anyway, there is no point in allocating
  // new array right away if the table is empty.
  BOOST_CHECK (!set.uses_dynamic_memory ());
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
