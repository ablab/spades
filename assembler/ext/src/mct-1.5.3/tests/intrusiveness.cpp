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


#include "tests/common.hpp"

#include <boost/mpl/vector.hpp>


namespace
{

  // Will define external use rules for this.
  struct ham_wrapper_1
  {
    ham  the_ham;
  };

  // But not for this.
  struct ham_wrapper_2
  {
    ham  the_ham;
  };

  enum enum_type
    {
      DUMMY_VALUE
    };

  union union_type
  {
    int  x;
  };

  struct class_type
  { };


  typedef  mpl::vector <bool, int, char, void*, int[42], enum_type, union_type,
                        void (*) (), void (class_type::*) ()>
           non_subclassable_default_constructible_types;

}


namespace mct
{
  template <>
  struct external_use <ham_wrapper_1>
    : extern_use_field <ham_wrapper_1, ham, &ham_wrapper_1::the_ham>
  { };
}


BOOST_AUTO_TEST_SUITE (intrusiveness)


BOOST_AUTO_TEST_CASE (test_simple_1)
{
  BOOST_CHECK (!supports_external_use <bool> ::value);
  BOOST_CHECK (!supports_external_use <char> ::value);
  BOOST_CHECK (!supports_external_use <int>  ::value);
  BOOST_CHECK (!supports_external_use <float>::value);

  // Type is defined in 'common.hpp'.
  BOOST_CHECK (supports_external_use <ham>::value);
}


BOOST_AUTO_TEST_CASE (test_implicit_propagation_1)
{
  // Need for commas in macros is always a pain in the ass.
  typedef  pair <bool, char>   pair_bool_char;
  typedef  pair <char, int>    pair_char_int;
  typedef  pair <int,  float>  pair_int_float;
  typedef  pair <int,  int>    pair_int_int;
  typedef  pair <ham,  bool>   pair_ham_bool;
  typedef  pair <bool, ham>    pair_bool_ham;
  typedef  pair <ham,  ham>    pair_ham_ham;

  BOOST_CHECK (!supports_external_use <pair_bool_char>::value);
  BOOST_CHECK (!supports_external_use <pair_char_int> ::value);
  BOOST_CHECK (!supports_external_use <pair_int_float>::value);
  BOOST_CHECK (!supports_external_use <pair_int_int>  ::value);

  BOOST_CHECK (supports_external_use <pair_ham_bool>::value);
  BOOST_CHECK (supports_external_use <pair_bool_ham>::value);
  BOOST_CHECK (supports_external_use <pair_ham_ham> ::value);
}


BOOST_AUTO_TEST_CASE (test_explicit_propagation_1)
{
  BOOST_CHECK (supports_external_use <ham_wrapper_1>::value);

  static  const bool  expected_type = impl::is_same <external_use <ham_wrapper_1>::value_type,
                                                     external_use <ham>::value_type>::value;
  BOOST_CHECK (expected_type);

  BOOST_CHECK (!supports_external_use <ham_wrapper_2>::value);
}


BOOST_AUTO_TEST_CASE (test_intrusive_storage_subclassable_1)
{
  struct size_check : pair <int, char>
  {
    char  extra;
  };

  static  const bool  is_same
    = impl::is_same <intrusive_storage <std::pair <int, char> >::type,
                     std::pair <int, char> >::value;

  BOOST_CHECK_EQUAL (is_same, sizeof (pair <int, char>) != sizeof (size_check));
}

BOOST_AUTO_TEST_CASE_TEMPLATE (test_intrusive_storage_non_subclassable_1,
                               type, non_subclassable_default_constructible_types)
{
  // This tests compilation more than the result, i.e. that the intrusiveness code
  // actually doesn't try to subclass such types.
  static  const bool  is_same
    = impl::is_same <typename intrusive_storage <type>::type, type>::value;
  BOOST_CHECK (is_same);
}


BOOST_AUTO_TEST_SUITE_END ()


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
