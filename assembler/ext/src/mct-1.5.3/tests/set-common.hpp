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


#ifndef MCT_TESTS_SET_COMMON_HPP
#define MCT_TESTS_SET_COMMON_HPP


#include "tests/common.hpp"

#if HAVE_UNORDERED_SET
# include <unordered_set>
using std::unordered_set;

#elif HAVE_TR1_UNORDERED_SET
# include <tr1/unordered_set>
using std::tr1::unordered_set;

#elif HAVE_BOOST_UNORDERED_SET_HPP
# include <boost/unordered_set.hpp>
using boost::unordered_set;
#endif

#include <boost/mpl/begin.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/vector.hpp>


template <typename implementation_type>
struct robust_for <implementation_type,
                   typename enable_if <is_set <implementation_type>::value>::type>
{
  // Note that we use simple default allocator: we don't need leak testing here and neither
  // require robust type to play nicely with advanced allocators like in Boost.Interprocess,
  // for example.
  typedef  unordered_set <typename implementation_type::value_type,
                          typename implementation_type::hasher,
                          typename implementation_type::key_equal,
                          allocator <typename implementation_type::value_type> >
           type;
};


template <typename implementation_type>
class test_container_interface <implementation_type,
                                typename enable_if <is_set <implementation_type>
                                                    ::value>::type>
  : public test_family_interface <implementation_type>
{
  typedef  test_family_interface <implementation_type>  base_type;

public:

  typedef  typename base_type::key_type        key_type;
  typedef  typename base_type::const_iterator  const_iterator;


protected:

  virtual  const key_type&
  extract_key (const_iterator iterator) const
  {  return *iterator;  }
};


template <typename tested_set_type, typename allocator_type>
bool
operator== (const tested_set_type& tested_set,
            const unordered_set <typename tested_set_type::value_type,
                                 typename tested_set_type::hasher,
                                 typename tested_set_type::key_equal,
                                 allocator_type>& robust_set)
{
  if (tested_set.size () != robust_set.size ())
    return false;

  typedef  unordered_set <typename tested_set_type::value_type,
                          typename tested_set_type::hasher,
                          typename tested_set_type::key_equal,
                          allocator_type>
           robust_set_type;

  for (typename robust_set_type::const_iterator scan = robust_set.begin ();
       scan != robust_set.end (); ++scan)
    {
      if (!contains (tested_set, *scan))
        return false;
    }

  return true;
}


namespace std
{

  template <typename set_type>
  typename enable_if <impl::is_same <typename set_type::key_type,
                                     typename set_type::value_type>::value, ostream&>::type
  operator<< (ostream& stream, const set_type& set)
  {
    return print_range (stream, "{}", set.begin (), set.end ());
  }

}


struct PlainSet
{
  template <typename Type,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Type>,
            typename Equal     = equal_to <Type>,
            typename Allocator = allocator <Type>,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  closed_hash_set <Type, Hash, Equal,
                              leak_test_allocator <Allocator>,
                              keep_hashes>
             set_type;
  };

  static  const bool  normal       = true;
  static  const bool  stable_order = false;

  // Hide difference in available methods: for non-linked sets 'before' is just ignored.
  template <typename Set>
  static  pair <typename Set::iterator, bool>
  insert_before (Set& set, typename Set::const_iterator /* before */,
                 const typename Set::value_type& value)
  {
    return set.insert (value);
  }
};


struct LinkedSet
{
  template <typename Type,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Type>,
            typename Equal     = equal_to <Type>,
            typename Allocator = allocator <Type>,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  linked_hash_set <Type, Hash, Equal,
                              leak_test_allocator <Allocator>,
                              keep_hashes>
             set_type;
  };

  static  const bool  normal       = true;
  static  const bool  stable_order = true;

  template <typename Set>
  static  pair <typename Set::iterator, bool>
  insert_before (Set& set, typename Set::const_iterator before,
                 const typename Set::value_type& value)
  {
    return set.insert (before, value);
  }
};


struct HugeLinkedSet : LinkedSet
{
  template <typename Type,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Type>,
            typename Equal     = equal_to <Type>,
            typename Allocator = allocator <Type>,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  huge_linked_hash_set <Type, Hash, Equal,
                                   leak_test_allocator <Allocator>,
                                   keep_hashes>
             set_type;
  };
};


struct ForwardSet
{
  template <typename Type,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Type>,
            typename Equal     = equal_to <Type>,
            typename Allocator = allocator <Type>,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  forward_hash_set <Type, Hash, Equal,
                               leak_test_allocator <Allocator>,
                               keep_hashes>
             set_type;
  };

  static  const bool  normal       = false;
  static  const bool  stable_order = true;

  template <typename Set>
  static  pair <typename Set::iterator, bool>
  insert_before (Set& set, typename Set::const_iterator before,
                 const typename Set::value_type& value)
  {
    return set.insert (before, value);
  }
};


struct HugeForwardSet : ForwardSet
{
  template <typename Type,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Type>,
            typename Equal     = equal_to <Type>,
            typename Allocator = allocator <Type>,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  huge_forward_hash_set <Type, Hash, Equal,
                                    leak_test_allocator <Allocator>,
                                    keep_hashes>
             set_type;
  };
};


template <typename Implementation,
          typename Type,
          typename Hash      = MCT_HASH_NAMESPACE::hash <Type>,
          typename Equal     = equal_to <Type>,
          typename Allocator = allocator <Type>,
          bool _keep_hashes  = false>
struct SetParameters
{
  typedef  Implementation  implementation;
  typedef  Type            type;
  typedef  Hash            hasher;
  typedef  Equal           equal_to;

  static  const bool  keep_hashes = _keep_hashes;

  typedef  typename Implementation::
           template resolve <Type, Hash, Equal, Allocator, _keep_hashes>::set_type
           tested_type;
  typedef  test_table <tested_type>  wrapper_type;
};


// To speed compilation up we list just a few possible combinations for each base type.

// 'std_equal' need to never compare different elements in 'data::valueN' equal and not
// compare different elements returned from data::generate() equal, at least for
// reasonably small iteration numbers.
typedef  mpl::vector <SetParameters <PlainSet, int>,
                      SetParameters <PlainSet, int, int_hasher_1>,
                      SetParameters <PlainSet, int, MCT_HASH_NAMESPACE::hash <int>,
                                     equal_to <int>, allocator <int>, true>,
                      SetParameters <LinkedSet, int>,
                      SetParameters <HugeLinkedSet, int>,
                      SetParameters <ForwardSet, int>,
                      SetParameters <HugeForwardSet, int>,
#                   if MCT_HAVE_LONG_LONG
                      SetParameters <PlainSet, unsigned long long>,
                      SetParameters <PlainSet, unsigned long long, unsigned_long_long_hasher>,
                      SetParameters <LinkedSet, unsigned long long, unsigned_long_long_hasher,
                                     equal_to <unsigned long long>, allocator <unsigned long long>,
                                     true>,
#                   endif
                      SetParameters <PlainSet, string>,
                      SetParameters <PlainSet, string, string_hasher_1>,
                      SetParameters <PlainSet, string, string_hasher_1, equal_to <string>,
                                     allocator <string>, true>,
                      SetParameters <LinkedSet, string, string_hasher_1>,
                      SetParameters <PlainSet, foo>,
                      SetParameters <LinkedSet, foo>,
                      SetParameters <PlainSet, bar>,
                      SetParameters <LinkedSet, bar>,
                      SetParameters <PlainSet, ham> >
         test_set_parameters_std_equal_std_allocator;

typedef  mpl::vector <
#                   if HAVE_BOOST_INTERPROCESS
                      SetParameters <PlainSet, int,
                                     MCT_HASH_NAMESPACE::hash <int>, equal_to <int>,
                                     interprocess_allocator <int> >,
                      SetParameters <LinkedSet, foo,
                                     MCT_HASH_NAMESPACE::hash <foo>, equal_to <foo>,
                                     interprocess_allocator <foo> >,
                      SetParameters <HugeLinkedSet, bar,
                                     MCT_HASH_NAMESPACE::hash <bar>, equal_to <bar>,
                                     interprocess_allocator <bar> >,
                      SetParameters <ForwardSet, bar,
                                     MCT_HASH_NAMESPACE::hash <bar>, equal_to <bar>,
                                     interprocess_allocator <bar> >,
                      SetParameters <HugeForwardSet, ham,
                                     MCT_HASH_NAMESPACE::hash <ham>, equal_to <ham>,
                                     interprocess_allocator <ham> >
#                   endif
                     >
         test_set_parameters_std_equal_interp_allocator;

typedef  mpl::insert_range <test_set_parameters_std_equal_std_allocator,
                            mpl::end <test_set_parameters_std_equal_std_allocator>::type,
                            test_set_parameters_std_equal_interp_allocator>::type
         test_set_parameters_std_equal;

typedef  mpl::push_back <test_set_parameters_std_equal,
                         SetParameters <PlainSet, int, int_hasher_special,
                                        int_comparator_special> >::type
         test_set_parameters;

typedef  mpl::copy_if <test_set_parameters_std_equal, normal_parameters <mpl::_1> >::type
         test_set_parameters_std_equal_normal;

typedef  mpl::copy_if <test_set_parameters_std_equal, stable_order_parameters <mpl::_1> >::type
         test_set_parameters_std_equal_stable_order;

typedef  mpl::copy_if <test_set_parameters_std_equal, linked_parameters <mpl::_1> >::type
         test_set_parameters_std_equal_linked;

typedef  mpl::copy_if <test_set_parameters, normal_parameters <mpl::_1> >::type
         test_set_parameters_normal;

typedef  mpl::copy_if <test_set_parameters, stable_order_parameters <mpl::_1> >::type
         test_set_parameters_stable_order;

typedef  mpl::copy_if <test_set_parameters, linked_parameters <mpl::_1> >::type
         test_set_parameters_linked;

typedef  mpl::copy_if <test_set_parameters, forward_parameters <mpl::_1> >::type
         test_set_parameters_forward;

typedef  mpl::copy_if <test_set_parameters, robust_iterator_validation <mpl::_1> >::type
         test_set_parameters_robust_iterator_validation;

typedef  mpl::copy_if <test_set_parameters_normal, robust_iterator_validation <mpl::_1> >::type
         test_set_parameters_normal_robust_iterator_validation;

typedef  mpl::vector <LinkedSet, HugeLinkedSet, ForwardSet, HugeForwardSet>
         test_set_implementations_stable_order;

typedef  mpl::insert <test_set_implementations_stable_order,
                      mpl::begin <test_set_implementations_stable_order>::type, PlainSet>::type
         test_set_implementations;


#define COMMON_SET_TEST_SETUP                                           \
  typedef  typename parameters::implementation     implementation;      \
  typedef  typename parameters::wrapper_type       set_type;            \
  typedef  typename set_type::implementation_type  set_base_type;       \
  typedef  typename set_type::iterator             iterator;            \
  typedef  typename set_type::const_iterator       const_iterator;      \
  typedef  typename parameters::type               type;                \
  typedef  typename parameters::hasher             hasher;              \
  typedef  typename parameters::equal_to           equal_to;            \
  typedef  Data <type>                             data;


#if HAVE_BOOST_INTERPROCESS

// We use these magic numbers to avoid having to combine all the various (sets, maps) type
// into the same MPL collection.
# define TEST_SET_TYPE_MAGIC_NUMBER     0x10000

typedef  mpl::transform <test_set_parameters_std_equal_interp_allocator,
                         get_tested <mpl::_1> >::type
         test_set_types_interp_allocator;

template <typename type>
struct external_validator
  <type,
   typename enable_if <is_set <type>::value
                       && is_interprocess_allocator <typename type::allocator_type>::value>
   ::type>
{
  static  void
  validate (const type& object)
  {
    typedef  typename mpl::find <test_set_types_interp_allocator, type>::type  position;
    typedef  mpl::begin <test_set_types_interp_allocator>::type                begin;

    validate_externally (&object, (TEST_SET_TYPE_MAGIC_NUMBER
                                   + mpl::distance <begin, position>::type::value));
  }
};

#endif  // HAVE_BOOST_INTERPROCESS


#endif  // Multi-inclusion guard.


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
