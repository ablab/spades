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


#ifndef MCT_TESTS_MAP_COMMON_HPP
#define MCT_TESTS_MAP_COMMON_HPP


#include "tests/common.hpp"

#if HAVE_UNORDERED_MAP
# include <unordered_map>
using std::unordered_map;

#elif HAVE_TR1_UNORDERED_MAP
# include <tr1/unordered_map>
using std::tr1::unordered_map;

#elif HAVE_BOOST_UNORDERED_MAP_HPP
# include <boost/unordered_map.hpp>
using boost::unordered_map;
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
                   typename enable_if <is_map <implementation_type>::value>::type>
{
  // Note that we use simple default allocator: we don't need leak testing here and neither
  // require robust type to play nicely with advanced allocators like in Boost.Interprocess,
  // for example.
  typedef  unordered_map <typename implementation_type::key_type,
                          typename implementation_type::mapped_type,
                          typename implementation_type::hasher,
                          typename implementation_type::key_equal,
                          allocator <typename implementation_type::value_type> >
           type;
};


template <typename implementation_type>
class test_container_interface <implementation_type,
                                typename enable_if <is_map <implementation_type>
                                                    ::value>::type>
  : public test_family_interface <implementation_type>
{
  typedef  test_family_interface <implementation_type>  base_type;

public:

  typedef  typename base_type::key_type               key_type;
  typedef  typename implementation_type::mapped_type  mapped_type;
  typedef  typename base_type::iterator               iterator;
  typedef  typename base_type::const_iterator         const_iterator;


  mapped_type&
  at (const key_type& key)
  {  return this->impl ().at (key);  }

  const mapped_type&
  at (const key_type& key) const
  {  return this->impl ().at (key);  }


  void
  assign_at_key (const key_type& key, const mapped_type& mapped)
  {
    this->impl  () [key] = mapped;
    this->valid () [key] = mapped;
    this->validate ();
  }

  // Assigning directly to iterator is of course possible, but will break all further
  // validate() calls.  Subclassing iterators sounds a bit too difficult to me.
  void
  assign_through_iterator (iterator position, mapped_type value)
  {
    BOOST_REQUIRE_EQUAL (position, this->find (position->first));
    position->second                              = value;
    this->valid ().find (position->first)->second = value;
    this->validate ();
  }


protected:

  virtual  const key_type&
  extract_key (const_iterator iterator) const
  {  return iterator->first;  }
};


template <typename tested_map_type, typename allocator_type>
bool
operator== (const tested_map_type& tested_map,
            const unordered_map <typename tested_map_type::key_type,
                                 typename tested_map_type::mapped_type,
                                 typename tested_map_type::hasher,
                                 typename tested_map_type::key_equal,
                                 allocator_type>& robust_map)
{
  if (tested_map.size () != robust_map.size ())
    return false;

  typedef  unordered_map <typename tested_map_type::key_type,
                          typename tested_map_type::mapped_type,
                          typename tested_map_type::hasher,
                          typename tested_map_type::key_equal,
                          allocator_type>
           robust_map_type;

  for (typename robust_map_type::const_iterator scan = robust_map.begin ();
       scan != robust_map.end (); ++scan)
    {
      typename tested_map_type::const_iterator  match = tested_map.find (scan->first);
      if (match == tested_map.end () || !(match->second == scan->second))
        return false;
    }

  return true;
}


namespace std
{

  template <typename map_type>
  typename enable_if <impl::is_same <typename map_type::value_type,
                                     pair <const typename map_type::key_type,
                                           typename map_type::mapped_type> >::value, ostream&>
           ::type
  operator<< (ostream& stream, const map_type& map)
  {
    return print_map_range (stream, "{}", map.begin (), map.end ());
  }

}


struct PlainMap
{
  template <typename Key, typename Mapped,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal     = equal_to <Key>,
            typename Allocator = allocator <pair <const Key, Mapped> >,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  closed_hash_map <Key, Mapped, Hash, Equal, leak_test_allocator <Allocator>,
                              keep_hashes>
             map_type;
  };

  static  const bool  normal       = true;
  static  const bool  stable_order = false;

  // Hide difference in available methods: for non-linked maps 'before' is just ignored.
  template <typename Map>
  static  pair <typename Map::iterator, bool>
  insert_before (Map& map, typename Map::const_iterator /* before */,
                 const typename Map::key_type& key, const typename Map::mapped_type& mapped)
  {
    return map.insert (make_pair (key, mapped));
  }
};


struct LinkedMap
{
  template <typename Key, typename Mapped,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal     = equal_to <Key>,
            typename Allocator = allocator <pair <const Key, Mapped> >,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  linked_hash_map <Key, Mapped, Hash, Equal, leak_test_allocator <Allocator>,
                              keep_hashes>
             map_type;
  };

  static  const bool  normal       = true;
  static  const bool  stable_order = true;

  // Hide difference in available methods: for non-linked maps 'before' is just ignored.
  template <typename Map>
  static  pair <typename Map::iterator, bool>
  insert_before (Map& map, typename Map::const_iterator before,
                 const typename Map::key_type& key, const typename Map::mapped_type& mapped)
  {
    return map.insert (before, make_pair (key, mapped));
  }
};


struct HugeLinkedMap : LinkedMap
{
  template <typename Key, typename Mapped,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal     = equal_to <Key>,
            typename Allocator = allocator <pair <const Key, Mapped> >,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  huge_linked_hash_map <Key, Mapped, Hash, Equal, leak_test_allocator <Allocator>,
                                   keep_hashes>
             map_type;
  };
};


struct ForwardMap
{
  template <typename Key, typename Mapped,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal     = equal_to <Key>,
            typename Allocator = allocator <pair <const Key, Mapped> >,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  forward_hash_map <Key, Mapped, Hash, Equal, leak_test_allocator <Allocator>,
                               keep_hashes>
             map_type;
  };

  static  const bool  normal       = false;
  static  const bool  stable_order = true;

  // Hide difference in available methods: for non-linked maps 'before' is just ignored.
  template <typename Map>
  static  pair <typename Map::iterator, bool>
  insert_before (Map& map, typename Map::const_iterator before,
                 const typename Map::key_type& key, const typename Map::mapped_type& mapped)
  {
    return map.insert (before, make_pair (key, mapped));
  }
};


struct HugeForwardMap : ForwardMap
{
  template <typename Key, typename Mapped,
            typename Hash      = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal     = equal_to <Key>,
            typename Allocator = allocator <pair <const Key, Mapped> >,
            bool keep_hashes   = false>
  struct resolve
  {
    typedef  huge_forward_hash_map <Key, Mapped, Hash, Equal, leak_test_allocator <Allocator>,
                                    keep_hashes>
             map_type;
  };
};


template <typename Implementation,
          typename Key,
          typename Mapped,
          typename Hash      = MCT_HASH_NAMESPACE::hash <Key>,
          typename Equal     = equal_to <Key>,
          typename Allocator = allocator <pair <const Key, Mapped> >,
          bool _keep_hashes  = false>
struct MapParameters
{
  typedef  Implementation                      implementation;
  typedef  Key                                 key_type;
  typedef  Mapped                              mapped_type;
  typedef  pair <const key_type, mapped_type>  value_type;
  typedef  pair <key_type, mapped_type>        assignable_value_type;
  typedef  Hash                                hasher;
  typedef  Equal                               equal_to;

  static  const bool  keep_hashes = _keep_hashes;

  typedef  typename Implementation::
           template resolve <Key, Mapped, Hash, Equal, Allocator, _keep_hashes>::map_type
           tested_type;
  typedef  test_table <tested_type>  wrapper_type;
};


// To speed compilation up we list just a few possible combinations for each base type.
typedef  mpl::vector <MapParameters <PlainMap, int, int>,
                      MapParameters <LinkedMap, int, int>,
                      MapParameters <HugeLinkedMap, int, int>,
                      MapParameters <ForwardMap, int, int>,
                      MapParameters <HugeForwardMap, int, int>,
                      MapParameters <PlainMap, int, string, int_hasher_1>,
                      MapParameters <PlainMap, int, string,
                                     MCT_HASH_NAMESPACE::hash <int>, equal_to <int>,
                                     allocator <pair <const int, string> >, true>,
                      MapParameters <PlainMap, int, ham>,
#                   if MCT_HAVE_LONG_LONG
                      MapParameters <PlainMap, unsigned long long, string>,
                      MapParameters <PlainMap, unsigned long long, int, unsigned_long_long_hasher>,
                      MapParameters <LinkedMap, unsigned long long, string,
                                     MCT_HASH_NAMESPACE::hash <unsigned long long>,
                                     equal_to <unsigned long long>,
                                     allocator <pair <const unsigned long long, string> >, true>,
#                   endif
                      MapParameters <PlainMap, string, int>,
                      MapParameters <PlainMap, string, string, string_hasher_1>,
                      MapParameters <PlainMap, string, int, string_hasher_1, equal_to <string>,
                                     allocator <pair <const string, int> >, true>,
                      MapParameters <LinkedMap, string, string>,
                      MapParameters <PlainMap, foo, int>,
                      MapParameters <LinkedMap, foo, string>,
                      MapParameters <PlainMap, bar, int>,
                      MapParameters <LinkedMap, bar, int>,
                      MapParameters <PlainMap, ham, int> >
         test_map_parameters_std_equal_std_allocator;

typedef  mpl::vector <
#                   if HAVE_BOOST_INTERPROCESS
                      MapParameters <PlainMap, int, foo,
                                     MCT_HASH_NAMESPACE::hash <int>, equal_to <int>,
                                     interprocess_allocator <pair <const int, foo> > >,
                      MapParameters <LinkedMap, foo, bar,
                                     MCT_HASH_NAMESPACE::hash <foo>, equal_to <foo>,
                                     interprocess_allocator <pair <const foo, bar> > >,
                      MapParameters <HugeLinkedMap, bar, int,
                                     MCT_HASH_NAMESPACE::hash <bar>, equal_to <bar>,
                                     interprocess_allocator <pair <const bar, int> > >,
                      MapParameters <ForwardMap, ham, int,
                                     MCT_HASH_NAMESPACE::hash <ham>, equal_to <ham>,
                                     interprocess_allocator <pair <const ham, int> > >,
                      MapParameters <HugeForwardMap, int, int,
                                     MCT_HASH_NAMESPACE::hash <int>, equal_to <int>,
                                     interprocess_allocator <pair <const int, int> > >
#                   endif
                     >
         test_map_parameters_std_equal_interp_allocator;

typedef  mpl::insert_range <test_map_parameters_std_equal_std_allocator,
                            mpl::end <test_map_parameters_std_equal_std_allocator>::type,
                            test_map_parameters_std_equal_interp_allocator>::type
         test_map_parameters_std_equal;

typedef  mpl::push_back <test_map_parameters_std_equal,
                         MapParameters <PlainMap, int, string, int_hasher_special,
                                        int_comparator_special> >::type
         test_map_parameters;

typedef  mpl::copy_if <test_map_parameters_std_equal, normal_parameters <mpl::_1> >::type
         test_map_parameters_std_equal_normal;

typedef  mpl::copy_if <test_map_parameters_std_equal, stable_order_parameters <mpl::_1> >::type
         test_map_parameters_std_equal_stable_order;

typedef  mpl::copy_if <test_map_parameters_std_equal, linked_parameters <mpl::_1> >::type
         test_map_parameters_std_equal_linked;

typedef  mpl::copy_if <test_map_parameters, normal_parameters <mpl::_1> >::type
         test_map_parameters_normal;

typedef  mpl::copy_if <test_map_parameters, stable_order_parameters <mpl::_1> >::type
         test_map_parameters_stable_order;

typedef  mpl::copy_if <test_map_parameters, linked_parameters <mpl::_1> >::type
         test_map_parameters_linked;

typedef  mpl::copy_if <test_map_parameters, forward_parameters <mpl::_1> >::type
         test_map_parameters_forward;

typedef  mpl::copy_if <test_map_parameters, robust_iterator_validation <mpl::_1> >::type
         test_map_parameters_robust_iterator_validation;

typedef  mpl::copy_if <test_map_parameters_normal, robust_iterator_validation <mpl::_1> >::type
         test_map_parameters_normal_robust_iterator_validation;

typedef  mpl::vector <LinkedMap, HugeLinkedMap, ForwardMap, HugeForwardMap>
         test_map_implementations_stable_order;

typedef  mpl::insert <test_map_implementations_stable_order,
                      mpl::begin <test_map_implementations_stable_order>::type, PlainMap>::type
         test_map_implementations;


#define COMMON_MAP_TEST_SETUP                                                   \
  typedef  typename parameters::implementation         implementation;          \
  typedef  typename parameters::wrapper_type           map_type;                \
  typedef  typename map_type::implementation_type      map_base_type;           \
  typedef  typename map_type::iterator                 iterator;                \
  typedef  typename map_type::const_iterator           const_iterator;          \
  typedef  typename parameters::key_type               key_type;                \
  typedef  typename parameters::value_type             value_type;              \
  typedef  typename parameters::assignable_value_type  assignable_value_type;   \
  typedef  typename parameters::mapped_type            mapped_type;             \
  typedef  typename parameters::hasher                 hasher;                  \
  typedef  typename parameters::equal_to               equal_to;                \
  typedef  Data <key_type>                             key_data;                \
  typedef  Data <assignable_value_type>                value_data;              \
  typedef  Data <mapped_type>                          mapped_data;


#if HAVE_BOOST_INTERPROCESS

// We use these magic numbers to avoid having to combine all the various (sets, maps) type
// into the same MPL collection.
# define TEST_MAP_TYPE_MAGIC_NUMBER     0x20000

typedef  mpl::transform <test_map_parameters_std_equal_interp_allocator,
                         get_tested <mpl::_1> >::type
         test_map_types_interp_allocator;

template <typename type>
struct external_validator
  <type,
   typename enable_if <is_map <type>::value
                       && is_interprocess_allocator <typename type::allocator_type>::value>
   ::type>
{
  static  void
  validate (const type& object)
  {
    typedef  typename mpl::find <test_map_types_interp_allocator, type>::type  position;
    typedef  mpl::begin <test_map_types_interp_allocator>::type                begin;

    validate_externally (&object, (TEST_MAP_TYPE_MAGIC_NUMBER
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
