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


#ifndef MCT_TESTS_COMMON_HPP
#define MCT_TESTS_COMMON_HPP


#include "tests/config.hpp"

#undef  MCT_ENABLE_DEBUGGING
#define MCT_ENABLE_DEBUGGING    1
#undef  MCT_CHECK_PRECONDITIONS
#define MCT_CHECK_PRECONDITIONS 1

// By default, MCT_CHECK_PRECONDITIONS disables iterators that don't refer to their table.
// However, in tests we want such iterators when possible, simply to make sure they work
// properly.
#undef  MCT_USE_EFFICIENT_ITERATORS
#define MCT_USE_EFFICIENT_ITERATORS     1

#include <mct/hash-map.hpp>
#include <mct/hash-set.hpp>

#include <algorithm>
#include <map>
#include <ostream>
#include <set>
#include <utility>
#include <vector>

#if MCT_CXX0X_SUPPORTED
# include <initializer_list>
#endif

// We use deprecated headers in order to work on ancient Boost versions.
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>

#if HAVE_BOOST_INTERPROCESS
# include <boost/interprocess/allocators/allocator.hpp>
# include <boost/interprocess/managed_shared_memory.hpp>
#endif  // HAVE_BOOST_INTERPROCESS


// A rather large portability section intended to make MCT testable with ancient Boost
// versions.  Currently I have access to a machine with GCC 4.1 and Boost 1.33.1 installed
// and interested in testing there.

// Don't care about nice failure formatting, just need it to compile.
#ifndef BOOST_CHECK_GE
# define BOOST_CHECK_GE(left, right)    BOOST_CHECK ((left) >= (right))
#endif

#ifndef BOOST_CHECK_LE
# define BOOST_CHECK_LE(left, right)    BOOST_CHECK ((left) <= (right))
#endif

#ifndef BOOST_CHECK_NE
# define BOOST_CHECK_NE(left, right)    BOOST_CHECK ((left) != (right))
#endif

#ifndef BOOST_TEST_CHECKPOINT
# define BOOST_TEST_CHECKPOINT(text)
#endif

#include <boost/version.hpp>

#if BOOST_VERSION / 100 == 1033

// Just add end semicolons.  I'd hate to have them in our source code.
#undef  BOOST_AUTO_TEST_SUITE
#define BOOST_AUTO_TEST_SUITE( suite_name )                             \
BOOST_AUTO_TC_REGISTRAR( suite_name )( BOOST_TEST_SUITE(                \
    BOOST_STRINGIZE( suite_name ) ) );

#undef  BOOST_AUTO_TEST_SUITE_END
#define BOOST_AUTO_TEST_SUITE_END()                                     \
BOOST_AUTO_TC_REGISTRAR( BOOST_JOIN( end_suite, __LINE__ ) )( 1 );

// Boost 1.33 has a typo there.
# undef  BOOST_AUTO_TEST_CASE_TEMPLATE
# define BOOST_AUTO_TEST_CASE_TEMPLATE( test_name, type_name, TL )      \
template<typename type_name>                                            \
void test_name( boost::type<type_name>* );                              \
                                                                        \
struct BOOST_AUTO_TC_INVOKER( test_name ) {                             \
    template<typename TestType>                                         \
    static void run( boost::type<TestType>* frwrd = 0 )                 \
    {                                                                   \
       test_name( frwrd );                                              \
    }                                                                   \
};                                                                      \
                                                                        \
BOOST_AUTO_TC_REGISTRAR( test_name )(                                   \
    boost::unit_test::ut_detail::template_test_case_gen<                \
        BOOST_AUTO_TC_INVOKER( test_name ),TL >(                        \
          BOOST_STRINGIZE( test_name ) ) );                             \
                                                                        \
template<typename type_name>                                            \
void test_name( boost::type<type_name>* )                               \
/**/

#endif


using mct::closed_hash_map;
using mct::closed_hash_set;
using mct::external_use;
using mct::forward_hash_map;
using mct::forward_hash_set;
using mct::huge_forward_hash_map;
using mct::huge_forward_hash_set;
using mct::huge_linked_hash_map;
using mct::huge_linked_hash_set;
using mct::intrusive_storage;
using mct::is_closed;
using mct::is_forward;
using mct::is_linked;
using mct::is_map;
using mct::is_set;
using mct::linked_hash_map;
using mct::linked_hash_set;
using mct::supports_external_use;
using mct::impl::enable_if;

using std::allocator;
using std::equal_to;
using std::logic_error;
using std::make_pair;
using std::map;
using std::max;
using std::min;
using std::ostream;
using std::pair;
using std::set;
using std::size_t;
using std::string;
using std::type_info;
using std::vector;

#if MCT_CXX0X_SUPPORTED
using std::forward;
using std::initializer_list;
using std::move;
#endif

namespace impl = mct::impl;
namespace mpl  = boost::mpl;


template <typename range_type>
size_t
range_size (const range_type& range)
{  return range.size ();  }

template <typename type, size_t n>
size_t
range_size (type (&) [n])
{  return n;  }

template <typename range_type>
typename range_type::iterator
range_begin (range_type& range)
{  return range.begin ();  }

template <typename range_type>
typename range_type::const_iterator
range_begin (const range_type& range)
{  return range.begin ();  }

template <typename type, size_t n>
type*
range_begin (type (& array) [n])
{  return array;  }

template <typename range_type>
typename range_type::iterator
range_end (range_type& range)
{  return range.end ();  }

template <typename range_type>
typename range_type::const_iterator
range_end (const range_type& range)
{  return range.end ();  }

template <typename type, size_t n>
type*
range_end (type (& array) [n])
{  return array + n;  }

#define RANGE(range)            range_begin (range), range_end (range)


// E.g. GCC 4.1 has 'long long' type, but no hash function specialization for it.
#if MCT_HAVE_LONG_LONG && !HAVE_LONG_LONG_HASH_SPECIALIZATION

MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <unsigned long long>
{
  size_t
  operator() (unsigned long long value) const
  {
    return static_cast <size_t> (value);
  }
};
MCT_HASH_NAMESPACE_LEAVE

#endif


// Saves typing.
template <typename Container, typename type>
inline  bool
contains (const Container& container, const type& value)
{
  return container.find (value) != container.end ();
}



// By default consider 'Map' to be just sequence of pairs.
template <typename Map, typename = void>
struct key_type_lookup
{
  typedef  typename Map::value_type::first_type  type;
};

// Actual maps.
template <typename Map>
struct key_type_lookup <Map, typename enable_if <!impl::is_same <typename Map::key_type,
                                                                 void>::value>::type>
{
  typedef  typename Map::key_type  type;
};

template <typename Map>
static  vector <typename key_type_lookup <Map>::type>
keys_of (const Map& map)
{
  vector <typename key_type_lookup <Map>::type>  keys;
  for (typename Map::const_iterator scan = map.begin (); scan != map.end (); ++scan)
    keys.push_back (scan->first);

  return keys;
}


template <typename InputIterator>
ostream&
print_range (ostream& stream, const char* delimiters, InputIterator begin, InputIterator end)
{
  stream << delimiters[0];
  if (begin != end)
    stream << ' ';

  bool first = true;
  for (; begin != end; ++begin)
    {
      if (first)
        first = false;
      else
        stream << ", ";

      stream << *begin;
    }

  if (!first)
    stream << ' ';
  stream << delimiters[1];

  return stream;
}

template <typename InputIterator>
ostream&
print_map_range (ostream& stream, const char* delimiters, InputIterator begin, InputIterator end)
{
  stream << delimiters[0];
  if (begin != end)
    stream << ' ';

  bool first = true;
  for (; begin != end; ++begin)
    {
      if (first)
        first = false;
      else
        stream << ", ";

      stream << begin->first << ": " << begin->second;
    }

  if (!first)
    stream << ' ';
  stream << delimiters[1];

  return stream;
}


template <typename type>
vector <type>
create_vector (type (* generator) (int), int stop, int start = 0)
{
  vector <type>  result;

  if (stop > start)
    {
      result.reserve (stop - start);
      for (int k = start; k < stop; ++k)
        result.push_back (generator (k));
    }

  return result;
}


namespace mct
{

  namespace impl
  {

    template <typename Bucket>
    ostream&
    operator<< (ostream& stream, const hash_table_const_iterator <Bucket>& iter)
    {
      // Don't use operator*() on the iterator, since that validates iterator and we don't
      // need exceptions here: presumably, if printing is going on, something is already
      // being reported.
      stream << "[iterator: " << &iter.bucket ()->value () << "]";
      return stream;
    }

  }

}


namespace std
{

  template <typename type>
  ostream&
  operator<< (ostream& stream, const vector <type>& array)
  {
    return print_range (stream, "[]", array.begin (), array.end ());
  }

  template <typename first_type, typename second_type>
  ostream&
  operator<< (ostream& stream, const pair <first_type, second_type>& pair)
  {
    return stream << "<" << pair.first << ", " << pair.second << ">";
  }

}


// This one can be used to test if item destructor is called when it should be.  In
// addition, it has fancy size (3 bytes normally) fancy comparator and fancy hash
// function.
struct foo
{
  static  size_t  num_alive;

  unsigned char  byte1;
  unsigned char  byte2;
  unsigned char  byte3;

  explicit
  foo (unsigned char byte1 = 0, unsigned char byte2 = 0, unsigned char byte3 = 0)
    : byte1 (byte1),
      byte2 (byte2),
      byte3 (byte3)
  {
    ++num_alive;
  }

  foo (const foo& that)
    : byte1 (that.byte1),
      byte2 (that.byte2),
      byte3 (that.byte3)
  {
    ++num_alive;
  }

  ~foo ()
  {
    --num_alive;
  }
};


inline  bool
operator== (const foo& foo1, const foo& foo2)
{
  return (min    (foo1.byte1, foo1.byte2) == min (foo2.byte1, foo2.byte2)
          && max (foo1.byte1, foo1.byte2) == max (foo2.byte1, foo2.byte2)
          && foo1.byte3 == foo2.byte3);
}


MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <foo>
{
  size_t
  operator() (const foo& _foo) const
  {
    return min (_foo.byte1, _foo.byte2) + 557 * (_foo.byte3 & 0x9b);
  }
};
MCT_HASH_NAMESPACE_LEAVE


inline  ostream&
operator<< (ostream& stream, const foo& _foo)
{
  stream << "("
         << static_cast <unsigned int> (_foo.byte1) << ", "
         << static_cast <unsigned int> (_foo.byte2) << ", "
         << static_cast <unsigned int> (_foo.byte3) << ")";
  return stream;
}


// This type is aimed at catching a specific bug discovered after 0.9.5 when compiling
// with GCC -O2 or -O3.
struct bar
{
  int  value : 16;

  explicit
  bar (int value = 0)
    : value (value)
  { }
};


inline  bool
operator== (const bar& bar1, const bar& bar2)
{
  return bar1.value == bar2.value;
}


MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <bar> : public hash <int>
{
  size_t
  operator() (const bar& _bar) const
  {
    return hash <int>::operator() (_bar.value);
  }
};
MCT_HASH_NAMESPACE_LEAVE


inline  ostream&
operator<< (ostream& stream, const bar& _bar)
{
  stream << '[' << _bar.value << ']';
  return stream;
}


// Structure to test optional intrusion support introduced in 1.1.3.
struct ham
{
  int   value_int;
  char  value_char;
  char  _extern_use_field;

  explicit
  ham (int value_int = 0, char value_char = 0)
    : value_int  (value_int),
      value_char (value_char)
  { }

  // We _must_ provide our own operator= that doesn't touch '_extern_use_field'.
  ham&
  operator= (const ham& that)
  {
    this->value_int  = that.value_int;
    this->value_char = that.value_char;
    return *this;
  }
};


inline  bool
operator== (const ham& ham1, const ham& ham2)
{
  return ham1.value_int == ham2.value_int && ham1.value_char == ham2.value_char;
}


MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <ham> : hash <int>, hash <char>
{
  size_t
  operator() (const ham& _ham) const
  {
    return hash <int>::operator() (_ham.value_int) ^ hash <char>::operator() (_ham.value_char);
  }
};
MCT_HASH_NAMESPACE_LEAVE


inline  ostream&
operator<< (ostream& stream, const ham& _ham)
{
  stream << '[' << _ham.value_int << ", " << +_ham.value_char << ']';
  return stream;
}


namespace mct
{
  template <>
  struct external_use <ham> : extern_use_field <ham, char, &ham::_extern_use_field>
  { };
}


struct int_wrapper
{
  // Not really needed.  Only used to make the type not have a trivial destructor.
  static  size_t  num_destroyed;

  int  value;

  explicit
  int_wrapper (int value = 0)
  {  this->value = value;  }

  ~int_wrapper ()
  {  ++num_destroyed;  }

  static  int_wrapper
  wrap (int value)
  {  return int_wrapper (value);  }
};


inline  bool
operator== (const int_wrapper& wrapper1, const int_wrapper& wrapper2)
{
  return wrapper1.value == wrapper2.value;
}


MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <int_wrapper> : public hash <int>
{
  size_t
  operator() (const int_wrapper& wrapper) const
  {
    return hash <int>::operator() (wrapper.value);
  }
};
MCT_HASH_NAMESPACE_LEAVE


inline  ostream&
operator<< (ostream& stream, const int_wrapper& wrapper)
{
  stream << '[' << wrapper.value << ']';
  return stream;
}


struct expected_exception : public std::exception
{
  virtual  const char*
  what () const throw ()
  {  return "an expected exception";  }

  static  void
  maybe_throw (int& countdown)
  {
    if (countdown > 0)
      --countdown;
    else if (countdown == 0)
      throw expected_exception ();
  }
};


template <int& countdown>
class countdown_setter
{
  const int  _revert_to;

public:

  countdown_setter (int value = 0)
    : _revert_to (countdown)
  {  countdown = value;  }

  ~countdown_setter ()
  {  countdown = _revert_to;  }
};


struct throws_on_copy : public int_wrapper
{
  static  int  countdown;

  typedef  countdown_setter <countdown>  enabler;

  explicit
  throws_on_copy (int value = 0)
    : int_wrapper (value)
  { }

  throws_on_copy (const throws_on_copy& that)
    : int_wrapper (that)
  {
    expected_exception::maybe_throw (countdown);
  }
};


MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <throws_on_copy> : public hash <int_wrapper>
{ };
MCT_HASH_NAMESPACE_LEAVE


struct throws_in_comparator : public int_wrapper
{
  static  int  countdown;

  typedef  countdown_setter <countdown>  enabler;

  explicit
  throws_in_comparator (int value = 0)
    : int_wrapper (value)
  { }
};


inline  bool
operator== (const throws_in_comparator& wrapper1, const throws_in_comparator& wrapper2)
{
  expected_exception::maybe_throw (throws_in_comparator::countdown);
  return wrapper1.value == wrapper2.value;
}


MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <throws_in_comparator> : public hash <int_wrapper>
{ };
MCT_HASH_NAMESPACE_LEAVE


struct throws_in_hasher : public int_wrapper
{
  static  int  countdown;

  typedef  countdown_setter <countdown>  enabler;

  explicit
  throws_in_hasher (int value = 0)
    : int_wrapper (value)
  { }
};


MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <throws_in_hasher> : public hash <int_wrapper>
{
  size_t
  operator() (const int_wrapper& wrapper) const
  {
    expected_exception::maybe_throw (throws_in_hasher::countdown);
    return hash <int>::operator() (wrapper.value);
  }
};
MCT_HASH_NAMESPACE_LEAVE


struct throwing_int_strict_weak_ordering
{
  static  int  countdown;

  typedef  countdown_setter <countdown>  enabler;

  bool
  operator() (int x, int y) const
  {
    expected_exception::maybe_throw (throwing_int_strict_weak_ordering::countdown);
    return x < y;
  }

  // For comparing elements of maps by their keys only.
  template <typename Mapped>
  bool
  operator() (const pair <const int, Mapped>& x, const pair <const int, Mapped>& y) const
  {
    expected_exception::maybe_throw (throwing_int_strict_weak_ordering::countdown);
    return x.first < y.first;
  }
};


#if MCT_CXX0X_SUPPORTED

struct must_not_be_copied : public int_wrapper
{
  explicit
  must_not_be_copied (int value = 0)
    : int_wrapper (value)
  { }

  must_not_be_copied (const must_not_be_copied& /* that */)
  {
    throw logic_error ("copy constructor used where it mustn't be");
  }

  must_not_be_copied (must_not_be_copied&& that)
    : int_wrapper (that)
  { }
};


MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <must_not_be_copied> : public hash <int_wrapper>
{ };
MCT_HASH_NAMESPACE_LEAVE


struct must_not_be_moved : public int_wrapper
{
  explicit
  must_not_be_moved (int value = 0)
    : int_wrapper (value)
  { }

  must_not_be_moved (const must_not_be_moved& that)
    : int_wrapper (that)
  { }

  must_not_be_moved (must_not_be_moved&& /* that */)
  {
    throw logic_error ("move constructor used where it mustn't be");
  }
};


MCT_HASH_NAMESPACE_ENTER
template <>
struct hash <must_not_be_moved> : public hash <int_wrapper>
{ };
MCT_HASH_NAMESPACE_LEAVE

#endif  // MCT_CXX0X_SUPPORTED


#if MCT_CXX0X_SUPPORTED


// Various types below should be subclasses of 'int_wrapper' or something similar.

#define INITIALIZERS_0_50(type)                                 \
  { type ( 0), type ( 1), type ( 2), type ( 3), type ( 4),      \
    type ( 5), type ( 6), type ( 7), type ( 8), type ( 9),      \
    type (10), type (11), type (12), type (13), type (14),      \
    type (15), type (16), type (17), type (18), type (19),      \
    type (20), type (21), type (22), type (23), type (24),      \
    type (25), type (26), type (27), type (28), type (29),      \
    type (30), type (31), type (32), type (33), type (34),      \
    type (35), type (36), type (37), type (38), type (39),      \
    type (40), type (41), type (42), type (43), type (44),      \
    type (45), type (46), type (47), type (48), type (49) }

#define INITIALIZERS_40_60(type)                                \
  { type (40), type (41), type (42), type (43), type (44),      \
    type (45), type (46), type (47), type (48), type (49),      \
    type (50), type (51), type (52), type (53), type (54),      \
    type (55), type (56), type (57), type (58), type (59) }

#define MAP_INITIALIZERS_0_50(key_type, mapped_type)                                            \
 { make_pair (key_type ( 0), mapped_type ( 0)), make_pair (key_type ( 1), mapped_type ( 1)),    \
   make_pair (key_type ( 2), mapped_type ( 2)), make_pair (key_type ( 3), mapped_type ( 3)),    \
   make_pair (key_type ( 4), mapped_type ( 4)), make_pair (key_type ( 5), mapped_type ( 5)),    \
   make_pair (key_type ( 6), mapped_type ( 6)), make_pair (key_type ( 7), mapped_type ( 7)),    \
   make_pair (key_type ( 8), mapped_type ( 8)), make_pair (key_type ( 9), mapped_type ( 9)),    \
   make_pair (key_type (10), mapped_type (10)), make_pair (key_type (11), mapped_type (11)),    \
   make_pair (key_type (12), mapped_type (12)), make_pair (key_type (13), mapped_type (13)),    \
   make_pair (key_type (14), mapped_type (14)), make_pair (key_type (15), mapped_type (15)),    \
   make_pair (key_type (16), mapped_type (16)), make_pair (key_type (17), mapped_type (17)),    \
   make_pair (key_type (18), mapped_type (18)), make_pair (key_type (19), mapped_type (19)),    \
   make_pair (key_type (20), mapped_type (20)), make_pair (key_type (21), mapped_type (21)),    \
   make_pair (key_type (22), mapped_type (22)), make_pair (key_type (23), mapped_type (23)),    \
   make_pair (key_type (24), mapped_type (24)), make_pair (key_type (25), mapped_type (25)),    \
   make_pair (key_type (26), mapped_type (26)), make_pair (key_type (27), mapped_type (27)),    \
   make_pair (key_type (28), mapped_type (28)), make_pair (key_type (29), mapped_type (29)),    \
   make_pair (key_type (30), mapped_type (30)), make_pair (key_type (31), mapped_type (31)),    \
   make_pair (key_type (32), mapped_type (32)), make_pair (key_type (33), mapped_type (33)),    \
   make_pair (key_type (34), mapped_type (34)), make_pair (key_type (35), mapped_type (35)),    \
   make_pair (key_type (36), mapped_type (36)), make_pair (key_type (37), mapped_type (37)),    \
   make_pair (key_type (38), mapped_type (38)), make_pair (key_type (39), mapped_type (39)),    \
   make_pair (key_type (40), mapped_type (40)), make_pair (key_type (41), mapped_type (41)),    \
   make_pair (key_type (42), mapped_type (42)), make_pair (key_type (43), mapped_type (43)),    \
   make_pair (key_type (44), mapped_type (44)), make_pair (key_type (45), mapped_type (45)),    \
   make_pair (key_type (46), mapped_type (46)), make_pair (key_type (47), mapped_type (47)),    \
   make_pair (key_type (48), mapped_type (48)), make_pair (key_type (49), mapped_type (49)) }

#define MAP_INITIALIZERS_40_60(key_type, mapped_type)                                           \
 { make_pair (key_type (40), mapped_type (40)), make_pair (key_type (41), mapped_type (41)),    \
   make_pair (key_type (42), mapped_type (42)), make_pair (key_type (43), mapped_type (43)),    \
   make_pair (key_type (44), mapped_type (44)), make_pair (key_type (45), mapped_type (45)),    \
   make_pair (key_type (46), mapped_type (46)), make_pair (key_type (47), mapped_type (47)),    \
   make_pair (key_type (48), mapped_type (48)), make_pair (key_type (49), mapped_type (49)),    \
   make_pair (key_type (50), mapped_type (50)), make_pair (key_type (51), mapped_type (51)),    \
   make_pair (key_type (52), mapped_type (52)), make_pair (key_type (53), mapped_type (53)),    \
   make_pair (key_type (54), mapped_type (54)), make_pair (key_type (55), mapped_type (55)),    \
   make_pair (key_type (56), mapped_type (56)), make_pair (key_type (57), mapped_type (57)),    \
   make_pair (key_type (58), mapped_type (58)), make_pair (key_type (59), mapped_type (59)) }


#endif  // MCT_CXX0X_SUPPORTED


#if HAVE_BOOST_INTERPROCESS

namespace interp = boost::interprocess;


struct interp_memory_chunk
{
  const bool                     main_user;
  interp::managed_shared_memory  chunk;

public:

  interp_memory_chunk (bool main_user)
    : main_user (main_user)
  {
    if (main_user)
      {
        interp::shared_memory_object::remove ("MCT_interprocess_test");
        chunk = interp::managed_shared_memory (interp::create_only,
                                               "MCT_interprocess_test", 1024U * 1024U);
      }
    else
      chunk = interp::managed_shared_memory (interp::open_only, "MCT_interprocess_test");
  }

  ~interp_memory_chunk ()
  {
    if (main_user)
      interp::shared_memory_object::remove ("MCT_interprocess_test");
  }
};

extern  interp_memory_chunk  interp_memory;


template <typename type>
class interprocess_allocator
  : public interp::allocator <type, interp::managed_shared_memory::segment_manager>
{
  typedef interp::allocator <type, interp::managed_shared_memory::segment_manager>
          base_type;


public:

  template <typename type2>
  struct rebind
  {
    typedef  interprocess_allocator <type2>  other;
  };


  interprocess_allocator ()
    : base_type (interp_memory.chunk.get_segment_manager ())
  { }

  explicit
  interprocess_allocator (const base_type& base)
    : base_type (base)
  { }

  interprocess_allocator (const interprocess_allocator& that)
    : base_type (that)
  { }

  template <typename base>
  interprocess_allocator (const interprocess_allocator <base>& that)
    : base_type (that)
  { }
};

#endif  // HAVE_BOOST_INTERPROCESS


class leak_test_data
{
  static  int  num_allocators;

protected:

  static  set <pair <const type_info*, pair <void*, size_t> > >  allocated_chunks;
  static  set <pair <const type_info*, void*> >                  constructed_objects;

  static  void
  new_allocator ()
  {  ++num_allocators;  }

  static  void  delete_allocator ();
};


template <typename Wrapped>
struct leak_test_allocator : public Wrapped, protected leak_test_data
{
  typedef  Wrapped  base_type;

  typedef  typename base_type::value_type  value_type;
  typedef  typename base_type::pointer     pointer;
  typedef  typename base_type::size_type   size_type;


  template <typename type>
  struct rebind
  {
    typedef  leak_test_allocator <typename base_type::template rebind <type>::other>  other;
  };


  leak_test_allocator ()
  {
    leak_test_data::new_allocator ();
  }

  explicit
  leak_test_allocator (const base_type& base)
    : base_type (base)
  {
    leak_test_data::new_allocator ();
  }

  leak_test_allocator (const leak_test_allocator& that)
    : base_type (that)
  {
    leak_test_data::new_allocator ();
  }

  template <typename base>
  leak_test_allocator (const leak_test_allocator <base>& that)
    : base_type (that)
  {
    leak_test_data::new_allocator ();
  }

  ~leak_test_allocator ()
  {
    leak_test_data::delete_allocator ();
  }


  pointer
  allocate (size_type size)
  {
    // Tests are not supposed to cause huge allocations and this check might prevent a
    // buggy implementation to go out of hands.
    if (size > 1000000)
      throw std::bad_alloc ();

    const pointer  chunk = base_type::allocate (size);

    allocated_chunks.insert (make_pair (&typeid (chunk), make_pair (&*chunk, size)));
    return chunk;
  }

  void
  deallocate (pointer chunk, size_type size)
  {
    allocated_chunks.erase (make_pair (&typeid (chunk), make_pair (&*chunk, size)));
    base_type::deallocate (chunk, size);
  }

  void
  construct (pointer element, const value_type& value)
  {
    base_type::construct (element, value);
    if (!impl::has_trivial_destructor <value_type>::value)
      constructed_objects.insert (make_pair (&typeid (element), &*element));
  }

#if MCT_CXX0X_SUPPORTED

  void
  construct (pointer element, value_type&& value)
  {
    base_type::construct (element, forward <value_type> (value));
    if (!impl::has_trivial_destructor <value_type>::value)
      constructed_objects.insert (make_pair (&typeid (element), &*element));
  }

#endif  // MCT_CXX0X_SUPPORTED

  void
  destroy (pointer element)
  {
    if (!impl::has_trivial_destructor <value_type>::value)
      {
        if (!constructed_objects.erase (make_pair (&typeid (element), &*element)))
          throw logic_error ("destroying an object that was never (successfully) constructed");
      }

    base_type::destroy (element);
  }
};


struct int_hasher_1
{
  size_t
  operator() (int value) const
  {  return value + 3;;  }
};

struct int_hasher_special
{
  // Intentionally 'bad' hash.
  size_t
  operator() (int value) const
  {  return value & 0x3;  }
};

struct int_comparator_special
{
  bool
  operator() (int a, int b) const
  {  return (a & 0x7) == (b & 0x7);  }
};

struct string_hasher_1 : public MCT_HASH_NAMESPACE::hash <string>
{
  size_t
  operator() (const string& value) const
  {
    return MCT_HASH_NAMESPACE::hash <string>::operator() (value) ^ 1337;
  }
};

#if MCT_HAVE_LONG_LONG

struct unsigned_long_long_hasher
{
  size_t
  operator() (unsigned long long value) const
  {
    return (value & 0xffffffffu) + (value >> 32);
  }
};

#endif  // MCT_HAVE_LONG_LONG


template <typename key_type, typename value_type>
struct key_extractor
{
  static  const key_type&
  extract (const value_type& value)
  {  return value.first;  }
};

template <typename key_type>
struct key_extractor <key_type, key_type>
{
  static  const key_type&
  extract (const key_type& value)
  {  return value;  }
};


template <typename Container, typename Initializer>
void
assert_identical_order (const Container& container, const Initializer& initializer)
{
  typedef  typename Container::key_type    key_type;
  typedef  typename Container::value_type  value_type;
  typedef  map <int, key_type>             order_map;

  vector <key_type>    actual_order;
  vector <key_type>    expected_order;
  order_map            order;

  for (typename Container::const_iterator scan = container.begin ();
       scan != container.end (); ++scan)
    {
      const key_type&  key = key_extractor <key_type, value_type>::extract (*scan);
      actual_order.push_back (key);

      int  index = 0;
      for (typename Initializer::const_iterator initializer_scan = initializer.begin ();
           initializer_scan != initializer.end (); ++initializer_scan)
        {
          if (container.key_eq () (key, *initializer_scan))
            {
              order.insert (make_pair (index, *initializer_scan));
              break;
            }

          ++index;
        }
    }

  for (typename order_map::const_iterator scan = order.begin (); scan != order.end (); ++scan)
    expected_order.push_back (scan->second);

  BOOST_CHECK_EQUAL (actual_order, expected_order);
}


template <typename type>
struct Data;


template <>
struct Data <int>
{
  static  const vector <int>&
  values1 ()
  {
    static  const int  as_array[] = { 10, -6, 1, 3 };
    static  const vector <int>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <int>&
  values2 ()
  {
    static  const int  as_array[] = { 0, -3, 10, 8, 12, 1, -9 };
    static  const vector <int>  as_vector (RANGE (as_array));

    return as_vector;
  }

  // Currently only used in sort() function tests and is valid only for this 'Data'
  // specialization.  The intention is to have many shuffled values.
  static  const vector <int>&
  values3 ()
  {
    static  const int  as_array[] = {  15, 48, 16, 11, 45, 41, 30, 35, 31, 43,
                                       22, 49, 18,  6, 32, 20, 19, 37, 44, 26,
                                        7,  1, 24, 29, 40, 13,  0, 42, 33,  4,
                                        5, 47, 23, 17, 38,  8, 46, 21, 14, 36,
                                       12, 34, 10,  3, 25,  9, 39,  2, 28, 27 };
    static  const  vector <int>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  int
  generate (int iteration)
  {  return iteration;  }
};


template <>
struct Data <string>
{
  static  const vector <string>&
  values1 ()
  {
    static  const string  as_array[] = { "mon", "tue", "wed", "thu", "fri", "sat", "sun" };
    static  const vector <string>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <string>&
  values2 ()
  {
    static  const string  as_array[] = { "jan", "feb", "mar", "apr", "may", "etc" };
    static  const vector <string>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <int>&
  values3 ()
  {  throw std::logic_error ("not implemented");  }

  static  string
  generate (int iteration)
  {
    string  result;
    result.resize (iteration, '*');
    return result;
  }
};


template <>
struct Data <foo>
{
  static  const vector <foo>&
  values1 ()
  {
    static  const foo  as_array[]
      = { foo (0, 15, 3), foo (3, 8, 1), foo (9, 2, 2), foo (6, 0, 0), foo (1, 0, 7) };
    static  const vector <foo>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <foo>&
  values2 ()
  {
    static  const foo  as_array[]
      = { foo (1, 2, 1), foo (2, 8, 17), foo (3, 3, 3), foo (19, 0, 2) };
    static  const vector <foo>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <int>&
  values3 ()
  {  throw std::logic_error ("not implemented");  }

  static  foo
  generate (int iteration)
  {
    unsigned int   first = (iteration / 0x100);
    unsigned char  byte3 = (iteration & 0xff);
    unsigned int   byte1 = 0;

    while (first > byte1)
      first -= ++byte1;

    unsigned int  byte2 = first;

    if (byte1 >= 0x100 || byte2 >= 0x100)
      throw logic_error ("iteration is too large");

    return foo (byte1, byte2, byte3);
  }
};


template <>
struct Data <bar>
{
  static  const vector <bar>&
  values1 ()
  {
    static  const bar  as_array[] = { bar (42), bar (0), bar (3), bar (-3) };
    static  const vector <bar>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <bar>&
  values2 ()
  {
    static  const bar  as_array[] = { bar (7), bar (13), bar (1), bar (0), bar (88), bar (14) };
    static  const vector <bar>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <int>&
  values3 ()
  {  throw std::logic_error ("not implemented");  }

  static  bar
  generate (int iteration)
  {  return bar (iteration);  }
};


template <>
struct Data <ham>
{
  static  const vector <ham>&
  values1 ()
  {
    static  const ham  as_array[]
      = { ham (113, 'x'), ham (0, '^'), ham (0, '\0'), ham (-37, 'x') };
    static  const vector <ham>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <ham>&
  values2 ()
  {
    static  const ham  as_array[]
      = { ham (66, '0'), ham (-4, '7'), ham (91, '~'), ham (91, '/'), ham (61, '/'),
          ham (812, '7') };
    static  const vector <ham>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <int>&
  values3 ()
  {  throw std::logic_error ("not implemented");  }

  static  ham
  generate (int iteration)
  {  return ham (iteration, iteration);  }
};


#if MCT_HAVE_LONG_LONG

template <>
struct Data <unsigned long long>
{
  static  const vector <unsigned long long>&
  values1 ()
  {
    static  const unsigned long long  as_array[] = { 10, -6, 1, 3 };
    static  const vector <unsigned long long>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <unsigned long long>&
  values2 ()
  {
    static  const unsigned long long  as_array[] = { 0, -3, 10, 8, 12, 1, -9 };
    static  const vector <unsigned long long>  as_vector (RANGE (as_array));

    return as_vector;
  }

  static  const vector <int>&
  values3 ()
  {  throw std::logic_error ("not implemented");  }

  static  unsigned long long
  generate (int iteration)
  {  return iteration;  }
};

#endif  // MCT_HAVE_LONG_LONG


template <typename key_type, typename mapped_type>
struct Data <pair <key_type, mapped_type> >
{
  static  const vector <pair <key_type, mapped_type> >&
  values1 ()
  {
    static  const vector <pair <key_type, mapped_type> >  as_vector
      (combine (Data <key_type>::values1 (), Data <mapped_type>::values1 ()));

    return as_vector;
  }

  static  const vector <pair <key_type, mapped_type> >&
  values2 ()
  {
    static  const vector <pair <key_type, mapped_type> >  as_vector
      (combine (Data <key_type>::values2 (), Data <mapped_type>::values2 ()));

    return as_vector;
  }

  static  const vector <pair <key_type, mapped_type> >&
  values3 ()
  {
    static  const vector <pair <key_type, mapped_type> >  as_vector
      (combine (Data <key_type>::values3 (), Data <mapped_type>::values3 ()));

    return as_vector;
  }

  static  const vector <key_type>&
  keys1 ()
  {
    static  const vector <key_type>  as_vector
      (truncate (Data <key_type>::values1 (), Data <mapped_type>::values1 ()));

    return as_vector;
  }

  static  const vector <key_type>&
  keys2 ()
  {
    static  const vector <key_type>  as_vector
      (truncate (Data <key_type>::values2 (), Data <mapped_type>::values2 ()));

    return as_vector;
  }

  static  const vector <key_type>&
  keys3 ()
  {
    static  const vector <key_type>  as_vector
      (truncate (Data <key_type>::values3 (), Data <mapped_type>::values2 ()));

    return as_vector;
  }

  static  pair <key_type, mapped_type>
  generate (int iteration)
  {
    return make_pair (Data <key_type>   ::generate (iteration),
                      Data <mapped_type>::generate (iteration));
  }


private:

  static  const vector <pair <key_type, mapped_type> >
  combine (const vector <key_type>& left, const vector <mapped_type>& right)
  {
    vector <pair <key_type, mapped_type> >  result;

    for (size_t k = 0, size = min (left.size (), right.size ()); k < size; ++k)
      result.push_back (make_pair (left[k], right[k]));

    return result;
  }

  static  const vector <key_type>
  truncate (const vector <key_type>& left, const vector <mapped_type>& right)
  {
    vector <key_type>  result (left);

    result.resize (min (left.size (), right.size ()));
    return result;
  }
};


// A function object that intentionally doesn't order many numbers against each other:
// meant for testing sort stability.
struct even_first
{
  bool
  operator() (int x, int y) const
  {  return (x & 1) < (y & 1);  }

  // For comparing elements of maps by their keys only.
  template <typename Mapped>
  bool
  operator() (const pair <const int, Mapped>& x, const pair <const int, Mapped>& y) const
  {  return (x.first & 1) < (y.first & 1);  }

  // Same for sequences of pairs.
  template <typename Mapped>
  bool
  operator() (const pair <int, Mapped>& x, const pair <int, Mapped>& y) const
  {  return (x.first & 1) < (y.first & 1);  }
};


template <typename type, typename = void>
struct named_object
{
  type* const  object;


  named_object ()
    : object (new type ())
  { }

#if MCT_CXX0X_SUPPORTED

  // Needed because "perfect forwarding" is not that perfect after all: troubles
  // forwarding initializer_list at least on GCC 4.4.  In our usecases that can come as
  // the first argument only.
  template <typename value_type, typename... Args>
  named_object (initializer_list <value_type> initializers, Args&&... args)
    : object (new type (initializers, std::forward <Args&&> (args)...))
  { }

  template <typename... Args>
  named_object (Args&&... args)
    : object (new type (std::forward <Args&&> (args)...))
  { }

#else  // not MCT_CXX0X_SUPPORTED

  template <typename Arg1>
  named_object (const Arg1& arg1)
    : object (new type (arg1))
  { }

  template <typename Arg1, typename Arg2>
  named_object (const Arg1& arg1, const Arg2& arg2)
    : object (new type (arg1, arg2))
  { }

  template <typename Arg1, typename Arg2, typename Arg3>
  named_object (const Arg1& arg1, const Arg2& arg2, const Arg3& arg3)
    : object (new type (arg1, arg2, arg3))
  { }

  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  named_object (const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4)
    : object (new type (arg1, arg2, arg3, arg4))
  { }

  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5>
  named_object (const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4,
                const Arg5& arg5)
    : object (new type (arg1, arg2, arg3, arg4, arg5))
  { }

  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4,
            typename Arg5, typename Arg6>
  named_object (const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4,
                const Arg5& arg5, const Arg6& arg6)
    : object (new type (arg1, arg2, arg3, arg4, arg5, arg6))
  { }

#endif  // not MCT_CXX0X_SUPPORTED


  ~named_object ()
  {  delete object;  }


private:

  // Copy (and move) constructors/operator= must not be used.
  named_object (const named_object&);
  void  operator= (const named_object&);

#if MCT_CXX0X_SUPPORTED
  named_object (named_object&&);
  void  operator= (named_object&&);
#endif
};


#if HAVE_BOOST_INTERPROCESS


template <typename type>
struct is_interprocess_allocator : impl::false_type
{ };

template <typename type>
struct is_interprocess_allocator <interprocess_allocator <type> > : impl::true_type
{ };

template <typename allocator_type>
struct is_interprocess_allocator <leak_test_allocator <allocator_type> >
  : is_interprocess_allocator <allocator_type>
{ };


string  generate_unique_object_name ();


template <typename type>
struct named_object <type,
                     typename enable_if <is_interprocess_allocator <typename type::allocator_type>
                                         ::value>::type>
{
  const string  name;
  type* const   object;


  named_object ()
    : name   (generate_unique_object_name ()),
      object (interp_memory.chunk.construct <type> (name.c_str ()) ())
  { }

#if MCT_CXX0X_SUPPORTED

  template <typename value_type, typename... Args>
  named_object (initializer_list <value_type> initializers, Args&&... args)
    : name   (generate_unique_object_name ()),
      object (interp_memory.chunk
              // This call will give warnings at least on GCC 4.4, but I don't see what I
              // can do about that, as the relevant definitions are in Boost.
              .construct <type> (name.c_str ()) (initializers, std::forward <Args&&> (args)...))
  { }

  template <typename... Args>
  named_object (Args&&... args)
    : name   (generate_unique_object_name ()),
      object (interp_memory.chunk
              .construct <type> (name.c_str ()) (std::forward <Args&&> (args)...))
  { }

#else  // not MCT_CXX0X_SUPPORTED

  template <typename Arg1>
  named_object (const Arg1& arg1)
    : name   (generate_unique_object_name ()),
      object (interp_memory.chunk.construct <type> (name.c_str ()) (arg1))
  { }

  template <typename Arg1, typename Arg2>
  named_object (const Arg1& arg1, const Arg2& arg2)
    : name   (generate_unique_object_name ()),
      object (interp_memory.chunk.construct <type> (name.c_str ()) (arg1, arg2))
  { }

  template <typename Arg1, typename Arg2, typename Arg3>
  named_object (const Arg1& arg1, const Arg2& arg2, const Arg3& arg3)
    : name   (generate_unique_object_name ()),
      object (interp_memory.chunk.construct <type> (name.c_str ()) (arg1, arg2, arg3))
  { }

  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  named_object (const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4)
    : name   (generate_unique_object_name ()),
      object (interp_memory.chunk.construct <type> (name.c_str ()) (arg1, arg2, arg3, arg4))
  { }

  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5>
  named_object (const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4,
                const Arg5& arg5)
    : name   (generate_unique_object_name ()),
      object (interp_memory.chunk.construct <type> (name.c_str ()) (arg1, arg2, arg3, arg4, arg5))
  { }

  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4,
            typename Arg5, typename Arg6>
  named_object (const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4,
                const Arg5& arg5, const Arg6& arg6)
    : name   (generate_unique_object_name ()),
      object (interp_memory.chunk.construct <type> (name.c_str ()) (arg1, arg2, arg3, arg4,
                                                                     arg5, arg6))
  { }

#endif  // not MCT_CXX0X_SUPPORTED


  ~named_object ()
  {
    interp_memory.chunk.destroy <type> (name.c_str ());
  }


private:

  // Copy (and move) constructors/operator= must not be used.
  named_object (const named_object&);
  void  operator= (const named_object&);

#if MCT_CXX0X_SUPPORTED
  named_object (named_object&&);
  void  operator= (named_object&&);
#endif
};


void  validate_externally (const void* object, int type_index);

#endif  // HAVE_BOOST_INTERPROCESS


template <typename implementation_type, typename = void>
struct robust_for;

template <typename implementation_type, typename = void>
struct external_validator
{
  static  void
  validate (const implementation_type& /* implementation */)
  { }
};


template <typename parameters>
struct get_tested
{
  typedef  typename parameters::tested_type  type;
};


// In tests small memory/performance hit of virtual functions doesn't matter.
template <typename Implementation>
class test_wrapper
{
public:

  typedef  Implementation  implementation_type;

  virtual  const implementation_type&  impl () const = 0;

  virtual  void  validate () const = 0;

protected:

  typedef  typename robust_for <implementation_type>::type  robust_type;

  virtual  implementation_type&  impl  () = 0;
  virtual  robust_type&          valid () = 0;
  virtual  const robust_type&    valid () const = 0;
};


template <typename Implementation>
class test_table_interface : public test_wrapper <Implementation>
{
  typedef  test_wrapper <Implementation>  wrapper_type;

public:

  typedef  typename wrapper_type::implementation_type    implementation_type;
  typedef  typename implementation_type::key_type        key_type;
  typedef  typename implementation_type::value_type      value_type;
  typedef  typename implementation_type::key_equal       key_equal;
  typedef  typename implementation_type::iterator        iterator;
  typedef  typename implementation_type::const_iterator  const_iterator;
  typedef  typename implementation_type::size_type       size_type;


  iterator
  begin ()
  {  return this->impl ().begin ();  }

  const_iterator
  begin () const
  {  return this->impl ().begin ();  }

  const_iterator
  cbegin () const
  {  return this->impl ().cbegin ();  }

  iterator
  end ()
  {  return this->impl ().end ();  }

  const_iterator
  end () const
  {  return this->impl ().end ();  }

  const_iterator
  cend () const
  {  return this->impl ().cend ();  }


  size_type
  bucket_count () const
  {  return this->impl ().bucket_count ();  }

  bool
  empty () const
  {  return this->impl ().empty ();  }

  size_type
  size () const
  {  return this->impl ().size ();  }


  key_equal
  key_eq () const
  {  return this->impl ().key_eq ();  }


  iterator
  find (const key_type& key)
  {  return this->impl ().find (key);  }

  const_iterator
  find (const key_type& key) const
  {  return this->impl ().find (key);  }

  size_type
  count (const key_type& key) const
  {  return this->impl ().count (key);  }

  pair <iterator, iterator>
  equal_range (const key_type& key)
  {  return this->impl ().equal_range (key);  }

  pair <const_iterator, const_iterator>
  equal_range (const key_type& key) const
  {  return this->impl ().equal_range (key);  }


  pair <iterator, bool>
  insert (const value_type& value)
  {
    const pair <iterator, bool>  result (this->impl ().insert (value));
    this->valid ().insert (value);
    this->validate ();
    return result;
  }

  iterator
  insert (const_iterator hint, const value_type& value)
  {
    const iterator  result (this->impl ().insert (hint, value));
    this->valid ().insert (this->valid ().end (), value);
    this->validate ();
    return result;
  }

  template <typename InputIterator>
  void
  insert (InputIterator begin, InputIterator end)
  {
    this->impl  ().insert (begin, end);
    this->valid ().insert (begin, end);
    this->validate ();
  }


  void
  clear ()
  {
    this->impl  ().clear ();
    this->valid ().clear ();
    this->validate ();
  }


  void
  swap (test_table_interface& that)
  {
    this->impl  ().swap (that.impl  ());
    this->valid ().swap (that.valid ());
    this->validate ();
    that .validate ();
  }


  void
  max_load_factor (float max_load_factor)
  {
    this->impl ().max_load_factor (max_load_factor);
    this->validate ();
  }

  void
  rehash (size_t num_buckets)
  {
    this->impl ().rehash (num_buckets);
    this->validate ();
  }


  size_type
  used_memory () const
  {  return this->impl ().used_memory ();  }


  bool
  uses_dynamic_memory () const
  {
    // _Strictly_ speaking this is worse than testing '_data.buckets' directly, but should
    // be good enough unless some drastic changes are made to the implementation.  The
    // "proper" way is not easily possible because the field is private.
    return this->impl ().used_memory () > sizeof (implementation_type);
  }


protected:

  virtual  const key_type&  extract_key (const_iterator iterator) const = 0;
};


// Exists for closed/linked/forward specialization.
template <typename implementation_type, typename = typename implementation_type::_container_family>
class test_family_interface;


template <typename implementation_type>
class test_family_interface <implementation_type, impl::closed_family>
  : public test_table_interface <implementation_type>
{
  typedef  test_table_interface <implementation_type>  base_type;

public:

  typedef  typename base_type::key_type        key_type;
  typedef  typename base_type::iterator        iterator;
  typedef  typename base_type::const_iterator  const_iterator;


  bool
  erase (const key_type& key)
  {
    const bool  result = this->impl ().erase (key);
    BOOST_REQUIRE_EQUAL (result, this->valid ().erase (key));
    this->validate ();
    return result;
  }

  iterator
  erase (const_iterator position)
  {
    this->valid ().erase (this->extract_key (position));
    const iterator  result = this->impl ().erase (position);
    this->validate ();
    return result;
  }

  iterator
  erase (const_iterator first, const_iterator last)
  {
    // Order is not required to be the same; emulate.
    for (const_iterator scan = first; scan != last; ++scan)
      this->valid ().erase (this->extract_key (scan));

    const iterator  result = this->impl ().erase (first, last);
    this->validate ();
    return result;
  }
};


template <typename implementation_type>
class test_family_interface <implementation_type, impl::linked_family>
  : public test_family_interface <implementation_type, impl::closed_family>
{
  typedef  test_family_interface <implementation_type, impl::closed_family>  base_type;

public:

  typedef  typename base_type::value_type      value_type;
  typedef  typename base_type::iterator        iterator;
  typedef  typename base_type::const_iterator  const_iterator;


  value_type&
  front ()
  {  return this->impl ().front ();  }

  const value_type&
  front () const
  {  return this->impl ().front ();  }

  value_type&
  back ()
  {  return this->impl ().back ();  }

  const value_type&
  back () const
  {  return this->impl ().back ();  }


  void
  pop_front ()
  {
    this->valid ().erase (this->extract_key (this->impl ().begin ()));
    this->impl  ().pop_front ();
    this->validate ();
  }

  void
  pop_back ()
  {
    this->valid ().erase (this->extract_key (--this->impl ().end ()));
    this->impl  ().pop_back ();
    this->validate ();
  }


  void
  relink (const_iterator before, iterator element)
  {
    this->impl ().relink (before, element);
    this->validate ();
  }

  void
  reverse ()
  {
    this->impl ().reverse ();
    this->validate ();
  }

  void
  sort ()
  {
    this->impl ().sort ();
    this->validate ();
  }
};


template <typename implementation_type>
class test_family_interface <implementation_type, impl::forward_family>
  : public test_table_interface <implementation_type>
{
  typedef  test_table_interface <implementation_type>  base_type;

public:

  typedef  typename base_type::value_type      value_type;
  typedef  typename base_type::iterator        iterator;
  typedef  typename base_type::const_iterator  const_iterator;


  value_type&
  front ()
  {  return this->impl ().front ();  }

  const value_type&
  front () const
  {  return this->impl ().front ();  }

  value_type&
  back ()
  {  return this->impl ().back ();  }

  const value_type&
  back () const
  {  return this->impl ().back ();  }


  iterator
  before_begin ()
  {  return this->impl ().before_begin ();  }

  const_iterator
  before_begin () const
  {  return this->impl ().before_begin ();  }

  const_iterator
  cbefore_begin () const
  {  return this->impl ().cbefore_begin ();  }

  iterator
  before_end ()
  {  return this->impl ().before_end ();  }

  const_iterator
  before_end () const
  {  return this->impl ().before_end ();  }

  const_iterator
  cbefore_end () const
  {  return this->impl ().cbefore_end ();  }


  void
  pop_front ()
  {
    this->valid ().erase (this->extract_key (this->impl ().begin ()));
    this->impl  ().pop_front ();
    this->validate ();
  }


  void
  erase_after (const_iterator after)
  {
    const_iterator  position = after;
    this->valid ().erase (this->extract_key (++position));
    this->impl  ().erase_after (after);
    this->validate ();
  }

  void
  erase_after (const_iterator after, const_iterator last)
  {
    // Order is not required to be the same; emulate.
    for (const_iterator scan = after; ++scan != last;)
      this->valid ().erase (this->extract_key (scan));

    this->impl ().erase_after (after, last);
    this->validate ();
  }


  void
  relink_after (const_iterator to_after, const_iterator from_after)
  {
    this->impl ().relink_after (to_after, from_after);
    this->validate ();
  }

  void
  reverse ()
  {
    this->impl ().reverse ();
    this->validate ();
  }

  void
  sort ()
  {
    this->impl ().sort ();
    this->validate ();
  }
};


// Exists for set/map specialization.
template <typename implementation_type, typename = void>
class test_container_interface;


template <typename Implementation>
class test_table : public test_container_interface <Implementation>
{
  typedef  test_wrapper <Implementation>  wrapper_type;

protected:

  typedef  typename wrapper_type::robust_type            robust_type;

public:

  typedef  typename wrapper_type::implementation_type    implementation_type;

  typedef  typename implementation_type::value_type      value_type;
  typedef  typename implementation_type::hasher          hasher;
  typedef  typename implementation_type::key_equal       key_equal;
  typedef  typename implementation_type::allocator_type  allocator_type;


protected:

  named_object <implementation_type>  _tested;
  robust_type                         _valid;


public:

  explicit
  test_table (size_t                num_buckets = 0,
              const hasher&         hash        = hasher (),
              const key_equal&      equal       = key_equal (),
              const allocator_type& allocator   = allocator_type ())
    : _tested (num_buckets, hash, equal, allocator),
      _valid  (num_buckets, hash, equal)
  {
    validate ();
  }

  test_table (const test_table& that)
    : _tested (that.impl  ()),
      _valid  (that.valid ())
  { }

# if MCT_CXX0X_SUPPORTED

  test_table (test_table&& that)
    : _tested (std::move (that.impl ())),
      _valid  (that.valid ())
  {
    that.valid ().clear ();
    that.validate ();
  }

  test_table (initializer_list <value_type> initializer,
              size_t                num_buckets = 0,
              const hasher&         hash        = hasher (),
              const key_equal&      equal       = key_equal (),
              const allocator_type& allocator   = allocator_type ())
    : _tested (initializer, num_buckets, hash, equal, allocator),
      _valid  (num_buckets, hash, equal)
  {
    // Not requiring that the robust implementation has this constructor.
    valid ().insert (initializer.begin (), initializer.end ());
    validate ();
  }

# endif  // MCT_CXX0X_SUPPORTED

  template <typename InputIterator>
  test_table (InputIterator begin, InputIterator end,
              size_t                num_buckets = 0,
              const hasher&         hash        = hasher (),
              const key_equal&      equal       = key_equal (),
              const allocator_type& allocator   = allocator_type ())
    : _tested (begin, end, num_buckets, hash, equal, allocator),
      _valid  (begin, end, num_buckets, hash, equal)
  {
    validate ();
  }

  ~test_table ()
  {
    validate ();
  }


  test_table&
  operator= (const test_table& that)
  {
    this->impl  () = that.impl  ();
    this->valid () = that.valid ();
    validate ();
    return *this;
  }

# if MCT_CXX0X_SUPPORTED

  test_table&
  operator= (test_table&& that)
  {
    this->impl  () = std::move (that.impl ());
    this->valid () = that.valid ();
    validate ();

    if (&that != this)
      that.valid ().clear ();
    that.validate ();

    return *this;
  }

  test_table&
  operator= (initializer_list <value_type> initializer)
  {
    impl  () = initializer;
    valid ().clear ();
    valid ().insert (initializer.begin (), initializer.end ());
    validate ();
    return *this;
  }

# endif  // MCT_CXX0X_SUPPORTED


  virtual  const implementation_type&
  impl () const
  {  return *_tested.object;  }


  virtual  void
  validate () const
  {
    this->impl ().validate_integrity ();
    external_validator <implementation_type>::validate (this->impl ());

    BOOST_REQUIRE_EQUAL (impl (), valid ());
  }


protected:

  virtual  implementation_type&
  impl ()
  {  return *_tested.object;  }

  virtual  robust_type&
  valid ()
  {  return _valid;  }

  virtual  const robust_type&
  valid () const
  {  return _valid;  }
};


template <typename Implementation>
inline  bool
operator== (const test_table <Implementation>& table1, const test_table <Implementation>& table2)
{
  return table1.impl () == table2.impl ();
}

template <typename Implementation1, typename Implementation2>
inline  bool
operator== (const Implementation1& table1, const test_table <Implementation2>& table2)
{
  return table1 == table2.impl ();
}

template <typename Implementation1, typename Implementation2>
inline  bool
operator== (const test_table <Implementation1>& table1, const Implementation2& table2)
{
  return table1.impl () == table2;
}


template <typename Implementation>
inline  void
swap (test_table <Implementation>& table1, test_table <Implementation>& table2)
{
  table1.swap (table2);
}


template <typename Implementation>
ostream&
operator<< (ostream& stream, const test_table <Implementation>& table)
{
  return stream << table.impl ();
}


template <typename parameters>
struct normal_parameters
{
  typedef  mpl::bool_<parameters::implementation::normal>  type;
};

template <typename parameters>
struct stable_order_parameters
{
  typedef  mpl::bool_<parameters::implementation::stable_order>  type;
};

template <typename parameters>
struct linked_parameters
{
  typedef  mpl::bool_<parameters::implementation::normal
                      && parameters::implementation::stable_order>
           type;
};

template <typename parameters>
struct forward_parameters
{
  typedef  mpl::bool_<!parameters::implementation::normal
                      && parameters::implementation::stable_order>
           type;
};

template <typename parameters>
struct robust_iterator_validation
{
  typedef  mpl::bool_<parameters::tested_type::iterator::ROBUST_PRECONDITION_CHECKING>  type;
};


#endif  // Multi-inclusion guard.


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
