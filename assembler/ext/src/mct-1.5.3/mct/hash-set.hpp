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


#ifndef MCT_HASH_SET_HPP
#define MCT_HASH_SET_HPP


#include <functional>

#include <mct/config.hpp>
#include <mct/intrusiveness.hpp>
#include <mct/impl/closed-hash-table.hpp>

#include MCT_HASH_HEADER


namespace mct
{

  namespace impl
  {
 
    struct set_tag
    { };


    template <typename Value, typename Allocator>
    struct set_bucket_traits
    {
      typedef  set_tag                                                  container_type;
      typedef  Value                                                    value_type;
      typedef  const value_type                                         assignable_value_type;
      typedef  value_type                                               key_type;
      typedef  typename Allocator::template rebind <value_type>::other  value_allocator_type;

      typedef  std::less <value_type>                                   key_compare;


      static  const key_type&
      extract_key (const value_type& value)
      {  return value;  }


      static  void
      construct_value (value_allocator_type& allocator, value_type& value,
                       const value_type& initializer)
      {  allocator.construct (&value, initializer);  }

#   if MCT_CXX0X_SUPPORTED

      static  void
      construct_value (value_allocator_type& allocator, value_type& value,
                       value_type&& initializer)
      {  allocator.construct (&value, std::forward <value_type> (initializer));  }

#   endif


      static  bool
      mapped_equal (const value_type& /* value1 */, const value_type& /* value2 */)
      {
        return true;
      }
    };


    template <typename Value, typename Allocator, bool keep_hashes>
    class set_bucket : public plain_bucket_base <set_bucket_traits <Value, Allocator>, keep_hashes>
    {
    public:

      typedef  typename Allocator::template rebind <set_bucket>::other  bucket_allocator_type;
      typedef  typename bucket_allocator_type::pointer                  bucket_pointer;
      typedef  typename bucket_allocator_type::const_pointer            const_bucket_pointer;
    };


    template <typename Value, typename Allocator, bool pointer_links, bool keep_hashes>
    class linked_set_bucket
      : public linked_bucket_base <set_bucket_traits <Value, Allocator>,
                                   pointer_links, keep_hashes>
    {
    public:

      typedef  typename Allocator::template rebind <linked_set_bucket>::other
               bucket_allocator_type;
      typedef  typename bucket_allocator_type::pointer        bucket_pointer;
      typedef  typename bucket_allocator_type::const_pointer  const_bucket_pointer;
    };


    template <typename Value, typename Allocator, bool pointer_links, bool keep_hashes>
    class forward_set_bucket
      : public forward_bucket_base <set_bucket_traits <Value, Allocator>,
                                    pointer_links, keep_hashes>
    {
    public:

      typedef  typename Allocator::template rebind <forward_set_bucket>::other
               bucket_allocator_type;
      typedef  typename bucket_allocator_type::pointer        bucket_pointer;
      typedef  typename bucket_allocator_type::const_pointer  const_bucket_pointer;
    };

  }


  template <typename type, typename = void>
  struct is_set : impl::false_type
  { };

  template <typename type>
  struct is_set <type, typename impl::enable_if <impl::is_same <typename type::_container_type,
                                                                impl::set_tag>::value>::type>
    : impl::true_type
  { };


  // FIXME: Is it possible to at least reduce the following boilerplate somehow other than
  //        with macros?  The latter would probably lead to confusing debugging.

  template <typename Value,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Value>,
            typename Equal       = std::equal_to <Value>,
            typename Allocator   = std::allocator <Value>,
            bool     keep_hashes = false>
  class closed_hash_set
    : public impl::closed_hash_table <impl::set_bucket <Value, Allocator, keep_hashes>,
                                      Hash, Equal>
  {
    typedef  impl::closed_hash_table <impl::set_bucket <Value, Allocator, keep_hashes>,
                                      Hash, Equal>
             table_type;

  public:

    typedef  typename table_type::value_type  value_type;


    explicit
    closed_hash_set (std::size_t      num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    closed_hash_set (const closed_hash_set& that)
      : table_type (that)
    { }

    closed_hash_set (const closed_hash_set& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    closed_hash_set (closed_hash_set&& that)
      : table_type (std::move (that))
    { }

    closed_hash_set (std::initializer_list <value_type> initializer,
                     std::size_t      num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif  // MCT_CXX0X_SUPPORTED

    template <typename InputIterator>
    closed_hash_set (InputIterator first, InputIterator last,
                     std::size_t      num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    closed_hash_set&
    operator= (const closed_hash_set& that)
    {
      return static_cast <closed_hash_set&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    closed_hash_set&
    operator= (closed_hash_set&& that)
    {
      return static_cast <closed_hash_set&> (table_type::operator= (std::move (that)));
    }

    closed_hash_set&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <closed_hash_set&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Value,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Value>,
            typename Equal       = std::equal_to <Value>,
            typename Allocator   = std::allocator <Value>,
            bool     keep_hashes = false>
  class linked_hash_set
    : public impl::linked_hash_table <impl::linked_set_bucket <Value, Allocator,
                                                               sizeof (void*) <= sizeof (int),
                                                               keep_hashes>,
                                      Hash, Equal>
  {
    typedef  impl::linked_hash_table <impl::linked_set_bucket <Value, Allocator,
                                                               sizeof (void*) <= sizeof (int),
                                                               keep_hashes>,
                                      Hash, Equal>
             table_type;

  public:

    typedef  typename table_type::value_type  value_type;


    explicit
    linked_hash_set (std::size_t      num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    linked_hash_set (const linked_hash_set& that)
      : table_type (that)
    { }

    linked_hash_set (const linked_hash_set& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    linked_hash_set (linked_hash_set&& that)
      : table_type (std::move (that))
    { }

    linked_hash_set (std::initializer_list <value_type> initializer,
                     std::size_t      num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif

    template <typename InputIterator>
    linked_hash_set (InputIterator first, InputIterator last,
                     std::size_t      num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    linked_hash_set&
    operator= (const linked_hash_set& that)
    {
      return static_cast <linked_hash_set&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    linked_hash_set&
    operator= (linked_hash_set&& that)
    {
      return static_cast <linked_hash_set&> (table_type::operator= (std::move (that)));
    }

    linked_hash_set&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <linked_hash_set&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Value,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Value>,
            typename Equal       = std::equal_to <Value>,
            typename Allocator   = std::allocator <Value>,
            bool     keep_hashes = false>
  class huge_linked_hash_set
    : public impl::linked_hash_table <impl::linked_set_bucket <Value, Allocator,
                                                               true, keep_hashes>,
                                      Hash, Equal>
  {
    typedef  impl::linked_hash_table <impl::linked_set_bucket <Value, Allocator,
                                                               true, keep_hashes>,
                                      Hash, Equal>
             table_type;

  public:

    typedef  typename table_type::value_type  value_type;


    explicit
    huge_linked_hash_set (std::size_t      num_buckets = 0,
                          const Hash&      hash        = Hash (),
                          const Equal&     equal       = Equal (),
                          const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    huge_linked_hash_set (const huge_linked_hash_set& that)
      : table_type (that)
    { }

    huge_linked_hash_set (const huge_linked_hash_set& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    huge_linked_hash_set (huge_linked_hash_set&& that)
      : table_type (std::move (that))
    { }

    huge_linked_hash_set (std::initializer_list <value_type> initializer,
                          std::size_t      num_buckets = 0,
                          const Hash&      hash        = Hash (),
                          const Equal&     equal       = Equal (),
                          const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif

    template <typename InputIterator>
    huge_linked_hash_set (InputIterator first, InputIterator last,
                          std::size_t      num_buckets = 0,
                          const Hash&      hash        = Hash (),
                          const Equal&     equal       = Equal (),
                          const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    huge_linked_hash_set&
    operator= (const huge_linked_hash_set& that)
    {
      return static_cast <huge_linked_hash_set&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    huge_linked_hash_set&
    operator= (huge_linked_hash_set&& that)
    {
      return static_cast <huge_linked_hash_set&> (table_type::operator= (std::move (that)));
    }

    huge_linked_hash_set&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <huge_linked_hash_set&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Value,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Value>,
            typename Equal       = std::equal_to <Value>,
            typename Allocator   = std::allocator <Value>,
            bool     keep_hashes = false>
  class forward_hash_set
    : public impl::forward_hash_table <impl::forward_set_bucket <Value, Allocator,
                                                                 sizeof (void*) <= sizeof (int),
                                                                 keep_hashes>,
                                      Hash, Equal>
  {
    typedef  impl::forward_hash_table <impl::forward_set_bucket <Value, Allocator,
                                                                 sizeof (void*) <= sizeof (int),
                                                                 keep_hashes>,
                                      Hash, Equal>
             table_type;

  public:

    typedef  typename table_type::value_type  value_type;


    explicit
    forward_hash_set (std::size_t      num_buckets = 0,
                      const Hash&      hash        = Hash (),
                      const Equal&     equal       = Equal (),
                      const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    forward_hash_set (const forward_hash_set& that)
      : table_type (that)
    { }

    forward_hash_set (const forward_hash_set& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    forward_hash_set (forward_hash_set&& that)
      : table_type (std::move (that))
    { }

    forward_hash_set (std::initializer_list <value_type> initializer,
                      std::size_t      num_buckets = 0,
                      const Hash&      hash        = Hash (),
                      const Equal&     equal       = Equal (),
                      const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif

    template <typename InputIterator>
    forward_hash_set (InputIterator first, InputIterator last,
                      std::size_t      num_buckets = 0,
                      const Hash&      hash        = Hash (),
                      const Equal&     equal       = Equal (),
                      const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    forward_hash_set&
    operator= (const forward_hash_set& that)
    {
      return static_cast <forward_hash_set&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    forward_hash_set&
    operator= (forward_hash_set&& that)
    {
      return static_cast <forward_hash_set&> (table_type::operator= (std::move (that)));
    }

    forward_hash_set&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <forward_hash_set&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Value,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Value>,
            typename Equal       = std::equal_to <Value>,
            typename Allocator   = std::allocator <Value>,
            bool     keep_hashes = false>
  class huge_forward_hash_set
    : public impl::forward_hash_table <impl::forward_set_bucket <Value, Allocator,
                                                                 true, keep_hashes>,
                                      Hash, Equal>
  {
    typedef  impl::forward_hash_table <impl::forward_set_bucket <Value, Allocator,
                                                                 true, keep_hashes>,
                                      Hash, Equal>
             table_type;

  public:

    typedef  typename table_type::value_type  value_type;


    explicit
    huge_forward_hash_set (std::size_t      num_buckets = 0,
                           const Hash&      hash        = Hash (),
                           const Equal&     equal       = Equal (),
                           const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    huge_forward_hash_set (const huge_forward_hash_set& that)
      : table_type (that)
    { }

    huge_forward_hash_set (const huge_forward_hash_set& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    huge_forward_hash_set (huge_forward_hash_set&& that)
      : table_type (std::move (that))
    { }

    huge_forward_hash_set (std::initializer_list <value_type> initializer,
                           std::size_t      num_buckets = 0,
                           const Hash&      hash        = Hash (),
                           const Equal&     equal       = Equal (),
                           const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif

    template <typename InputIterator>
    huge_forward_hash_set (InputIterator first, InputIterator last,
                           std::size_t      num_buckets = 0,
                           const Hash&      hash        = Hash (),
                           const Equal&     equal       = Equal (),
                           const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    huge_forward_hash_set&
    operator= (const huge_forward_hash_set& that)
    {
      return static_cast <huge_forward_hash_set&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    huge_forward_hash_set&
    operator= (huge_forward_hash_set&& that)
    {
      return static_cast <huge_forward_hash_set&> (table_type::operator= (std::move (that)));
    }

    huge_forward_hash_set&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <huge_forward_hash_set&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Value, typename Hash, typename Equal, typename Allocator, bool keep_hashes>
  inline  void
  swap (closed_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set1,
        closed_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set2)
  {
    set1.swap (set2);
  }

  template <typename Value, typename Hash, typename Equal, typename Allocator, bool keep_hashes>
  inline  void
  swap (linked_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set1,
        linked_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set2)
  {
    set1.swap (set2);
  }

  template <typename Value, typename Hash, typename Equal, typename Allocator, bool keep_hashes>
  inline  void
  swap (huge_linked_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set1,
        huge_linked_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set2)
  {
    set1.swap (set2);
  }

  template <typename Value, typename Hash, typename Equal, typename Allocator, bool keep_hashes>
  inline  void
  swap (forward_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set1,
        forward_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set2)
  {
    set1.swap (set2);
  }

  template <typename Value, typename Hash, typename Equal, typename Allocator, bool keep_hashes>
  inline  void
  swap (huge_forward_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set1,
        huge_forward_hash_set <Value, Hash, Equal, Allocator, keep_hashes>& set2)
  {
    set1.swap (set2);
  }

}


#endif  // Multi-inclusion guard.


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
