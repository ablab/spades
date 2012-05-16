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


#ifndef MCT_HASH_MAP_HPP
#define MCT_HASH_MAP_HPP


#include <mct/config.hpp>
#include <mct/intrusiveness.hpp>
#include <mct/impl/closed-hash-table.hpp>

#include MCT_HASH_HEADER


namespace mct
{

  namespace impl
  {

    struct map_tag
    { };


    template <typename Key, typename Mapped, typename Allocator>
    struct map_bucket_traits
    {
      typedef  map_tag                                                  container_type;
      typedef  Key                                                      key_type;
      typedef  Mapped                                                   mapped_type;
      typedef  std::pair <const Key, Mapped>                            value_type;
      typedef  value_type                                               assignable_value_type;
      typedef  typename Allocator::template rebind <value_type>::other  value_allocator_type;

      struct key_compare
      {
        bool
        operator() (const value_type& a, const value_type& b) const
        {  return a.first < b.first;  }
      };


      static  const key_type&
      extract_key (const value_type& value)
      {  return value.first;  }

      static  const key_type&
      extract_key (const key_type& key)
      {  return key;  }


      static  void
      construct_value (value_allocator_type& allocator, value_type& value,
                       const value_type& initializer)
      {  allocator.construct (&value, initializer);  }

      static  void
      construct_value (value_allocator_type& allocator, value_type& value, const key_type& key)
      {  allocator.construct (&value, value_type (key, mapped_type ()));  }

#   if MCT_CXX0X_SUPPORTED

      static  void
      construct_value (value_allocator_type& allocator, value_type& value,
                       value_type&& initializer)
      {
        allocator.construct (&value, std::forward <value_type> (initializer));
      }

      static  void
      construct_value (value_allocator_type& allocator, value_type& value, key_type&& key)
      {
        allocator.construct (&value, value_type (std::forward <key_type> (key), mapped_type ()));
      }

#   endif


      static  bool
      mapped_equal (const value_type& value1, const value_type& value2)
      {
        return value1.second == value2.second;
      }
    };


    template <typename Key, typename Mapped, typename Allocator, bool keep_hashes>
    class map_bucket : public plain_bucket_base <map_bucket_traits <Key, Mapped, Allocator>,
                                                 keep_hashes>
    {
    public:

      typedef  typename Allocator::template rebind <map_bucket>::other  bucket_allocator_type;
      typedef  typename bucket_allocator_type::pointer                  bucket_pointer;
      typedef  typename bucket_allocator_type::const_pointer            const_bucket_pointer;
    };


    template <typename Key, typename Mapped, typename Allocator,
              bool pointer_links, bool keep_hashes>
    class linked_map_bucket
      : public linked_bucket_base <map_bucket_traits <Key, Mapped, Allocator>,
                                   pointer_links, keep_hashes>
    {
    public:

      typedef  typename Allocator::template rebind <linked_map_bucket>::other
               bucket_allocator_type;
      typedef  typename bucket_allocator_type::pointer        bucket_pointer;
      typedef  typename bucket_allocator_type::const_pointer  const_bucket_pointer;
    };


    template <typename Key, typename Mapped, typename Allocator,
              bool pointer_links, bool keep_hashes>
    class forward_map_bucket
      : public forward_bucket_base <map_bucket_traits <Key, Mapped, Allocator>,
                                    pointer_links, keep_hashes>
    {
    public:

      typedef  typename Allocator::template rebind <forward_map_bucket>::other
               bucket_allocator_type;
      typedef  typename bucket_allocator_type::pointer        bucket_pointer;
      typedef  typename bucket_allocator_type::const_pointer  const_bucket_pointer;
    };

  }


  template <typename type, typename = void>
  struct is_map : impl::false_type
  { };

  template <typename type>
  struct is_map <type, typename impl::enable_if <impl::is_same <typename type::_container_type,
                                                                impl::map_tag>::value>::type>
    : impl::true_type
  { };


  // FIXME: Is it possible to at least reduce the following boilerplate somehow other than
  //        with macros?  The latter would probably lead to confusing debugging.

  template <typename Key,
            typename Mapped,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal       = std::equal_to <Key>,
            typename Allocator   = std::allocator <std::pair <const Key, Mapped> >,
            bool     keep_hashes = false>
  class closed_hash_map
    : public impl::closed_hash_table <impl::map_bucket <Key, Mapped, Allocator, keep_hashes>,
                                      Hash, Equal>
  {
  protected:

    typedef  impl::closed_hash_table <impl::map_bucket <Key, Mapped, Allocator, keep_hashes>,
                                      Hash, Equal>
             table_type;

    typedef  typename table_type::bucket_pointer  bucket_pointer;

  public:

    typedef  typename table_type::key_type        key_type;
    typedef  Mapped                               mapped_type;
    typedef  typename table_type::value_type      value_type;


    explicit
    closed_hash_map (size_t           num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    closed_hash_map (const closed_hash_map& that)
      : table_type (that)
    { }

    closed_hash_map (const closed_hash_map& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    closed_hash_map (closed_hash_map&& that)
      : table_type (std::move (that))
    { }

    closed_hash_map (std::initializer_list <value_type> initializer,
                     std::size_t      num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif  // MCT_CXX0X_SUPPORTED

    template <typename InputIterator>
    closed_hash_map (InputIterator first, InputIterator last,
                     size_t           num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    mapped_type&
    operator[] (const key_type& key)
    {
      return this->template lookup_or_insert <const key_type&> (key).first->value ().second;
    }

# if MCT_CXX0X_SUPPORTED

    mapped_type&
    operator[] (key_type&& key)
    {
      return (this->template lookup_or_insert <key_type&&> (std::forward <key_type> (key))
              .first->value ().second);
    }

# endif

    mapped_type&
    at (const key_type& key)
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("closed_hash_map::at");

      return bucket->value ().second;
    }

    const mapped_type&
    at (const key_type& key) const
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("closed_hash_map::at");

      return bucket->value ().second;
    }


    closed_hash_map&
    operator= (const closed_hash_map& that)
    {
      return static_cast <closed_hash_map&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    closed_hash_map&
    operator= (closed_hash_map&& that)
    {
      return static_cast <closed_hash_map&> (table_type::operator= (std::move (that)));
    }

    closed_hash_map&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <closed_hash_map&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Key,
            typename Mapped,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal       = std::equal_to <Key>,
            typename Allocator   = std::allocator <std::pair <const Key, Mapped> >,
            bool     keep_hashes = false>
  class linked_hash_map
    : public impl::linked_hash_table <impl::linked_map_bucket <Key, Mapped, Allocator,
                                                               sizeof (void*) <= sizeof (int),
                                                               keep_hashes>,
                                      Hash, Equal>
  {
  protected:

    typedef  impl::linked_hash_table <impl::linked_map_bucket <Key, Mapped, Allocator,
                                                               sizeof (void*) <= sizeof (int),
                                                               keep_hashes>,
                                      Hash, Equal>
             table_type;

    typedef  typename table_type::bucket_pointer  bucket_pointer;

  public:

    typedef  typename table_type::key_type        key_type;
    typedef  Mapped                               mapped_type;
    typedef  typename table_type::value_type      value_type;


    explicit
    linked_hash_map (size_t           num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    linked_hash_map (const linked_hash_map& that)
      : table_type (that)
    { }

    linked_hash_map (const linked_hash_map& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    linked_hash_map (linked_hash_map&& that)
      : table_type (std::move (that))
    { }

    linked_hash_map (std::initializer_list <value_type> initializer,
                     std::size_t      num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif  // MCT_CXX0X_SUPPORTED

    template <typename InputIterator>
    linked_hash_map (InputIterator first, InputIterator last,
                     size_t           num_buckets = 0,
                     const Hash&      hash        = Hash (),
                     const Equal&     equal       = Equal (),
                     const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    mapped_type&
    operator[] (const key_type& key)
    {
      return (this->template lookup_or_insert <const key_type&> (key, this->_data.end ())
              .first->value ().second);
    }

# if MCT_CXX0X_SUPPORTED

    mapped_type&
    operator[] (key_type&& key)
    {
      return (this->template lookup_or_insert <key_type&&> (std::forward <key_type> (key),
                                                            this->_data.end ())
              .first->value ().second);
    }

# endif


    mapped_type&
    at (const key_type& key)
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("linked_hash_map::at");

      return bucket->value ().second;
    }

    const mapped_type&
    at (const key_type& key) const
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("linked_hash_map::at");

      return bucket->value ().second;
    }


    linked_hash_map&
    operator= (const linked_hash_map& that)
    {
      return static_cast <linked_hash_map&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    linked_hash_map&
    operator= (linked_hash_map&& that)
    {
      return static_cast <linked_hash_map&> (table_type::operator= (std::move (that)));
    }

    linked_hash_map&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <linked_hash_map&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Key,
            typename Mapped,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal       = std::equal_to <Key>,
            typename Allocator   = std::allocator <std::pair <const Key, Mapped> >,
            bool     keep_hashes = false>
  class huge_linked_hash_map
    : public impl::linked_hash_table <impl::linked_map_bucket <Key, Mapped, Allocator,
                                                               true, keep_hashes>,
                                      Hash, Equal>
  {
  protected:

    typedef  impl::linked_hash_table <impl::linked_map_bucket <Key, Mapped, Allocator,
                                                               true, keep_hashes>,
                                      Hash, Equal>
             table_type;

    typedef  typename table_type::bucket_pointer  bucket_pointer;

  public:

    typedef  typename table_type::key_type        key_type;
    typedef  Mapped                               mapped_type;
    typedef  typename table_type::value_type      value_type;


    explicit
    huge_linked_hash_map (size_t           num_buckets = 0,
                          const Hash&      hash        = Hash (),
                          const Equal&     equal       = Equal (),
                          const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    huge_linked_hash_map (const huge_linked_hash_map& that)
      : table_type (that)
    { }

    huge_linked_hash_map (const huge_linked_hash_map& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    huge_linked_hash_map (huge_linked_hash_map&& that)
      : table_type (std::move (that))
    { }

    huge_linked_hash_map (std::initializer_list <value_type> initializer,
                          std::size_t      num_buckets = 0,
                          const Hash&      hash        = Hash (),
                          const Equal&     equal       = Equal (),
                          const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif  // MCT_CXX0X_SUPPORTED

    template <typename InputIterator>
    huge_linked_hash_map (InputIterator first, InputIterator last,
                          size_t           num_buckets = 0,
                          const Hash&      hash        = Hash (),
                          const Equal&     equal       = Equal (),
                          const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    mapped_type&
    operator[] (const key_type& key)
    {
      return (this->template lookup_or_insert <const key_type&> (key, this->_data.end ())
              .first->value ().second);
    }

# if MCT_CXX0X_SUPPORTED

    mapped_type&
    operator[] (key_type&& key)
    {
      return (this->template lookup_or_insert <key_type&&> (std::forward <key_type> (key),
                                                            this->_data.end ())
              .first->value ().second);
    }

# endif


    mapped_type&
    at (const key_type& key)
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("huge_linked_hash_map::at");

      return bucket->value ().second;
    }

    const mapped_type&
    at (const key_type& key) const
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("huge_linked_hash_map::at");

      return bucket->value ().second;
    }


    huge_linked_hash_map&
    operator= (const huge_linked_hash_map& that)
    {
      return static_cast <huge_linked_hash_map&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    huge_linked_hash_map&
    operator= (huge_linked_hash_map&& that)
    {
      return static_cast <huge_linked_hash_map&> (table_type::operator= (std::move (that)));
    }

    huge_linked_hash_map&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <huge_linked_hash_map&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Key,
            typename Mapped,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal       = std::equal_to <Key>,
            typename Allocator   = std::allocator <std::pair <const Key, Mapped> >,
            bool     keep_hashes = false>
  class forward_hash_map
    : public impl::forward_hash_table <impl::forward_map_bucket <Key, Mapped, Allocator,
                                                                 sizeof (void*) <= sizeof (int),
                                                                 keep_hashes>,
                                       Hash, Equal>
  {
  protected:

    typedef  impl::forward_hash_table <impl::forward_map_bucket <Key, Mapped, Allocator,
                                                                 sizeof (void*) <= sizeof (int),
                                                                 keep_hashes>,
                                       Hash, Equal>
             table_type;

    typedef  typename table_type::bucket_pointer  bucket_pointer;

  public:

    typedef  typename table_type::key_type        key_type;
    typedef  Mapped                               mapped_type;
    typedef  typename table_type::value_type      value_type;


    explicit
    forward_hash_map (size_t           num_buckets = 0,
                      const Hash&      hash        = Hash (),
                      const Equal&     equal       = Equal (),
                      const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    forward_hash_map (const forward_hash_map& that)
      : table_type (that)
    { }

    forward_hash_map (const forward_hash_map& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    forward_hash_map (forward_hash_map&& that)
      : table_type (std::move (that))
    { }

    forward_hash_map (std::initializer_list <value_type> initializer,
                      std::size_t      num_buckets = 0,
                      const Hash&      hash        = Hash (),
                      const Equal&     equal       = Equal (),
                      const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif  // MCT_CXX0X_SUPPORTED

    template <typename InputIterator>
    forward_hash_map (InputIterator first, InputIterator last,
                      size_t           num_buckets = 0,
                      const Hash&      hash        = Hash (),
                      const Equal&     equal       = Equal (),
                      const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    mapped_type&
    operator[] (const key_type& key)
    {
      return (this->template lookup_or_insert <const key_type&> (key, this->_data.back)
              .first->value ().second);
    }

# if MCT_CXX0X_SUPPORTED

    mapped_type&
    operator[] (key_type&& key)
    {
      return (this->template lookup_or_insert <key_type&&> (std::forward <key_type> (key),
                                                            this->_data.back)
              .first->value ().second);
    }

# endif


    mapped_type&
    at (const key_type& key)
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("forward_hash_map::at");

      return bucket->value ().second;
    }

    const mapped_type&
    at (const key_type& key) const
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("forward_hash_map::at");

      return bucket->value ().second;
    }


    forward_hash_map&
    operator= (const forward_hash_map& that)
    {
      return static_cast <forward_hash_map&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    forward_hash_map&
    operator= (forward_hash_map&& that)
    {
      return static_cast <forward_hash_map&> (table_type::operator= (std::move (that)));
    }

    forward_hash_map&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <forward_hash_map&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Key,
            typename Mapped,
            typename Hash        = MCT_HASH_NAMESPACE::hash <Key>,
            typename Equal       = std::equal_to <Key>,
            typename Allocator   = std::allocator <std::pair <const Key, Mapped> >,
            bool     keep_hashes = false>
  class huge_forward_hash_map
    : public impl::forward_hash_table <impl::forward_map_bucket <Key, Mapped, Allocator,
                                                                 true, keep_hashes>,
                                       Hash, Equal>
  {
  protected:

    typedef  impl::forward_hash_table <impl::forward_map_bucket <Key, Mapped, Allocator,
                                                                 true, keep_hashes>,
                                       Hash, Equal>
             table_type;

    typedef  typename table_type::bucket_pointer  bucket_pointer;

  public:

    typedef  typename table_type::key_type        key_type;
    typedef  Mapped                               mapped_type;
    typedef  typename table_type::value_type      value_type;


    explicit
    huge_forward_hash_map (size_t           num_buckets = 0,
                           const Hash&      hash        = Hash (),
                           const Equal&     equal       = Equal (),
                           const Allocator& allocator   = Allocator ())
      : table_type (num_buckets, hash, equal, allocator)
    { }

    huge_forward_hash_map (const huge_forward_hash_map& that)
      : table_type (that)
    { }

    huge_forward_hash_map (const huge_forward_hash_map& that, const Allocator& allocator)
      : table_type (that, allocator)
    { }

# if MCT_CXX0X_SUPPORTED

    huge_forward_hash_map (huge_forward_hash_map&& that)
      : table_type (std::move (that))
    { }

    huge_forward_hash_map (std::initializer_list <value_type> initializer,
                           std::size_t      num_buckets = 0,
                           const Hash&      hash        = Hash (),
                           const Equal&     equal       = Equal (),
                           const Allocator& allocator   = Allocator ())
      : table_type (initializer, num_buckets, hash, equal, allocator)
    { }

# endif  // MCT_CXX0X_SUPPORTED

    template <typename InputIterator>
    huge_forward_hash_map (InputIterator first, InputIterator last,
                           size_t           num_buckets = 0,
                           const Hash&      hash        = Hash (),
                           const Equal&     equal       = Equal (),
                           const Allocator& allocator   = Allocator ())
      : table_type (first, last, num_buckets, hash, equal, allocator)
    { }


    mapped_type&
    operator[] (const key_type& key)
    {
      return (this->template lookup_or_insert <const key_type&> (key, this->_data.back)
              .first->value ().second);
    }

# if MCT_CXX0X_SUPPORTED

    mapped_type&
    operator[] (key_type&& key)
    {
      return (this->template lookup_or_insert <key_type&&> (std::forward <key_type> (key),
                                                            this->_data.back)
              .first->value ().second);
    }

# endif


    mapped_type&
    at (const key_type& key)
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("huge_forward_hash_map::at");

      return bucket->value ().second;
    }

    const mapped_type&
    at (const key_type& key) const
    {
      const bucket_pointer  bucket (this->lookup (key));
      if (MCT_OPTIMIZATION_UNLIKELY (bucket == this->_data.end ()))
        throw std::out_of_range ("huge_forward_hash_map::at");

      return bucket->value ().second;
    }


    huge_forward_hash_map&
    operator= (const huge_forward_hash_map& that)
    {
      return static_cast <huge_forward_hash_map&> (table_type::operator= (that));
    }

# if MCT_CXX0X_SUPPORTED

    huge_forward_hash_map&
    operator= (huge_forward_hash_map&& that)
    {
      return static_cast <huge_forward_hash_map&> (table_type::operator= (std::move (that)));
    }

    huge_forward_hash_map&
    operator= (std::initializer_list <value_type> initializer)
    {
      return static_cast <huge_forward_hash_map&> (table_type::operator= (initializer));
    }

# endif  // MCT_CXX0X_SUPPORTED
  };


  template <typename Key, typename Mapped, typename Hash, typename Equal, typename Allocator,
            bool keep_hashes>
  inline  void
  swap (closed_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map1,
        closed_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map2)
  {
    map1.swap (map2);
  }

  template <typename Key, typename Mapped, typename Hash, typename Equal, typename Allocator,
            bool keep_hashes>
  inline  void
  swap (linked_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map1,
        linked_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map2)
  {
    map1.swap (map2);
  }

  template <typename Key, typename Mapped, typename Hash, typename Equal, typename Allocator,
            bool keep_hashes>
  inline  void
  swap (huge_linked_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map1,
        huge_linked_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map2)
  {
    map1.swap (map2);
  }

  template <typename Key, typename Mapped, typename Hash, typename Equal, typename Allocator,
            bool keep_hashes>
  inline  void
  swap (forward_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map1,
        forward_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map2)
  {
    map1.swap (map2);
  }

  template <typename Key, typename Mapped, typename Hash, typename Equal, typename Allocator,
            bool keep_hashes>
  inline  void
  swap (huge_forward_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map1,
        huge_forward_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map2)
  {
    map1.swap (map2);
  }

}


#endif  // Multi-inclusion guard.


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
