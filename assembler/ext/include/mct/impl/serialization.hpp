// This file is part of Miscellaneous Container Templates.
//
//             https://launchpad.net/libmct/
//
// Copyright (c) 2012 Paul Pogonyshev.
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


#ifndef MCT_IMPL_SERIALIZATION_HPP
#define MCT_IMPL_SERIALIZATION_HPP


#include <boost/aligned_storage.hpp>
#include <boost/serialization/collection_size_type.hpp>
#include <boost/serialization/item_version_type.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/version.hpp>


namespace mct
{

  namespace impl
  {

    template <typename type>
    class object_unserializer
    {
      // Given that there is Boost.Serialization library, we sure can assume
      // 'aligned_storage' structure presence.
      boost::aligned_storage <sizeof (type)>  _storage;

    public:

      template <typename Archive>
      object_unserializer (Archive& archive, unsigned file_version)
      {
        boost::serialization::load_construct_data (archive, address (), file_version);
      }

      ~object_unserializer ()
      {
        address ()->~type ();
      }

      type*
      address ()
      {  return static_cast <type*> (_storage.address ());  }
    };


    template <typename table_type, typename archive_type>
    inline  typename enable_if <archive_type::is_saving::value>::type
    do_serialize (archive_type& archive, const table_type& table, unsigned)
    {
      boost::serialization::collection_size_type     count (table.size ());
      const float                                    max_load_factor = table.max_load_factor ();
      const boost::serialization::item_version_type  item_version
        (boost::serialization::version <typename table_type::value_type>::value);

      archive << BOOST_SERIALIZATION_NVP (count)
              << BOOST_SERIALIZATION_NVP (max_load_factor)
              << BOOST_SERIALIZATION_NVP (item_version);

      for (typename table_type::const_iterator iterator = table.begin (); count-- > 0;)
        {
          boost::serialization::save_construct_data (archive, &*iterator, item_version);
          archive << boost::serialization::make_nvp ("item", *iterator++);
        }
    }


    template <typename table_type, typename archive_type>
    inline  typename enable_if <mct::is_set <table_type>::value>::type
    reset_item_address (archive_type& archive,
                        const typename table_type::value_type& stored_at,
                        const typename table_type::value_type& copied_from)
    {
      archive.reset_object_address (&stored_at, &copied_from);
    }

    template <typename table_type, typename archive_type>
    inline  typename enable_if <mct::is_map <table_type>::value>::type
    reset_item_address (archive_type& archive,
                        const typename table_type::value_type& stored_at,
                        const typename table_type::value_type& copied_from)
    {
      archive.reset_object_address (&stored_at.second, &copied_from.second);
    }

    // Logically this should be 'unserialize', but we override the same name, so that the
    // callers can use the function regardless of archive type.
    template <typename table_type, typename archive_type>
    inline  typename enable_if <archive_type::is_loading::value>::type
    do_serialize (archive_type& archive, table_type& table, unsigned file_version)
    {
      boost::serialization::collection_size_type  count;
      float                                       max_load_factor;
      boost::serialization::item_version_type     item_version;

      archive >> BOOST_SERIALIZATION_NVP (count)
              >> BOOST_SERIALIZATION_NVP (max_load_factor)
              >> BOOST_SERIALIZATION_NVP (item_version);

      // I'm not sure if a call to clear() is needed (isn't the object freshly constructed
      // at this point?), but this is what Boost itself does.  Let's assume they know
      // better, not like it's a big deal anyway.
      table.clear ();
      table.max_load_factor (max_load_factor);
      table.rehash_for_insertion (count);

#   if MCT_SELF_VALIDATION
      // Not a proper type because it doesn't matter and because I don't want to sprinkle
      // the code with more friend declarations.
      void*  buckets = table._data.buckets;
#   endif

      for (typename table_type::size_type items = 0; items < count; ++items)
        {
          object_unserializer <typename table_type::value_type>  item (archive, file_version);
          archive >> boost::serialization::make_nvp ("item", *item.address ());

          std::pair <typename table_type::iterator, bool>  result
            = table.insert (*item.address ());

#       if MCT_SELF_VALIDATION
          if (!buckets)
            buckets = table._data.buckets;
          else if (table._data.buckets != buckets)
            throw new std::logic_error ("table bucket array reallocated during unserialization");
#       endif

          reset_item_address <table_type, archive_type> (archive, *result.first, *item.address ());
        }

      // Strictly speaking not precondition, but should be good enough.
      MCT_PRECONDITION (table.size () == count,
                        "serialized representation contained items with equal keys");
    }

  }

}


#endif  // Multi-inclusion guard.


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
