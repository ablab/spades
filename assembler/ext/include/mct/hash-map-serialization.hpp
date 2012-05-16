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


#ifndef MCT_HASH_MAP_SERIALIZATION_HPP
#define MCT_HASH_MAP_SERIALIZATION_HPP


// It's up to the includer to make sure that Boost.Serialization is accessible.  The main
// reason this doesn't go directly to 'hash-map.hpp' is to avoid any not really necessary
// dependency and leave such configuration to the user.

// For 'std::pair' serialization function.
#include <boost/serialization/utility.hpp>

#include <mct/hash-map.hpp>
#include <mct/impl/serialization.hpp>


// Per library suggestion we define our functions in 'boost::serialization' namespace.
namespace boost
{

  namespace serialization
  {

    template <typename Archive, typename Key, typename Mapped, typename Hash, typename Equal,
              typename Allocator, bool keep_hashes>
    inline  void
    serialize (Archive& archive,
               mct::closed_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map,
               unsigned file_version)
    {
      mct::impl::do_serialize (archive, map, file_version);
    }

    template <typename Archive, typename Key, typename Mapped, typename Hash, typename Equal,
              typename Allocator, bool keep_hashes>
    inline  void
    serialize (Archive& archive,
               mct::linked_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map,
               unsigned file_version)
    {
      mct::impl::do_serialize (archive, map, file_version);
    }

    template <typename Archive, typename Key, typename Mapped, typename Hash, typename Equal,
              typename Allocator, bool keep_hashes>
    inline  void
    serialize (Archive& archive,
               mct::huge_linked_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map,
               unsigned file_version)
    {
      mct::impl::do_serialize (archive, map, file_version);
    }

    template <typename Archive, typename Key, typename Mapped, typename Hash, typename Equal,
              typename Allocator, bool keep_hashes>
    inline  void
    serialize (Archive& archive,
               mct::forward_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map,
               unsigned file_version)
    {
      mct::impl::do_serialize (archive, map, file_version);
    }

    template <typename Archive, typename Key, typename Mapped, typename Hash, typename Equal,
              typename Allocator, bool keep_hashes>
    inline  void
    serialize (Archive& archive,
               mct::huge_forward_hash_map <Key, Mapped, Hash, Equal, Allocator, keep_hashes>& map,
               unsigned file_version)
    {
      mct::impl::do_serialize (archive, map, file_version);
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
