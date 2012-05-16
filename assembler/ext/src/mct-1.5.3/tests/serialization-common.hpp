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


#ifndef MCT_TESTS_SERIALIZATION_COMMON_HPP
#define MCT_TESTS_SERIALIZATION_COMMON_HPP

#include <sstream>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>


// This header is supposed to be included directly in '.cpp' files.
namespace
{

  template <typename type>
  void
  test_serialization (const type& object)
  {
    // Hardly worth it to make yet another template for this.

    {
      std::stringstream  buffer;

      {
        boost::archive::xml_oarchive  out (buffer);
        out << boost::serialization::make_nvp ("object", object);
      }

      {
        boost::archive::xml_iarchive  in (buffer);
        type                          object2;

        in >> boost::serialization::make_nvp ("object", object2);
        BOOST_CHECK_EQUAL (object, object2);
      }
    }

    {
      std::stringstream  buffer;

      {
        boost::archive::text_oarchive  out (buffer);
        out << boost::serialization::make_nvp ("object", object);
      }

      {
        boost::archive::text_iarchive  in (buffer);
        type                           object2;

        in >> boost::serialization::make_nvp ("object", object2);
        BOOST_CHECK_EQUAL (object, object2);
      }
    }

    {
      std::stringstream  buffer;

      {
        boost::archive::binary_oarchive  out (buffer);
        out << boost::serialization::make_nvp ("object", object);
      }

      {
        boost::archive::binary_iarchive  in (buffer);
        type                             object2;

        in >> boost::serialization::make_nvp ("object", object2);
        BOOST_CHECK_EQUAL (object, object2);
      }
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
