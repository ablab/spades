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


#include "tests/map-common.hpp"
#include "tests/set-common.hpp"
#include "tests/common.hpp"
#include "tests/external-validator.hpp"

#include <sstream>


namespace
{

  int
  error_state (const string& message)
  {
    interp_memory.chunk.construct <external_validation_error> ("error") (message);

    // It seems any nonzero return state makes Boost.Test mad.  We rely on the error
    // message alone for this reason.
    return 0;
  }


  class object_validator
  {
    const external_validation_data  _data;
    // Referencing an external class rather than encapsulating it because Boost's for_each
    // creates copies of the functor rather than reusing the same object.
    bool&                           _validated;


  public:

    object_validator (const external_validation_data& data, bool& validated)
      : _data      (data),
        _validated (validated)
    { }


    template <typename type>
    void
    operator() (const type&)
    {
      if (is_set <type>::value)
        {
          typedef  typename mpl::find <test_set_types_interp_allocator, type>::type  position;
          typedef  mpl::begin <test_set_types_interp_allocator>::type                begin;

          if (_data.second
              != TEST_SET_TYPE_MAGIC_NUMBER + mpl::distance <begin, position>::type::value)
            return;
        }

      if (is_map <type>::value)
        {
          typedef  typename mpl::find <test_map_types_interp_allocator, type>::type  position;
          typedef  mpl::begin <test_map_types_interp_allocator>::type                begin;

          if (_data.second
              != TEST_MAP_TYPE_MAGIC_NUMBER + mpl::distance <begin, position>::type::value)
            return;
        }

      if (_validated)
        throw logic_error ("implementation type matches two validators");

      _validated = true;
      static_cast <const type*> (_data.first.get ())->validate_integrity ();
    }
  };

}


// FIXME: How to get rid of all this stuff?  Needed in the main program, but not really
//        here...
size_t  foo                               ::num_alive     =  0;
size_t  int_wrapper                       ::num_destroyed =  0;
int     throws_on_copy                    ::countdown     = -1;
int     throws_in_comparator              ::countdown     = -1;
int     throws_in_hasher                  ::countdown     = -1;
int     throwing_int_strict_weak_ordering ::countdown     = -1;

int                                                    leak_test_data::num_allocators = 0;
set <pair <const type_info*, pair <void*, size_t> > >  leak_test_data::allocated_chunks;
set <pair <const type_info*, void*> >                  leak_test_data::constructed_objects;


interp_memory_chunk  interp_memory (false);


void
leak_test_data::delete_allocator ()
{ }


int
main ()
{
  const external_validation_data* const  to_validate
    = interp_memory.chunk.find <external_validation_data> ("to_validate").first;

  if (!to_validate)
    return error_state ("no validation information found");

  try
    {
      bool  validated = false;

      mpl::for_each <test_set_types_interp_allocator> (object_validator (*to_validate, validated));
      mpl::for_each <test_map_types_interp_allocator> (object_validator (*to_validate, validated));

      if (!validated)
        {
          std::ostringstream  buffer;
          buffer << "type " << to_validate->second << " didn't match any validators";
          throw logic_error (buffer.str ());
        }
    }
  catch (logic_error& exception)
    {
      return error_state (exception.what ());
    }

  return 0;
}


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
