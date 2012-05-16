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


#define BOOST_AUTO_TEST_MAIN
#include "tests/common.hpp"
#include "tests/external-validator.hpp"

#include <cstdlib>
#include <sstream>


size_t  foo                               ::num_alive     =  0;
size_t  int_wrapper                       ::num_destroyed =  0;
int     throws_on_copy                    ::countdown     = -1;
int     throws_in_comparator              ::countdown     = -1;
int     throws_in_hasher                  ::countdown     = -1;
int     throwing_int_strict_weak_ordering ::countdown     = -1;

int                                                    leak_test_data::num_allocators = 0;
set <pair <const type_info*, pair <void*, size_t> > >  leak_test_data::allocated_chunks;
set <pair <const type_info*, void*> >                  leak_test_data::constructed_objects;


#if HAVE_BOOST_INTERPROCESS
interp_memory_chunk  interp_memory (true);
#endif


void
leak_test_data::delete_allocator ()
{
  if (!--num_allocators)
    {
      BOOST_CHECK_EQUAL (allocated_chunks.size (), 0u);
      BOOST_CHECK_EQUAL (constructed_objects.size (), 0u);

      allocated_chunks.clear ();
      constructed_objects.clear ();
    }
}


#if HAVE_BOOST_INTERPROCESS

void
validate_externally (const void* object, int type_index)
{
  struct cleanup
  {
    ~cleanup ()
    {
      interp_memory.chunk.destroy <external_validation_data>  ("to_validate");
      interp_memory.chunk.destroy <external_validation_error> ("error");
    }
  };

  {
    // Boost test framework seems to somehow handle std::system() failures in such a way
    // our cleanup code won't even run.  How very useful of them.  Make sure the shared
    // objects are not there.
    cleanup  _;
    MCT_UNUSED (_);
  }

  static  bool  calling_external_validator = false;

  BOOST_REQUIRE_MESSAGE (!calling_external_validator,
                         "external validator program is not available");

  interp_memory.chunk.construct <external_validation_data> ("to_validate") (object, type_index);

  cleanup  _;
  MCT_UNUSED (_);

  calling_external_validator = true;
  std::system ("./tests/external-validator");
  calling_external_validator = false;

  const external_validation_error* const  error
    = interp_memory.chunk.find <external_validation_error> ("error").first;

  if (error)
    BOOST_ERROR ("in external validation: " + error->message ());
}


string
generate_unique_object_name ()
{
  static  int  counter = 0;

  std::ostringstream  buffer;

  buffer << "obj" << ++counter;
  return buffer.str ();
}

#endif  // HAVE_BOOST_INTERPROCESS


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
