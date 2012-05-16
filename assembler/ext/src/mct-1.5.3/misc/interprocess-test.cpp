// This test fails on GCC 4.6 with Boost 1.46 when compiled with optimization:
//
//     g++ interprocess.cpp -lboost_thread -o interprocess
//
// I don't see what's the problem, but it must be in either compiler or Boost.  With GCC
// 4.5 or 4.4 or with Boost 1.48 the test passes.  Therefore, MCT disables all
// Interprocess tests if Boost 1.48 or later is not found --- this is kinda unlikely a
// problem with MCT itself, but I don't want to see failing tests either way.


#include <iostream>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>

namespace interp = boost::interprocess;

struct interp_memory_chunk
{
  interp::managed_shared_memory  chunk;

  interp_memory_chunk ()
  {
    interp::shared_memory_object::remove ("MCT_interprocess_test");
    chunk = interp::managed_shared_memory (interp::create_only, "MCT_interprocess_test", 0x10000);
  }

  ~interp_memory_chunk ()
  {
    interp::shared_memory_object::remove ("MCT_interprocess_test");
  }
};

typedef  interp::allocator <int, interp::managed_shared_memory::segment_manager>  allocator_type;

inline  void
create (allocator_type& allocator, allocator_type::value_type& at, int value)
{
  allocator.construct (allocator.address (at), value);
}

int
main ()
{
  interp_memory_chunk      memory;
  allocator_type           allocator (memory.chunk.get_segment_manager ());
  allocator_type::pointer  data = allocator.allocate (1);

  // Simply create an integer initialized to a specific value.
  create (allocator, *data, 0xdeadbeef);

  // Immediately test it; with setup mentioned above the value will be zero...
  const int  value = *data;

  if (value == 0xdeadbeef)
    {
      std::cout << "pass\n";
      return 0;
    }
  else
    {
      std::cout << "FAIL (" << std::hex << value << ")\n";
      return 1;
    }
}
