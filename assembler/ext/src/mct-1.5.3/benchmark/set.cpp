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


// This is a benchmark, we don't want assertions to slow anything down.
#ifndef NDEBUG
# define NDEBUG
#endif

#include "benchmark/config.hpp"

#include <algorithm>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>

#include <mct/hash-set.hpp>

#if HAVE_UNORDERED_SET
# include <unordered_set>
#endif

#if HAVE_TR1_UNORDERED_SET
# include <tr1/unordered_set>
#endif

#if HAVE_BOOST_UNORDERED_SET_HPP
# include <boost/unordered_set.hpp>
#endif

#if HAVE_GOOGLE_DENSE_HASH_SET
# include <google/dense_hash_set>
#endif


using namespace std;


namespace
{

  struct error_exit
  { };


  // MCT's internal "huge prime multiplication" trick seems to give it a solid benefit in
  // cases of bad hash functions.  So, let's use something not as bad for benchmarking and
  // comparing with other implementations.  You can test that this indeed disadvantages
  // MCT on average by removing 31337 below.
  struct rehash_int : MCT_HASH_NAMESPACE::hash <int>
  {
    size_t
    operator() (int value) const
    {
      return 31337 * MCT_HASH_NAMESPACE::hash <int>::operator() (value);
    }
  };


  template <size_t size>
  struct byte_buffer
  {
    unsigned char  value[size];

    void
    fill_with (size_t bytes)
    {
      for (size_t k = 0, shift = 0; k != size; ++k)
        {
          value[k] = (bytes & (0xff << shift)) >> shift;
          if ((shift += 8) >= 8 * sizeof (size_t))
            shift = 0;
        }
    }
  };

  template <size_t size>
  inline  bool
  operator== (const byte_buffer <size>& buffer1, const byte_buffer <size>& buffer2)
  {
    return equal (buffer1.value, buffer1.value + size, buffer2.value);
  }

  template <size_t size>
  struct buffer_hash
  {
    size_t
    operator() (const byte_buffer <size>& buffer) const
    {
      size_t  hash = 0;
      for (const unsigned char* first = buffer.value, * last = first + size;
           first != last; ++first)
        hash = (37 * hash) + *first;

      return hash;
    }
  };


  template <typename Type, typename Hash, typename Equal, typename Allocator, bool keep_hashes>
  class mct_closed_hash_set
    : public mct::closed_hash_set <Type, Hash, Equal, Allocator, keep_hashes>
  {
  public:

    ~mct_closed_hash_set ()
    {
      print_statistics (*this);
    }
  };

  template <typename Type, typename Hash, typename Equal, typename Allocator, bool keep_hashes>
  class mct_linked_hash_set
    : public mct::linked_hash_set <Type, Hash, Equal, Allocator, keep_hashes>
  {
  public:

    ~mct_linked_hash_set ()
    {
      print_statistics (*this);
    }
  };

  template <typename Type, typename Hash, typename Equal, typename Allocator, bool keep_hashes>
  class mct_forward_hash_set
    : public mct::forward_hash_set <Type, Hash, Equal, Allocator, keep_hashes>
  {
  public:

    ~mct_forward_hash_set ()
    {
      print_statistics (*this);
    }
  };


#if HAVE_GOOGLE_DENSE_HASH_SET

  template <typename Type, typename Hash, typename Equal>
  class google_dense_hash_set;

  template <typename Hash, typename Equal>
  class google_dense_hash_set <int, Hash, Equal> : public google::dense_hash_set <int, Hash, Equal>
  {
  public:

    google_dense_hash_set ()
    {
      // Benchmarks are written to never use negative values.
      this->set_empty_key (-1);
      this->set_deleted_key (-2);
    }
  };

  template <typename Hash, typename Equal>
  class google_dense_hash_set <byte_buffer <100>, Hash, Equal>
    : public google::dense_hash_set <byte_buffer <100>, Hash, Equal>
  {
  public:

    google_dense_hash_set ()
    {
      // Benchmarks are written to never use buffers filled with the same bytes.
      byte_buffer <100>  buffer;

      buffer.fill_with (0);
      this->set_empty_key (buffer);

      buffer.fill_with (static_cast <size_t> (-1));
      this->set_deleted_key (buffer);
    }
  };

#endif


  template <bool keep_hashes>
  struct mct_closed
  {
    static  const bool  normal = true;

    template <typename Type, typename Hash, typename Equal>
    struct resolve
    {
      typedef  mct_closed_hash_set <Type, Hash, Equal, allocator <Type>, keep_hashes>  set_type;
    };
  };


  template <bool keep_hashes>
  struct mct_linked
  {
    static  const bool  normal = true;

    template <typename Type, typename Hash, typename Equal>
    struct resolve
    {
      typedef  mct_linked_hash_set <Type, Hash, Equal, allocator <Type>, keep_hashes>  set_type;
    };
  };


  template <bool keep_hashes>
  struct mct_forward
  {
    static  const bool  normal = false;

    template <typename Type, typename Hash, typename Equal>
    struct resolve
    {
      typedef  mct_forward_hash_set <Type, Hash, Equal, allocator <Type>, keep_hashes>  set_type;
    };
  };


#if HAVE_UNORDERED_SET

  struct std_unordered
  {
    static  const bool  normal = true;

    template <typename Type, typename Hash, typename Equal>
    struct resolve
    {
      typedef  std::unordered_set <Type, Hash, Equal>  set_type;
    };
  };

#endif


#if HAVE_TR1_UNORDERED_SET

  struct std_tr1_unordered
  {
    static  const bool  normal = true;

    template <typename Type, typename Hash, typename Equal>
    struct resolve
    {
      typedef  std::tr1::unordered_set <Type, Hash, Equal>  set_type;
    };
  };

#endif


#if HAVE_BOOST_UNORDERED_SET_HPP

  struct boost_unordered
  {
    static  const bool  normal = true;

    template <typename Type, typename Hash, typename Equal>
    struct resolve
    {
      typedef  boost::unordered_set <Type, Hash, Equal>  set_type;
    };
  };

#endif


#if HAVE_GOOGLE_DENSE_HASH_SET

  struct google_dense
  {
    static  const bool  normal = true;

    template <typename Type, typename Hash, typename Equal>
    struct resolve
    {
      typedef  google_dense_hash_set <Type, Hash, Equal>  set_type;
    };
  };

#endif


  template <size_t approximate_size>
  struct insert_erase_int_benchmark
  {
    template <typename Implementation>
    struct applicable_to
    {
      static  const bool  value = Implementation::normal;
    };

    typedef  int              Type;
    typedef  rehash_int       Hash;
    typedef  equal_to <Type>  Equal;


    template <typename set_type>
    void
    prepare (set_type& set)
    {
      // Level playing field a bit.
      set.max_load_factor (0.7f);
      set.rehash (approximate_size);
    }

    template <typename set_type>
    size_t
    do_run (set_type& set, size_t loop_size)
    {
      size_t  size_sum = 0;

      for (size_t k = 0; k < loop_size; ++k)
        {
          // This is tested to give various intermediate sizes up to 'approximate_size' or
          // a little smaller.  Surely biased, but should be good enough for a benchmark.
          const int  number = static_cast <int> ((k * (1 + 2 * k)) % approximate_size);
          if (k % (7 * approximate_size) <= 4 * approximate_size)
            set.insert (number);
          else
            set.erase (number);

#       if MCT_ENABLE_DEBUGGING
          size_sum += set.size ();
#       endif
        }

      return size_sum;
    }
  };


  template <size_t approximate_size>
  struct lookup_100_bytes_benchmark
  {
    template <typename Implementation>
    struct applicable_to
    {
      static  const bool  value = true;
    };

    typedef  byte_buffer <100>  Type;
    typedef  buffer_hash <100>  Hash;
    typedef  equal_to <Type>    Equal;


    vector <Type>  buffers;


    template <typename set_type>
    void
    prepare (set_type& set)
    {
      // Level playing field a bit.
      set.max_load_factor (0.7f);

      buffers.reserve (2 * approximate_size);

      Type  buffer;
      for (size_t k = 0; k != approximate_size; ++k)
        {
          buffer.fill_with (1337 + 47 * k);
          buffers.push_back (buffer);
          set.insert (buffer);

          buffer.fill_with (42 + 3 * k);
          buffers.push_back (buffer);
        }
    }

    template <typename set_type>
    size_t
    do_run (set_type& set, size_t loop_size)
    {
      size_t  check_sum = 0;

      for (size_t k = 0, num_buffers = buffers.size (), index = 0; k < loop_size; ++k)
        {
          typename set_type::const_iterator  entry = set.find (buffers[index]);
          MCT_UNUSED (entry);

#       if MCT_ENABLE_DEBUGGING
          if (entry != set.end ())
            check_sum += k;
#       endif

          if (++index == num_buffers)
            index = 0;
        }

      return check_sum;
    }
  };


  struct benchmark_loop
  {
    size_t           loop_size;
    size_t           num_tries;
# if MCT_ENABLE_DEBUGGING
    size_t           check_sum;
# endif
    vector <double>  used_times;


    explicit
    benchmark_loop (size_t loop_size, size_t num_tries)
      : loop_size (loop_size),
        num_tries (num_tries)
    { }

    template <typename Implementation, typename Benchmark>
    void
    run ()
    {
      for (size_t pass = 0; pass < num_tries; ++pass)
        {
          typedef  typename Implementation::template resolve <typename Benchmark::Type,
                                                              typename Benchmark::Hash,
                                                              typename Benchmark::Equal>::set_type
                   set_type;

          set_type   set;
          Benchmark  benchmark;

          benchmark.prepare (set);

          clock_t  start     = clock ();
          size_t   check_sum = benchmark.do_run (set, loop_size);
          clock_t  finish    = clock ();
          double   used_time = static_cast <double> (finish - start) / CLOCKS_PER_SEC;

          MCT_UNUSED (check_sum);

#       if MCT_ENABLE_DEBUGGING
          if (pass > 0 && check_sum != this->check_sum)
            {
              cerr << "\n*** error ***: check sum differs between passes: "
                   << this->check_sum << " and " << check_sum << "\n";
              throw error_exit ();
            }

          this->check_sum = check_sum;
#       endif

          used_times.push_back (used_time);

          if (num_tries > 1)
            cout << used_time << (pass + 1 < num_tries ? ' ' : '\n');
        }
    }

    double
    used_time () const
    {
      return *min_element (used_times.begin (), used_times.end ());
    }
  };


  template <typename Implementation, typename Benchmark,
            bool applicable = Benchmark::template applicable_to <Implementation>::value>
  struct benchmark_runner
  {
    static  void
    run (benchmark_loop& loop, const string& implementation, const string& /* benchmark_name */,
         bool internal)
    {
      if (internal)
        cout << "Benchmarking implementation '" << implementation << "'\n";

      loop.run <Implementation, Benchmark> ();

      cout << "Used time: " << loop.used_time () << " s\n";
#   if MCT_ENABLE_DEBUGGING
      cout << "Check sum: " << loop.check_sum << "\n";
#   endif
      if (internal)
        cout << "\n";
    }
  };

  template <typename Implementation, typename Benchmark>
  struct benchmark_runner <Implementation, Benchmark, false>
  {
    static  void
    run (benchmark_loop& /* loop */, const string& implementation, const string& benchmark_name,
         bool internal)
    {
      if (!internal)
        {
          cerr << "Benchmark '" << benchmark_name << "' is not applicable to implementation '"
               << implementation << "'\n";
          throw error_exit ();
        }
    }
  };


  template <typename Set>
  void
  print_statistics (const Set& set)
  {
    MCT_UNUSED (set);

    // This can be enabled with SCons command line; see 'variables.template'.
# if MCT_ENABLE_DEBUGGING
    typename Set::statistics  stats (set.collect_statistics ());

    cerr << "MCT closed_hash_table statistics:\n";
    cerr << "  number of buckets:   " << set.bucket_count () << "\n";
    cerr << "  debris ratio:        " << stats.debris_ratio << "\n";
    cerr << "  present item lookup: " << stats.avg_present_lookup
         << " (max. " << stats.max_present_lookup << ")\n";
    cerr << "  absent item lookup:  " << stats.avg_absent_lookup
         << " (max. " << stats.max_absent_lookup << ")\n";
    cerr << "  used memory:         " << set.used_memory () / 1024.0 << " KB\n";
# endif
  }


  template <typename Implementation>
  void
  run_benchmark (benchmark_loop& loop, const string& implementation, const string& benchmark_name,
                 bool internal)
  {
    if (benchmark_name == "ie-int-small")
      {
        benchmark_runner <Implementation, insert_erase_int_benchmark <100> >::run
          (loop, implementation, benchmark_name, internal);
      }
    else if (benchmark_name == "ie-int-large")
      {
        benchmark_runner <Implementation, insert_erase_int_benchmark <10000> >::run
          (loop, implementation, benchmark_name, internal);
      }
    else if (benchmark_name == "ie-int-huge")
      {
        benchmark_runner <Implementation, insert_erase_int_benchmark <1000000> >::run
          (loop, implementation, benchmark_name, internal);
      }
    else if (benchmark_name == "lookup-100-small")
      {
        benchmark_runner <Implementation, lookup_100_bytes_benchmark <100> >::run
          (loop, implementation, benchmark_name, internal);
      }
    else if (benchmark_name == "lookup-100-large")
      {
        benchmark_runner <Implementation, lookup_100_bytes_benchmark <10000> >::run
          (loop, implementation, benchmark_name, internal);
      }
    else if (benchmark_name == "lookup-100-huge")
      {
        benchmark_runner <Implementation, lookup_100_bytes_benchmark <1000000> >::run
          (loop, implementation, benchmark_name, internal);
      }
    else
      {
        cerr << "Unknown benchmark '" << benchmark_name << "'\n";
        throw error_exit ();
      }
  }

  void
  run_benchmark (const string& implementation, const string& benchmark_name,
                 size_t loop_size, size_t num_tries, bool internal = false)
  {
    if (implementation == "all")
      {
        run_benchmark ("mct",         benchmark_name, loop_size, num_tries, true);
        run_benchmark ("mct-kh",      benchmark_name, loop_size, num_tries, true);
        run_benchmark ("mct-link",    benchmark_name, loop_size, num_tries, true);
        run_benchmark ("mct-link-kh", benchmark_name, loop_size, num_tries, true);
        run_benchmark ("mct-fwd",     benchmark_name, loop_size, num_tries, true);
        run_benchmark ("mct-fwd-kh",  benchmark_name, loop_size, num_tries, true);
#     if HAVE_UNORDERED_SET
        run_benchmark ("std", benchmark_name, loop_size, num_tries, true);
#     endif
#     if HAVE_TR1_UNORDERED_SET
        run_benchmark ("tr1", benchmark_name, loop_size, num_tries, true);
#     endif
#     if HAVE_BOOST_UNORDERED_SET_HPP
        run_benchmark ("boost", benchmark_name, loop_size, num_tries, true);
#     endif
#     if HAVE_GOOGLE_DENSE_HASH_SET
        run_benchmark ("google", benchmark_name, loop_size, num_tries, true);
#     endif
        return;
      }

    benchmark_loop  loop (loop_size, num_tries);

    if (implementation == "mct")
      run_benchmark <mct_closed <false> > (loop, implementation, benchmark_name, internal);
    else if (implementation == "mct-kh")
      run_benchmark <mct_closed <true> > (loop, implementation, benchmark_name, internal);
    else if (implementation == "mct-link")
      run_benchmark <mct_linked <false> > (loop, implementation, benchmark_name, internal);
    else if (implementation == "mct-link-kh")
      run_benchmark <mct_linked <true> > (loop, implementation, benchmark_name, internal);
    else if (implementation == "mct-fwd")
      run_benchmark <mct_forward <false> > (loop, implementation, benchmark_name, internal);
    else if (implementation == "mct-fwd-kh")
      run_benchmark <mct_forward <true> > (loop, implementation, benchmark_name, internal);

# if HAVE_UNORDERED_SET
    else if (implementation == "std")
      run_benchmark <std_unordered> (loop, implementation, benchmark_name, internal);
# endif

# if HAVE_TR1_UNORDERED_SET
    else if (implementation == "tr1")
      run_benchmark <std_tr1_unordered> (loop, implementation, benchmark_name, internal);
# endif

# if HAVE_BOOST_UNORDERED_SET_HPP
    else if (implementation == "boost")
      run_benchmark <boost_unordered> (loop, implementation, benchmark_name, internal);
# endif

# if HAVE_GOOGLE_DENSE_HASH_SET
    else if (implementation == "google")
      run_benchmark <google_dense> (loop, implementation, benchmark_name, internal);
# endif

    else
      {
        cerr << "Unsupported set implementation '" << implementation << "'\n";
        throw error_exit ();
      }
  }

  void
  do_parse (const string& input, size_t& result)
  {
    istringstream  parser (input);
    parser >> result;

    if (parser.fail ())
      {
        cerr << "Cannot parse '" << input << "' as integer\n";
        throw error_exit ();
      }
  }

}


int
main (int argc, char** argv)
{
  size_t  loop_size = 100000000;
  size_t  num_tries = 1;

  if (argc < 3 || 5 < argc)
    {
      cerr << "Usage: " << argv[0] << " IMPLEMENTATION BENCHMARK-NAME [LOOP-SIZE [NUM-TRIES]]\n";

      cerr << "\nSupported implementations are:\n";
      cerr << "  mct         -- mct::closed_hash_set <...>\n";
      cerr << "  mct-kh      -- mct::closed_hash_set <..., true>\n";
      cerr << "  mct-link    -- mct::linked_hash_set <...>\n";
      cerr << "  mct-link-kh -- mct::linked_hash_set <..., true>\n";
      cerr << "  mct-fwd     -- mct::forward_hash_set <...>       (not for all benchmarks)\n";
      cerr << "  mct-fwd-kh  -- mct::forward_hash_set <..., true> (not for all benchmarks)\n";
#   if HAVE_UNORDERED_SET
      cerr << "  std         -- std::unordered_set <...>\n";
#   endif
#   if HAVE_TR1_UNORDERED_SET
      cerr << "  tr1         -- std::tr1::unordered_set <...>\n";
#   endif
#   if HAVE_BOOST_UNORDERED_SET_HPP
      cerr << "  boost       -- boost::unordered_set <...>\n";
#   endif
#   if HAVE_GOOGLE_DENSE_HASH_SET
      cerr << "  google      -- google::dense_hash_set <...>\n";
#   endif
      cerr << "  all         -- all of the above in that order\n";

      cerr << "\nKnown benchmark names:\n";
      cerr << "  ie-int-(small|large|huge) -- insert/erase in a set of approx. 100, 10000\n"
           << "                               or 1000000 integers correspondingly\n";
      cerr << "  lookup-100-(small|large|huge) -- create a set of 100, 10000, 1000000 elements\n"
           << "                                   100 bytes each and just lookup\n";

      cerr << "\nLOOP-SIZE is optional and defaults to " << loop_size << ".\n";

      cerr << "\nNUM-TRIES is optional and defaults to " << num_tries << ". When greater than 1,\n"
           << "that many passes will be made and the best time reported.\n";

      return 2;
    }

  try
    {
      if (argc >= 4)
        do_parse (argv[3], loop_size);

      if (argc >= 5)
        do_parse (argv[4], num_tries);

      run_benchmark (argv[1], argv[2], loop_size, num_tries);
    }
  catch (error_exit&)
    {
      return 1;
    }

  return 0;
}


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
