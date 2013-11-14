// brute_force_clean.hpp Copyright (c) 10.11.2013 Kuprashevich Maksim
// This file is under MIT license.
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
// Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
// OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef BRUTE_FORCE_CLEAN_HPP
#define BRUTE_FORCE_CLEAN_HPP

#include <iostream>

#include "adapter_index.hpp"
#include "io/ireadstream.hpp"
#include "sequence/sequence_tools.hpp"
#include "io/read.hpp"
#include <config_struct_cclean.hpp>
#include <ssw/ssw_cpp.h> // Striped Smith-Waterman aligner

#include "additional.cpp"

using std::string;
using cclean::AdapterIndex;
using additional::WORK_MODE_TYPE;
using additional::NONE;
using additional::SIMPLE;
using additional::BRUTE_SIMPLE;
using additional::BRUTE_DEEP;

class BruteForceClean
{
  // Class that get read with oper() and clean it, if that possible
  public:
    BruteForceClean(std::ostream& output, std::ostream& bed, const string &db,
                    const AdapterIndex &gen, const WORK_MODE_TYPE &brute)
      : output_stream_(output), bed_stream_(bed), adap_gen_(gen), brute_(brute),
        threshold_(cfg::get().mismatch_threshold), db_name_(db),
        aligned_part_fraction_(cfg::get().aligned_part_fraction) {}

    // ReadProcessor class put each read in this operator
    bool operator()(const Read &read);
    inline int aligned() { return cuted_; }

  private:
    inline double GetScoreWithQuality(const StripedSmithWaterman::Alignment &a,
                                      const Quality &qual);

    static int cuted_;

    const WORK_MODE_TYPE &brute_;
    const AdapterIndex &adap_gen_;
    const double aligned_part_fraction_;
    const std::string &db_name_;

    std::ostream &output_stream_;
    std::ostream &bed_stream_;
    int threshold_;
};

#endif // BRUTE_FORCE_CLEAN_HPP
