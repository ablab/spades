// test_generator.hpp Copyright (c) 14.11.2013 Kuprashevich Maksim
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

#ifndef TEST_GENERATOR_HPP
#define TEST_GENERATOR_HPP

#include <iostream>
#include <string>
#include <map>
#include <exception>
#include <fstream>

#include "io/read_processor.hpp"
#include "io/read.hpp"
#include "logger/logger.hpp"
#include "brute_force_clean.hpp"
#include "simple_tools.hpp"
#include "logger/log_writers.hpp"

constexpr int READS_IN_TEST = 2000000;
constexpr int NTH = 5;

enum TEST_GENERATOR_TYPE {
  EVERY_NTH_ADAPTER = 0,
  EVERY_READ = 1
};

namespace cclean_test {


void GenerateDataSet(const string &input_file, const string &output_file,
                     TEST_GENERATOR_TYPE mode = EVERY_READ) {
  // Generates dataset with TEST_GENERATOR_TYPE option
  // EVERY_READ: takes to dataset every read in range 0 to READS_IN_TEST
  if (!FileExists(input_file)) {
      std::cout << "File " << input_file << " doesn't exists." << std::endl;
      return;
  }

  if (!FileExists(output_file)) {
      std::cout << "File " << output_file << " doesn't exists." << std::endl;
      return;
  }

  ireadstream input_stream(input_file);
  std::ofstream output_stream(output_file, std::ofstream::out);
  Read next_read;

  if (mode == EVERY_READ) {
    int index = 0;
    std::cout << "Processing " << READS_IN_TEST << " reads from " << input_file
              << " to " << output_file << "..." << std::endl;

    while (index < READS_IN_TEST && !input_stream.eof()) {
      input_stream >> next_read;
      next_read.print(output_stream, 0);  // What is that - offset?
      ++index;
    }

    std::cout << "Reads processed " << index << std::endl;
    output_stream.close();
  }

  if (mode == EVERY_NTH_ADAPTER) {
    // TODO: Develop this mode, when bruteforce lerans to output cuted reads
  }
}

bool IsFilesDiff(const string &new_data, const string &old_data)
{
  if (!FileExists(new_data)) {
      std::cout << "File " << new_data << " doesn't exists." << std::endl;
      return false;
  }

  if (!FileExists(old_data)) {
      std::cout << "File " << old_data << " doesn't exists." << std::endl;
      return false;
  }

  ireadstream new_stream(new_data);
  ireadstream old_stream(old_data);
  Read next_new_read;
  Read next_old_read;

  bool match = true;
  std::string seq_new;
  std::string seq_old;

  while (!new_stream.eof() && !old_stream.eof()){
    new_stream >> next_new_read;
    old_stream >> next_old_read;

    // Thats means they didn't match
    seq_new = next_new_read.getSequenceString();
    seq_old = next_old_read.getSequenceString();
    if (strcmp(seq_new.c_str(), seq_old.c_str())) {
      std::cout << "Read " << next_new_read.getName() << " didn't match in "
                << "old file " << old_data << " and new file " << new_data
                << endl;
      match = false;
    }
  }

  if (!(new_stream.eof() && old_stream.eof())) {
    std::cout << "Files size didn't match!" << std::endl;
    match = false;
  }

  return match;
}

// End of namespace
}
#endif // TEST_GENERATOR_HPP
