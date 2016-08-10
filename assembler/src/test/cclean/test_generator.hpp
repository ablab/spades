//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

﻿#ifndef TEST_GENERATOR_HPP
#define TEST_GENERATOR_HPP

#include <iostream>
#include <string>
#include <map>
#include <exception>
#include <fstream>

#include "running_modes.hpp"
#include "io/read_processor.hpp"
#include "io/read.hpp"
#include "utils/logger/logger.hpp"
#include "brute_force_clean.hpp"
#include "utils/simple_tools.hpp"
#include "utils/logger/log_writers.hpp"
#include "adapter_index.hpp"

constexpr int READS_IN_TEST = 700000;
constexpr int NTH = 5;

enum TestGeneratorType {
  EVERY_NTH_ADAPTER = 0,
  EVERY_READ = 1
};

namespace cclean_test {


void GenerateDataSet(const std::string &input_file, const std::string &output_file,
                     TestGeneratorType mode = EVERY_READ) {
  // Generates dataset with TEST_GENERATOR_TYPE option
  // EVERY_READ: takes to dataset every read in range 0 to READS_IN_TEST
  if (!FileExists(input_file)) {
      std::cout << "File " << input_file << " doesn't exists." << std::endl;
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
      next_read.print(output_stream, Read::PHRED_OFFSET);  // What is that - offset?
      ++index;
    }

    std::cout << "Reads processed " << index << std::endl;
    output_stream.close();
  }

  if (mode == EVERY_NTH_ADAPTER) {
    // TODO: Develop this mode, when bruteforce lerans to output cuted reads
  }
}

void CompareAndPrintFastqFiles(const std::string &new_data,
                               const std::string &old_data) {

  if (!FileExists(new_data)) {
      std::cout << "File " << new_data << " doesn't exists." << std::endl;
      return;
  }

  if (!FileExists(old_data)) {
      std::cout << "File " << old_data << " doesn't exists." << std::endl;
      return;
  }

  ireadstream new_stream(new_data);
  ireadstream old_stream(old_data);
  Read next_new_read;
  Read next_old_read;

  std::unordered_map<std::string, Read> new_reads;
  std::unordered_map<std::string, Read> old_reads;

  while (!new_stream.eof()) {
    new_stream >> next_new_read;
    std::pair<std::string, Read> next_read(next_new_read.getName(), next_new_read);
    new_reads.insert(next_read);
  }

  while (!old_stream.eof()) {
    old_stream >> next_old_read;
    std::pair<std::string, Read> next_read(next_old_read.getName(), next_old_read);
    old_reads.insert(next_read);
  }

  std::cout << "Get " << old_reads.size() << " reads from " << old_data << " and "
            << new_reads.size() << " reads from " << new_data << std::endl;
  std::cout << "--------------------" << std::endl;

  bool match = true;

  for (auto iter = new_reads.begin(); iter != new_reads.end(); ++iter) {

    if (old_reads.count(iter->first)) {
      std::string new_read_seq = iter->second.getSequenceString();
      std::string old_read_seq = old_reads[iter->first].getSequenceString();

      if (new_read_seq != old_read_seq) {
        std::cout << "Read " << iter->first << " in " << new_data << " and " <<
                     old_data << " are mismatched. " << std::endl;
        std::cout << new_data << " :\n" << new_read_seq << std::endl;
        std::cout << old_data << " :\n" << old_read_seq << std::endl;
        std::cout << "--------------------" << std::endl;
        match = false;
      }

      old_reads.erase(iter->first);
    }
    else {
      std::cout << "Read " << iter->first << " is in " << new_data <<
                   " but no in " << old_data << std::endl;
      std::cout << "--------------------" << std::endl;
      match = false;
    }

  }

  for (auto iter = old_reads.begin(); iter != old_reads.end(); ++iter) {
    std::cout << "Read " << iter->first << " is in " << old_data <<
                 " but no in " << new_data << std::endl;
    std::cout << "--------------------" << std::endl;
    match = false;
  }

  if(match)
    std::cout << "Files matchs!" << std::endl;

}

void TestMultithreadsCorrect(const std::string &aligned_output,
                             const std::string &bed, const std::string &input,
                             const std::string &output, const std::string &db,
                             const WorkModeType &mode, int times) {
  // check algorithm on multithreads correctness and reads repeats
  std::ofstream reads_text_output(aligned_output);
  std::ofstream reads_bed(aligned_output);
  ireadstream in(input);
  std::ofstream reads_output(output);
  cclean::AdapterIndex index;
//  cclean::AdapterIndexBuilder().FillAdapterIndex(db, index);
  for (int i = 0; i < times; ++i) {
//    ExactAndAlign(reads_text_output, reads_bed, &in, reads_output, db, index,
//                  mode);
    ireadstream test_input(output);
    std::set<std::string> hash_table;
    Read read;
    while (!test_input.eof()) {
      test_input >> read;
      if (!hash_table.count(read.getName()))  {
        hash_table.insert(read.getName());
      }
      else  {
        std::cout << "Repeated read " << read.getName() << " test " << i
                     << " failed!\n";
        return;
      }
    }
  }
}

bool AreTextFilesDifferent(const std::string &new_data, const std::string &old_data) {
  if (!FileExists(new_data)) {
      std::cout << "File " << new_data << " doesn't exists." << std::endl;
      return false;
  }

  if (!FileExists(old_data)) {
      std::cout << "File " << old_data << " doesn't exists." << std::endl;
      return false;
  }

  std::ifstream old_stream(old_data);
  if (!old_stream.is_open()) {
    std::cout << "Can't open " << old_data << "!" << std::endl;
    return false;
  }

  std::ifstream new_stream(new_data);
  if (!new_stream.is_open()) {
    std::cout << "Can't open " << new_data << "!" << std::endl;
    return false;
  }

  std::string old_line;
  std::string new_line;
  int i = 0;
  bool match = true;

  while ( getline(old_stream, old_line) && getline(new_stream, new_line) ) {
        if (strcmp(old_line.c_str(), new_line.c_str())) {
          std::cout << "File " << new_data << " didn't match with " << old_data
                    << " on line " << i << std::endl;
          match = false;
        }
        ++i;
      }

  return match;
}

// End of namespace
}
#endif // TEST_GENERATOR_HPP
