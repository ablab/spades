//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <iostream>
#include "utils.cpp"
#include "test_generator.hpp"

int main()
{
    cclean_test::GenerateDataSet("/home/undead/cclean_data/ecoli_mda_lane1.fastq", "/home/undead/cclean_data/dataset_500.fastq");
//    cclean_test::CompareAndPrintFastqFiles("/home/undead/cclean_data/output_bruteforce.fastq",
//                                       "/home/undead/cclean_data/output_simple.fastq");
//    cclean_test::CompareAndPrintFastqFiles("/home/undead/cclean_data/output_bruteforce.fastq",
//                                       "/home/undead/cclean_data/output_100_trimm.fastq");
  return EXIT_SUCCESS;
}
