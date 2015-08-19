//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef EDIT_H
#define EDIT_H

#include <string>
#include <vector>
#include <iostream>

using namespace::std;

////////////////////////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////////////////////////
extern char* fastqf;
extern char* file_of_fastqf;
extern bool zip_output;
extern int threads;
extern unsigned int chunks_per_thread;
extern int trimq;

////////////////////////////////////////////////////////////////////////////////
// methods
////////////////////////////////////////////////////////////////////////////////
void combine_output(string fqf, string mid_ext, bool uncorrected_out);
void combine_output_paired(string fqf1, string fqf2, string mid_ext, bool uncorrected_out);
void chunkify_fastq(string fqf, vector<streampos> & starts, vector<unsigned long long> & counts);
void guess_quality_scale(string fqf);
vector<string> parse_fastq(vector<string> & fastqfs, vector<int> & pairedend_codes);
void unzip_fastq(string & fqf);
void zip_fastq(string fqf);
vector<string> split(string s, char c);
vector<string> split(string);
int quick_trim(string strqual, vector<int> & untrusted);
#endif
