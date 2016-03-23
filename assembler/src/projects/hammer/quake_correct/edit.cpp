//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <sys/stat.h>
#include <fstream>
#include "omp_wrapper.h"
#include <iostream>
#include <sstream>
#include <cstring>
#include "gzstream.h"
#include <vector>
#include "Read.h"

////////////////////////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////////////////////////
// -r, fastq file of reads
char* fastqf = NULL;
// -f, file of fastq files of reads
char* file_of_fastqf = NULL;

// -z, zip output files
bool zip_output = false;

// -q
int Read::quality_scale = -1;
// -p, number of threads
int threads = 4;

// -t
int trimq = 3;

unsigned int chunks_per_thread = 200;


////////////////////////////////////////////////////////////////////////////////
// split
//
// Split on whitespace
////////////////////////////////////////////////////////////////////////////////
vector<string> split(string s) {
  vector<string> splits;
  int split_num = 0;
  bool last_space = true;

  for(int i = 0; i < s.size(); i++) {
    if(s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r') {
      if(!last_space)
    split_num++;
      last_space = true;
    } else {
      if(split_num == splits.size())
    splits.push_back("");
      splits[split_num] += s[i];
      last_space = false;
    }
  }

  return splits;
}


////////////////////////////////////////////////////////////////////////////////
// split
//
// Split on the character c, trying to match Python's split method
////////////////////////////////////////////////////////////////////////////////
vector<string> split(string s, char c)
{
  vector<string> splits;
  splits.push_back("");
  int split_num = 0;

  for(int i = 0; i < s.size(); i++) {
    if(s[i] == c) {
      split_num++;
      splits.push_back("");
    } else {
      splits[split_num] += s[i];
    }
  }
  
  return splits;
}


////////////////////////////////////////////////////////////////////////////////
// unzip_fastq
//
// Unzip read file and remove ".gz" suffix from 'fqf'.
////////////////////////////////////////////////////////////////////////////////
void unzip_fastq(string & fqf) {
  char mycmd[500];

  // rename
  string fqf_zip(fqf);
  fqf.erase(fqf.size()-3);

  // unzip but leave original file
  strcpy(mycmd, "gunzip -c ");
  strcat(mycmd, fqf_zip.c_str());
  strcat(mycmd, " > ");
  strcat(mycmd, fqf.c_str());
  system(mycmd);
}


////////////////////////////////////////////////////////////////////////////////
// zip_fastq
//
// Zip the original read file as well as the corrected
// read files.
////////////////////////////////////////////////////////////////////////////////
void zip_fastq(string fqf) {
  char mycmd[100];

  // gzip fqf
  //strcpy(mycmd, "gzip ");
  //strcat(mycmd, fqf.c_str());
  //system(mycmd);

  // remove unzipped fqf, leaving only zipped
  remove(fqf.c_str());  

  // determine output file
  /*
  string fqf_str(fqf);
  int suffix_index = fqf_str.rfind(".");
  string prefix = fqf_str.substr(0,suffix_index);
  string suffix = fqf_str.substr(suffix_index, fqf_str.size()-suffix_index);
  string pairf = prefix + string(".cor") + suffix;
  string singlef = prefix + string(".cor.single") + suffix;

  // gzip pair
  strcpy(mycmd, "gzip ");
  strcat(mycmd, pairf.c_str());
  system(mycmd);

  // gzip single
  struct stat st_file_info;
  if(stat(singlef.c_str(), &st_file_info) == 0) {
    strcpy(mycmd, "gzip ");
    strcat(mycmd, singlef.c_str());
    system(mycmd);
  }
  */
}


////////////////////////////////////////////////////////////////////////////////
// combine_logs
//
// Combine log files that may be in out_dir into a single log file named
// using fqf.
////////////////////////////////////////////////////////////////////////////////
void combine_logs(string fqf, string out_dir) {
  struct stat st_file_info;
  string log1 = out_dir+"/0.log";
  if(stat(log1.c_str(), &st_file_info) == 0) {
    // format log output
    string logf = fqf + ".log";
    ofstream corlog_out(logf.c_str());

    // combine
    string line;
    for(int t = 0; t < threads*chunks_per_thread; t++) {    
      string tc_file(out_dir+"/");
      stringstream tc_convert;
      tc_convert << t;
      tc_file += tc_convert.str();
      tc_file += ".log";
      
      if(stat(tc_file.c_str(), &st_file_info) == 0) {
    ifstream tc_out(tc_file.c_str());
    while(getline(tc_out, line)) {
      corlog_out << line << endl;
    }
    tc_out.close();
    remove(tc_file.c_str());
      }
    }
    corlog_out.close();
  }
}


////////////////////////////////////////////////////////////////////////////////
// combine_output_stream
//
// Combine output files in 'out_dir' into a single file defined by the given
// stream, and remove those files along the way.
////////////////////////////////////////////////////////////////////////////////
void combine_output_stream(ostream & combine_out, ostream & err_out, string out_dir) {
  string header, seq, mid, qual;
  struct stat st_file_info;
  for(int t = 0; t < threads*chunks_per_thread; t++) {
    string tc_file(out_dir+"/");
    stringstream tc_convert;
    tc_convert << t;
    tc_file += tc_convert.str();

    // if file exists, add to single output
    if(stat(tc_file.c_str(), &st_file_info) == 0) {
      ifstream tc_out(tc_file.c_str());
      while(getline(tc_out, header)) {
    getline(tc_out, seq);
    getline(tc_out, mid);
    getline(tc_out, qual);

    if(!err_out.good() || header.find("error") == -1)
      combine_out << header << endl << seq << endl << mid << endl << qual << endl;
    else
      err_out << header.substr(0,header.find("error")) << endl << seq << endl << mid << endl << qual << endl;
      }
      tc_out.close();
      remove(tc_file.c_str());
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
// combine_output
//
// Combine output files in 'out_dir' into a single file and remove 'out_dir'
////////////////////////////////////////////////////////////////////////////////
void combine_output(string fqf, string mid_ext, bool uncorrected_out) {
  // format output directory
  string path_suffix = split(fqf,'/').back();
  string out_dir("."+path_suffix);

  // format output file
  int suffix_index = fqf.rfind(".");
  string prefix, suffix;
  if(suffix_index == -1) {
    prefix = fqf+".";
    suffix = "";
  } else {
    prefix = fqf.substr(0,suffix_index+1);
    suffix = fqf.substr(suffix_index, fqf.size()-suffix_index);
  }

  string outf;
  string errf;
  if(zip_output) {
    // zipped
    outf = prefix + mid_ext + suffix + ".gz";
    ogzstream combine_out(outf.c_str());    
    ogzstream err_out;
    if(uncorrected_out) {
      errf = prefix + "err" + suffix + ".gz";
      err_out.open(errf.c_str());
    }
    combine_output_stream(combine_out, err_out, out_dir);
    combine_out.close();
    err_out.close();
  } else {
    // normal
    outf = prefix + mid_ext + suffix;
    ofstream combine_out(outf.c_str());
    ofstream err_out;
    if(uncorrected_out) {
      errf = prefix + "err" + suffix;
      err_out.open(errf.c_str());
    }
    combine_output_stream(combine_out, err_out, out_dir);
    combine_out.close();
    err_out.close();
  }

  // log
  combine_logs(fqf, out_dir);
    
  // remove output directory
  rmdir(out_dir.c_str());
}


////////////////////////////////////////////////////////////////////////////////
// combine_output_paired_stream
////////////////////////////////////////////////////////////////////////////////
void combine_output_paired_stream(string fqf1, string fqf2, ostream & pair_out1, ostream & single_out1, ostream & single_err_out1, ostream & err_out1, ostream & pair_out2, ostream & single_out2, ostream & single_err_out2, ostream & err_out2) {
  // format output directories
  string path_suffix1 = split(fqf1, '/').back();
  string out_dir1("."+path_suffix1);
  string path_suffix2 = split(fqf2, '/').back();
  string out_dir2("."+path_suffix2);

  string header1, seq1, mid1, qual1, header2, seq2, mid2, qual2;
  struct stat st_file_info;
  for(int t = 0; t < threads*chunks_per_thread; t++) {
    // format thread-chunk output files
    string tc_file1(out_dir1+"/");
    stringstream tc_convert1;
    tc_convert1 << t;
    tc_file1 += tc_convert1.str();

    string tc_file2(out_dir2+"/");
    stringstream tc_convert2;
    tc_convert2 << t;
    tc_file2 += tc_convert2.str();

    // if file exists, both must
    if(stat(tc_file1.c_str(), &st_file_info) == 0) {
      ifstream tc_out1(tc_file1.c_str());
      ifstream tc_out2(tc_file2.c_str());
     
      while(getline(tc_out1, header1)) {
    // get read1
    getline(tc_out1, seq1);
    getline(tc_out1, mid1);
    getline(tc_out1, qual1);

    // get read2
    if(!getline(tc_out2, header2)) {
      cerr << "Uneven number of reads in paired end read files " << tc_file1.c_str() << " and " << tc_file2.c_str() << endl;
      exit(EXIT_FAILURE);
    }    
    getline(tc_out2, seq2);
    getline(tc_out2, mid2);
    getline(tc_out2, qual2);

    if(header1.find("error") == -1) {
      if(header2.find("error") == -1) {
        // no errors
        pair_out1 << header1 << endl << seq1 << endl << mid1 << endl << qual1 << endl;
        pair_out2 << header2 << endl << seq2 << endl << mid2 << endl << qual2 << endl;
      } else {
        // error in 2        
        single_out1 << header1 << endl << seq1 << endl << mid1 << endl << qual1 << endl;
        if(single_err_out2.good())
          single_err_out2 << header2.substr(0,header2.find("error")) << endl << seq2 << endl << mid2 << endl << qual2 << endl;
      }
    } else {
      if(header2.find("error") == -1) {
        // error in 1
        if(single_err_out1.good())
          single_err_out1 << header1.substr(0,header1.find("error")) << endl << seq1 << endl << mid1 << endl << qual1 << endl;
        single_out2 << header2 << endl << seq2 << endl << mid2 << endl << qual2 << endl;
      } else {
        // error in 1,2
        if(err_out1.good()) {
          err_out1 << header1.substr(0,header1.find("error")) << endl << seq1 << endl << mid1 << endl << qual1 << endl;
          err_out2 << header2.substr(0,header2.find("error")) << endl << seq2 << endl << mid2 << endl << qual2 << endl;
        }
      }
    }
      }
      tc_out1.close();
      tc_out2.close();
      remove(tc_file1.c_str());
      remove(tc_file2.c_str());
    }
  }

  // logs
  combine_logs(fqf1, out_dir1);
  combine_logs(fqf2, out_dir2);  

  // remove output directory
  rmdir(out_dir1.c_str());
  rmdir(out_dir2.c_str());
}


////////////////////////////////////////////////////////////////////////////////
// combine_output_paired
//
// Combine output files in 'out_dir' into a single file and remove 'out_dir'
////////////////////////////////////////////////////////////////////////////////
void combine_output_paired(string fqf1, string fqf2, string mid_ext, bool uncorrected_out) {
  string prefix, suffix;

  if(zip_output) {
    // format output pair file1
    int suffix_index = fqf1.rfind(".");
    if(suffix_index == -1) {
      prefix = fqf1+".";
      suffix = "";
    } else {
      prefix = fqf1.substr(0,suffix_index+1);
      suffix = fqf1.substr(suffix_index, fqf1.size()-suffix_index);
    }
    string outf = prefix + mid_ext + suffix + ".gz";
    ogzstream pair_out1(outf.c_str());
          
    // and single file1
    outf = prefix + mid_ext + "_single" + suffix + ".gz";
    ogzstream single_out1(outf.c_str());

    // and error file1
    ogzstream single_err_out1;
    ogzstream err_out1;
    if(uncorrected_out) {
      outf = prefix + "err_single" + suffix + ".gz";
      single_err_out1.open(outf.c_str());
      outf = prefix + "err" + suffix + ".gz";
      err_out1.open(outf.c_str());
    }
    
    // format output pair file2
    suffix_index = fqf2.rfind(".");
    if(suffix_index == -1) {
      prefix = fqf2+".";
      suffix = "";
    } else {
      prefix = fqf2.substr(0,suffix_index+1);
      suffix = fqf2.substr(suffix_index, fqf2.size()-suffix_index);
    }
    outf = prefix + mid_ext + suffix + ".gz";
    ogzstream pair_out2(outf.c_str());

    // and single file2
    outf = prefix + mid_ext + "_single" + suffix + ".gz";
    ogzstream single_out2(outf.c_str());    

    // and error file1
    ogzstream single_err_out2;
    ogzstream err_out2;
    if(uncorrected_out) {
      outf = prefix + "err_single" + suffix + ".gz";
      single_err_out2.open(outf.c_str());
      outf = prefix + "err" + suffix + ".gz";
      err_out2.open(outf.c_str());
    }

    combine_output_paired_stream(fqf1, fqf2, pair_out1, single_out1, single_err_out1, err_out1, pair_out2, single_out2, single_err_out2, err_out2);

    pair_out1.close();
    pair_out2.close();
    single_out1.close();
    single_out2.close();
  
  } else {
    // format output pair file1
    int suffix_index = fqf1.rfind(".");
    if(suffix_index == -1) {
      prefix = fqf1+".";
      suffix = "";
    } else {
      prefix = fqf1.substr(0,suffix_index+1);
      suffix = fqf1.substr(suffix_index, fqf1.size()-suffix_index);
    }
    string outf = prefix + mid_ext + suffix;
    ofstream pair_out1(outf.c_str());
    
    // and single file1
    outf = prefix + mid_ext + "_single" + suffix;
    ofstream single_out1(outf.c_str());

    // and error file1
    ofstream single_err_out1;
    ofstream err_out1;
    if(uncorrected_out) {
      outf = prefix + "err_single" + suffix;
      single_err_out1.open(outf.c_str());
      outf = prefix + "err" + suffix;
      err_out1.open(outf.c_str());
    }
    
    // format output pair file2
    suffix_index = fqf2.rfind(".");
    if(suffix_index == -1) {
      prefix = fqf1+".";
      suffix = "";
    } else {
      prefix = fqf2.substr(0,suffix_index+1);
      suffix = fqf2.substr(suffix_index, fqf2.size()-suffix_index);
    }
    outf = prefix + mid_ext + suffix;
    ofstream pair_out2(outf.c_str());

    // and single file2
    outf = prefix + mid_ext + "_single" + suffix;
    ofstream single_out2(outf.c_str());

    // and error file2
    ofstream single_err_out2;
    ofstream err_out2;    
    if(uncorrected_out) {
      outf = prefix + "err_single" + suffix;
      single_err_out2.open(outf.c_str());
      outf = prefix + "err" + suffix;
      err_out2.open(outf.c_str());
    }

    combine_output_paired_stream(fqf1, fqf2, pair_out1, single_out1, single_err_out1, err_out1, pair_out2, single_out2, single_err_out2, err_out2);
    
    pair_out1.close();
    pair_out2.close();
    single_out1.close();
    single_out2.close();
  }
}


////////////////////////////////////////////////////////////////////////////////
// chunkify_fastq
//
// Determine start points and sequence counts for all
// chunks to be processed in parallel.
////////////////////////////////////////////////////////////////////////////////
void chunkify_fastq(string fqf, vector<streampos> & starts, vector<unsigned long long> & counts) {
  // count number of sequences
  unsigned long long N = 0;
  ifstream reads_in(fqf.c_str());
  string toss;
  while(getline(reads_in, toss))
    N++;
  reads_in.close();
  N /= 4ULL;

  if(threads*chunks_per_thread > N) {
    // use 1 thread for everything
    counts.push_back(N);
    starts.push_back(0);   
    omp_set_num_threads(1);

  } else {
    // determine counts per thread
    unsigned long long sum = 0;
    for(int i = 0; i < threads*chunks_per_thread-1; i++) {
      counts.push_back(N / (threads*chunks_per_thread));
      sum += counts.back();
    }
    counts.push_back(N - sum);

    // find start points
    reads_in.open(fqf.c_str());
    starts.push_back(reads_in.tellg());
    unsigned long long s = 0;
    unsigned int t = 0;
    while(getline(reads_in,toss)) {
      // sequence
      getline(reads_in, toss);
      // +
      getline(reads_in, toss);
      // quality
      getline(reads_in, toss);
      
      if(++s == counts[t] && t < counts.size()-1) {
    starts.push_back(reads_in.tellg());
    s = 0;
    t++;
      }
      
      // set up parallelism
      omp_set_num_threads(threads);
    }
  }
}


////////////////////////////////////////////////////////////
// guess_quality_scale
//
// Guess at ascii scale of quality values by examining
// a bunch of reads and looking for quality values < 64,
// in which case we set it to 33.
//
// Assuming the file is unzipped.
////////////////////////////////////////////////////////////
void guess_quality_scale(string fqf) {
  string header, seq, mid, strqual;
  int reads_to_check = 10000;
  int reads_checked = 0;
  ifstream reads_in(fqf.c_str());
  while(getline(reads_in, header)) {
    getline(reads_in, seq);
    getline(reads_in, mid);
    getline(reads_in, strqual);
    
    for(int i = 0; i < strqual.size(); i++) {
      if(strqual[i] < 64) {
    cerr << "Guessing quality values are on ascii 33 scale" << endl;
    Read::quality_scale = 33;
    reads_in.close();
    return;
      }
    }

    if(++reads_checked >= reads_to_check)
      break;
  }
  reads_in.close();
  cerr << "Guessing quality values are on ascii 64 scale" << endl;
  Read::quality_scale = 64;
}


////////////////////////////////////////////////////////////////////////////////
// parse_fastq
//
// Accept a single fastq file from input, or parse a file with names of fastq
// files.  For multiple files, attached a paired end code to tell the correction
// method how to handle the file.
////////////////////////////////////////////////////////////////////////////////
vector<string> parse_fastq(vector<string> & fastqfs, vector<int> & pairedend_codes) {
  if(file_of_fastqf != NULL) {
    ifstream ff(file_of_fastqf);
    vector<string> next_fastqf;
    string line;

    while(getline(ff, line) && line.size() > 0) {
      next_fastqf = split(line);

      if(next_fastqf.size() == 1) {
    fastqfs.push_back(next_fastqf[0]);
    pairedend_codes.push_back(0);

      } else if(next_fastqf.size() == 2) {
    fastqfs.push_back(next_fastqf[0]);
    fastqfs.push_back(next_fastqf[1]);
    pairedend_codes.push_back(1);
    pairedend_codes.push_back(2);

      } else {
    cerr << "File of fastq file names must have a single fastq file per line for single reads or two fastqf files per line separated by a space for paired end reads " << endl;
    exit(EXIT_FAILURE);
      }
    }

  } else {
    fastqfs.push_back(string(fastqf));
    pairedend_codes.push_back(0);
  }

  return fastqfs;
}

////////////////////////////////////////////////////////////////////////////////
// quick_trim
//
// Trim the end of the read the way BWA does it.
// Removes affected untrusted k-mers.
// Returns the trimmed length.
////////////////////////////////////////////////////////////////////////////////
int quick_trim(string strqual, vector<int> & untrusted) {
  // find trim index
  int phredq;
  int current_trimfunc = 0;
  int max_trimfunc = 0;
  int trim_length = strqual.size();
  for(int i = strqual.size()-1; i >= 0; i--) {
    //phredq = floor(.5-10*log(1.0 - prob[i])/log(10));
    phredq = strqual[i] - Read::quality_scale;
    current_trimfunc += (trimq - phredq);
    if(current_trimfunc > max_trimfunc) {
      max_trimfunc = current_trimfunc;
      trim_length = i;
    }
  }

  // update untrusted
  for(int i = untrusted.size()-1; i >= 0; i--) {
    if(untrusted[i] > trim_length - bithash::k)
      untrusted.pop_back();
  }

  return trim_length;
}
