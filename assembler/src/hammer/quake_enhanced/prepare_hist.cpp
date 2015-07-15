//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <string>
#include "quake_enhanced/quake.hpp"
#include <cassert>
using std::string;
using quake_enhanced::Quake;

void Quake::AddToHist(double freq) {
  uint64_t c = static_cast<uint64_t>(freq + 0.5);
  if (real_hist_.size() <= c + 1) {
    real_hist_.resize(c + 2);
  }
  ++real_hist_[c];
}

void Quake::PrintRealHist(string hist_file) {
 if (hist_file != "") {
    FILE *ofile = fopen(hist_file.c_str(), "w");
    for (uint32_t i = 0; i < real_hist_.size(); ++i) {
      fprintf(ofile, "%d %d\n",
              i, real_hist_[i]);
    }
    fclose(ofile);
  }
  cur_state_ = kRealHistPrinted;
}

void Quake::PrepareRealHist() {
  FILE *ifile = fopen(kmer_count_file_.c_str(), "r");
  char kmer[kK + 1];
  char format[20];
  snprintf(format, sizeof(format), "%%%ds%%d%%f%%f", kK + 1);
  double freq = -1;
  int count;
  float q_count;
  while (fscanf(ifile, format, kmer, &count, &q_count, &freq) != EOF) {
    AddToHist(freq);
  }
  fclose(ifile);
  cur_state_ = kRealHistPrepared;
}

void Quake::PrepareTrustedHist(string trusted_hist_file, 
                               string bad_hist_file, uint32_t top_threshold, 
                               double average_min) {
  int fmin = -1;
  for (uint32_t i = 1; i < real_hist_.size() - 1; ++i) {
    if (real_hist_[i + 1] > real_hist_[i]) {
      fmin = i;
      break;
    }
  }
  if (fmin == -1) {
    printf("Bad histogram\n");
    return;
  }
  int fmax = -1;
  for (uint32_t i = fmin; i < real_hist_.size() - 1; ++i) {
    if (real_hist_[i + 1] < real_hist_[i] 
        && real_hist_[i] > real_hist_[fmin] * top_threshold) {
      fmax = i;
      break;
    }
  }
  if (fmax == -1) {
    printf("Bad histogram\n");
    return;
  }
  int lborder = fmax;
  int rborder = fmax;
  while (real_hist_[lborder] > real_hist_[fmax] * average_min) {
    --lborder;
  }
  while (real_hist_[rborder] > real_hist_[fmax] * average_min) {
    ++rborder;
  }

  uint32_t mass = 0;
  uint64_t mass_pos = 0;

  for (int i = lborder; i <= rborder; ++i) {
    mass += real_hist_[i];
    mass_pos += real_hist_[i] * i;
  }
  float average = mass_pos / static_cast<double>(mass);
  printf("Gauss median is at %f\n", average);
  trusted_hist_ = real_hist_;
  for (uint32_t i = 0; static_cast<int>(average - i) >=0; ++i) {
    uint32_t where = static_cast<uint32_t>(average - i);
    uint32_t from = static_cast<uint32_t>(average + i + 0.5);
    if (where == from || where + 1 >= trusted_hist_.size()) {
      continue;
    }
    uint32_t value = 0;
    if (from < trusted_hist_.size()) {
      value = trusted_hist_[from];
    }
    trusted_hist_[where] = min(trusted_hist_[where + 1], value);
  }
  if (trusted_hist_file != "") {
    FILE *trusted_hist_out = fopen(trusted_hist_file.c_str(), "w");
    for (uint32_t i = 0; i < trusted_hist_.size(); ++i) {
      if (trusted_hist_[i] > 0) {
        fprintf(trusted_hist_out, "%d %d\n", i, trusted_hist_[i]);
      }
    }    
    fclose(trusted_hist_out);
  }
  if (bad_hist_file != "") {
    FILE *bad_hist_out = fopen(bad_hist_file.c_str(), "w");
    for (uint32_t i = 0; i < trusted_hist_.size(); ++i) {
      if (real_hist_[i] - trusted_hist_[i] > 0) {
        fprintf(bad_hist_out, "%d %d\n", i, real_hist_[i] - trusted_hist_[i]);
      }
    }
    fclose(bad_hist_out);
  }
  cur_state_ = kTrustedHistPrepared;
}

void Quake::PrepareHists(string hist_file, string trusted_hist_file,
                         string bad_hist_file, uint32_t top_threshold,
                         double average_min) {
  assert(cur_state_ >= kCountDone);
  if (cur_state_ < kRealHistPrepared) {
    PrepareRealHist();    
  }
  if (cur_state_ < kRealHistPrinted) {
    PrintRealHist(hist_file);
  }
  PrepareTrustedHist(trusted_hist_file, bad_hist_file, top_threshold,
                     average_min);
}
