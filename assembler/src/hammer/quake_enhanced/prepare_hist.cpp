#include <string>
#include "hammer/quake_enhanced/quake.hpp"
using std::string;
using quake_enhanced::Quake;

void Quake::AddToHist(double freq) {
  uint64_t c = static_cast<uint64_t>(freq + 0.5);
  if (real_hist.size() <= c) {
    real_hist.resize(c + 1);
  }
  ++real_hist[c];
}

void Quake::PrintRealHist(string hist_file) {
 if (hist_file != "") {
    FILE *ofile = fopen(hist_file.c_str(), "w");
    for (uint32_t i = 0; i < real_hist.size(); ++i) {
      fprintf(ofile, "%d %d\n",
              i, real_hist[i]);
    }
    fclose(ofile);
  }
  cur_state_ = kRealHistPrinted;
}

void Quake::PrepareRealHist() {
  FILE *ifile = fopen(kmer_count_file.c_str(), "r");
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
  for (uint32_t i = 1; i < real_hist.size() - 1; ++i) {
    if (real_hist[i + 1] > real_hist[i]) {
      fmin = i;
      break;
    }
  }
  if (fmin == -1) {
    printf("Bad histogram\n");
    return;
  }
  int fmax = -1;
  for (uint32_t i = fmin; i < real_hist.size() - 1; ++i) {
    if (real_hist[i + 1] < real_hist[i] 
        && real_hist[i] > real_hist[fmin] * top_threshold) {
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
  while (real_hist[lborder] > real_hist[fmax] * average_min) {
    --lborder;
  }
  while (real_hist[rborder] > real_hist[fmax] * average_min) {
    ++rborder;
  }

  uint32_t mass = 0;
  uint64_t mass_pos = 0;

  for (int i = lborder; i <= rborder; ++i) {
    mass += real_hist[i];
    mass_pos += real_hist[i] * i;
  }
  float average = mass_pos / static_cast<double>(mass);
  printf("Gauss median is at %f\n", average);
  vector<uint32_t> hist_trusted(real_hist);
  for (uint32_t i = 0; static_cast<int>(average - i) >=0; ++i) {
    int where = static_cast<int>(average - i);
    int from = static_cast<int>(average + i + 0.5);
    if (where == from) {
      continue;
    }
    trusted_hist[where] = min(trusted_hist[where + 1], trusted_hist[from]);
  }
  if (trusted_hist_file != "") {
    FILE *trusted_hist_out = fopen(trusted_hist_file.c_str(), "w");
    for (uint32_t i = 0; i < trusted_hist.size(); ++i) {
      fprintf(trusted_hist_out, "%d %d\n", i, trusted_hist[i]);
    }    
    fclose(trusted_hist_out);
  }
  if (bad_hist_file != "") {
    FILE *bad_hist_out = fopen(bad_hist_file.c_str(), "w");
    for (uint32_t i = 0; i < trusted_hist.size(); ++i) {
      fprintf(bad_hist_out, "%d %d\n", i, real_hist[i] - trusted_hist[i]);
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
  if (cur_state_ < kTrustedHistPrepared) {
    PrepareTrustedHist(trusted_hist_file, bad_hist_file, top_threshold,
                       average_min);
  }
}
