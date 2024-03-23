//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "hammer_tools.hpp"

#include "valid_kmer_generator.hpp"
#include "globals.hpp"
#include "kmer_data.hpp"
#include "read_corrector.hpp"

#include "io/reads/ireadstream.hpp"
#include "io/kmers/mmapped_writer.hpp"
#include "utils/filesystem/path_helper.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "config_struct_hammer.hpp"

using namespace std;

namespace hammer {

void InitializeSubKMerPositions(int tau) {
  std::ostringstream log_sstream;
  log_sstream.str("");
  Globals::subKMerPositions = new std::vector<uint32_t>(tau + 2);
  for (uint32_t i=0; i < (uint32_t)(tau + 1); ++i) {
    Globals::subKMerPositions->at(i) = (i * K / (tau + 1) );
    log_sstream << Globals::subKMerPositions->at(i) << " ";
  }
  Globals::subKMerPositions->at(tau + 1) = K;
  INFO("Hamming graph threshold tau=" << cfg::get().general_tau << ", k=" << K << ", subkmer positions = [ " << log_sstream.str() << "]" );
}

filesystem::path getFilename(const filesystem::path & dirprefix,
                                  const string & suffix) {
  return dirprefix / suffix;
}

filesystem::path getFilename(const filesystem::path & dirprefix,
                             unsigned iter_count, const string & suffix ) {
  string iter_count_str = to_string(iter_count);
  iter_count_str = iter_count_str.length()==2 ? iter_count_str : '0' + iter_count_str;
  return dirprefix / (iter_count_str + '.' + suffix);
}

filesystem::path getReadsFilename(const filesystem::path & dirprefix, const filesystem::path &fname,
                        unsigned iter_no, const string & suffix) {
  string iter_no_str = to_string(iter_no);
  iter_no_str = iter_no_str.length()==2 ? iter_no_str : "0" + iter_no_str;
  return dirprefix / fname.stem().concat(iter_no_str + "." + suffix);
}

filesystem::path getFilename( const filesystem::path & dirprefix, const string & suffix,
                              int suffix_num) {
  return dirprefix / (suffix + "." + to_string(suffix_num));
}

filesystem::path getFilename( const filesystem::path & dirprefix, int iter_count,
                              const string & suffix, int suffix_num ) {
  string iter_count_str = to_string(iter_count);
  iter_count_str = iter_count_str.length()==2 ? iter_count_str : "0" + iter_count_str;
  return dirprefix / (iter_count_str + "." + suffix + "." + to_string(suffix_num));
}

filesystem::path getFilename( const filesystem::path & dirprefix, int iter_count,
                    const string & suffix, int suffix_num, const string & suffix2 ) {
  string iter_count_str = to_string(iter_count);
  iter_count_str = iter_count_str.length()==2 ? iter_count_str : "0" + iter_count_str;
  return dirprefix / (iter_count_str + "." + suffix + "." + to_string(suffix_num) + "." + suffix2);
}

CorrectionStats CorrectReadsBatch(std::vector<bool> &res,
                       std::vector<Read> &reads, size_t buf_size,
                       const KMerData &data) {
  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  bool discard_singletons = cfg::get().bayes_discard_only_singletons;
  bool correct_threshold = cfg::get().correct_use_threshold;
  bool discard_bad = cfg::get().correct_discard_bad;

  ReadCorrector corrector(data, cfg::get().correct_stats);
# pragma omp parallel for shared(reads, res, data) num_threads(correct_nthreads)
  for (size_t i = 0; i < buf_size; ++i) {
    if (reads[i].size() >= K) {
      res[i] =
          corrector.CorrectOneRead(reads[i],
                                   correct_threshold, discard_singletons, discard_bad);
    } else
      res[i] = false;
  }

  CorrectionStats stats;

  stats.changedReads += corrector.changed_reads();
  stats.changedNucleotides += corrector.changed_nucleotides();
  stats.uncorrectedNucleotides += corrector.uncorrected_nucleotides();
  stats.totalNucleotides += corrector.total_nucleotides();
  return stats;
}

CorrectionStats CorrectReadFile(const KMerData &data,
                     const std::filesystem::path &fname,
                     std::ofstream *outf_good, std::ofstream *outf_bad) {
  int qvoffset = cfg::get().input_qvoffset;
  int trim_quality = cfg::get().input_trim_quality;

  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  size_t read_buffer_size = correct_nthreads * cfg::get().correct_readbuffer;
  std::vector<Read> reads(read_buffer_size);
  std::vector<bool> res(read_buffer_size, false);

  ireadstream irs(fname, qvoffset);
  VERIFY(irs.is_open());

  unsigned buffer_no = 0;
  CorrectionStats stats;
  while (!irs.eof()) {
    size_t buf_size = 0;
    for (; buf_size < read_buffer_size && !irs.eof(); ++buf_size) {
      irs >> reads[buf_size];
      reads[buf_size].trimNsAndBadQuality(trim_quality);
    }
    INFO("Prepared batch " << buffer_no << " of " << buf_size << " reads.");

    stats += CorrectReadsBatch(res, reads, buf_size,
                               data);

    INFO("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) {
      reads[i].print(*(res[i] ? outf_good : outf_bad), qvoffset);
    }
    INFO("Written batch " << buffer_no);
    ++buffer_no;
  }
  return stats;
}

CorrectionStats CorrectPairedReadFiles(const KMerData &data,
                                       const std::filesystem::path &fnamel,
                                       const std::filesystem::path &fnamer,
                                       const std::filesystem::path &fnamea,
                                       ofstream *ofbadl,
                                       ofstream *ofcorl,
                                       ofstream *ofbadr,
                                       ofstream *ofcorr,
                                       ofstream *ofunp,
                                       ofstream *ofcora,
                                       bool has_aux) {
  int qvoffset = cfg::get().input_qvoffset;
  int trim_quality = cfg::get().input_trim_quality;

  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  size_t read_buffer_size = correct_nthreads * cfg::get().correct_readbuffer;
  std::vector<Read> l(read_buffer_size);
  std::vector<Read> r(read_buffer_size);
  std::vector<Read> a(read_buffer_size);
  std::vector<bool> left_res(read_buffer_size, false);
  std::vector<bool> right_res(read_buffer_size, false);

  unsigned buffer_no = 0;

  ireadstream irsl(fnamel, qvoffset), irsr(fnamer, qvoffset);
  VERIFY(irsl.is_open()); VERIFY(irsr.is_open());
  ireadstream irsa(fnamea, qvoffset);
  if (has_aux) {
      VERIFY(irsa.is_open());
  }
  CorrectionStats stats;

  while (!irsl.eof() && !irsr.eof()) {
    size_t buf_size = 0;
    for (; buf_size < read_buffer_size && !irsl.eof() && !irsr.eof(); ++buf_size) {
      irsl >> l[buf_size]; irsr >> r[buf_size];
      l[buf_size].trimNsAndBadQuality(trim_quality);
      r[buf_size].trimNsAndBadQuality(trim_quality);
      if (has_aux) {
          irsa >> a[buf_size];
      }
    }
    INFO("Prepared batch " << buffer_no << " of " << buf_size << " reads.");

    stats += CorrectReadsBatch(left_res, l, buf_size,
                      data);
    stats += CorrectReadsBatch(right_res, r, buf_size,
                      data);

    INFO("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) {
      if (left_res[i] && right_res[i]) {
        l[i].print(*ofcorl, qvoffset);
        r[i].print(*ofcorr, qvoffset);
        if (has_aux) {
          a[i].print(*ofcora, qvoffset);
        }
      } else {
        l[i].print(*(left_res[i] ? ofunp : ofbadl), qvoffset);
        r[i].print(*(right_res[i] ? ofunp : ofbadr), qvoffset);
      }
    }
    INFO("Written batch " << buffer_no);
    ++buffer_no;
  }
  if (!irsl.eof() || !irsr.eof())
      FATAL_ERROR("Pair of read files " << fnamel << " and " << fnamer << " contain unequal amount of reads");
  return stats;
}

std::string getLargestPrefix(const std::string &str1, const std::string &str2) {
  string substr = "";
  for (size_t i = 0; i != str1.size() && i != str2.size(); ++i) {
    if (str1[i] == str2[i])
      substr += str1[i];
    else
      break;
  }
  return substr;
}

std::filesystem::path CorrectSingleReadSet(size_t ilib, size_t iread, const std::filesystem::path &fn,
                                 CorrectionStats &stats) {
  std::filesystem::path usuffix = std::to_string(ilib) + "_" +
                        std::to_string(iread) + ".cor.fastq";

  std::filesystem::path outcor = getReadsFilename(cfg::get().output_dir, fn, Globals::iteration_no, usuffix);
  std::ofstream ofgood(outcor);
  std::ofstream ofbad(getReadsFilename(cfg::get().output_dir, fn, Globals::iteration_no, "bad.fastq"),
                      std::ios::out | std::ios::ate);
  stats += CorrectReadFile(*Globals::kmer_data, fn, &ofgood, &ofbad);
  return outcor;
}

size_t CorrectAllReads() {
  // Now for the reconstruction step; we still have the reads in rv, correcting them in place.
  int correct_nthreads = std::min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);

  INFO("Starting read correction in " << correct_nthreads << " threads.");

  CorrectionStats stats;

  const io::DataSet<> &dataset = cfg::get().dataset;
  io::DataSet<> outdataset;
  size_t ilib = 0;
  for (const auto& lib : dataset.libraries()) {
    auto outlib = lib;
    outlib.clear();

    size_t iread = 0;
    auto aux_iter = lib.aux_begin();
    for (auto I = lib.paired_begin(), E = lib.paired_end(); I != E; ++I, ++iread) {
      INFO("Correcting pair of reads: " << I->first << " and " << I->second);
      std::filesystem::path usuffix =  std::to_string(ilib) + "_" +
                             std::to_string(iread) + ".cor.fastq";

      std::filesystem::path unpaired = getLargestPrefix(I->first, I->second) + "_unpaired.fastq";

      std::filesystem::path outcorl = getReadsFilename(cfg::get().output_dir, I->first,  Globals::iteration_no, usuffix);
      std::filesystem::path outcorr = getReadsFilename(cfg::get().output_dir, I->second, Globals::iteration_no, usuffix);
      std::filesystem::path outcoru = getReadsFilename(cfg::get().output_dir, unpaired,  Globals::iteration_no, usuffix);

      std::ofstream ofcorl(outcorl);
      std::ofstream ofbadl(getReadsFilename(cfg::get().output_dir, I->first,  Globals::iteration_no, "bad.fastq"),
                           std::ios::out | std::ios::ate);
      std::ofstream ofcorr(outcorr);
      std::ofstream ofbadr(getReadsFilename(cfg::get().output_dir, I->second, Globals::iteration_no, "bad.fastq"),
                           std::ios::out | std::ios::ate);
      std::ofstream ofunp (outcoru);

      std::filesystem::path aux_path = getLargestPrefix(I->first, I->second) + "_mock.fastq";
      if (lib.has_aux()) {
        aux_path = *aux_iter;
        ++aux_iter;
      }
      std::filesystem::path outcora = getReadsFilename(cfg::get().output_dir, aux_path,  Globals::iteration_no, usuffix);
      std::ofstream ofcora(outcora);

      stats += CorrectPairedReadFiles(*Globals::kmer_data,
                             I->first, I->second, aux_path,
                             &ofbadl, &ofcorl, &ofbadr, &ofcorr, &ofunp, &ofcora, lib.has_aux());
      outlib.push_back_paired(outcorl, outcorr);
      outlib.push_back_single(outcoru);
      if (lib.has_aux()) {
          outlib.push_back_aux(outcora);
      }
    }

    for (auto I = lib.merged_begin(), E = lib.merged_end(); I != E; ++I, ++iread) {
      INFO("Correcting merged reads: " << *I);
      outlib.push_back_merged(CorrectSingleReadSet(ilib, iread, *I, stats));
    }

    for (auto I = lib.single_begin(), E = lib.single_end(); I != E; ++I, ++iread) {
      INFO("Correcting single reads: " << *I);
      outlib.push_back_single(CorrectSingleReadSet(ilib, iread, *I, stats));
    }
    outdataset.push_back(outlib);
    ilib += 1;
  }

  cfg::get_writable().dataset = outdataset;

  INFO("Correction done. Changed " << stats.changedNucleotides << " bases in " << stats.changedReads << " reads.");
  INFO("Failed to correct " << stats.uncorrectedNucleotides << " bases out of " << stats.totalNucleotides << ".");
  return stats.changedReads;
}

};
