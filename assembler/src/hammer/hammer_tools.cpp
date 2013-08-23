//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * hammer_tools.cpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#include "standard.hpp"
#include <iostream>
#include <fstream>
#include <boost/format.hpp>

#include <time.h>
#include <sys/resource.h>
#include <iomanip>

#include "io/ireadstream.hpp"
#include "mathfunctions.hpp"
#include "valid_kmer_generator.hpp"
#include "globals.hpp"
#include "config_struct_hammer.hpp"
#include "hammer_tools.hpp"
#include "kmer_data.hpp"
#include "read_corrector.hpp"

#include "io/mmapped_writer.hpp"

#include <sys/types.h>
#include <sys/wait.h>

using namespace std;
using namespace hammer;

void HammerTools::InitializeSubKMerPositions() {
  ostringstream log_sstream;
  log_sstream.str("");
  Globals::subKMerPositions = new std::vector<uint32_t>(cfg::get().general_tau + 2);
  for (uint32_t i=0; i < (uint32_t)(cfg::get().general_tau + 1); ++i) {
    Globals::subKMerPositions->at(i) = (i * K / (cfg::get().general_tau + 1) );
    log_sstream << Globals::subKMerPositions->at(i) << " ";
  }
  Globals::subKMerPositions->at(cfg::get().general_tau + 1) = K;
  INFO("Hamming graph threshold tau=" << cfg::get().general_tau << ", k=" << K << ", subkmer positions = [ " << log_sstream.str() << "]" );
}

string HammerTools::getFilename(const string & dirprefix, const string & suffix) {
  ostringstream tmp;
  tmp.str(""); tmp << dirprefix.data() << "/" << suffix.data();
  return tmp.str();
}

string HammerTools::getFilename(const string & dirprefix, unsigned iter_count, const string & suffix ) {
  ostringstream tmp;
  tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data();
  return tmp.str();
}

string HammerTools::getReadsFilename(const std::string & dirprefix, const std::string &fname, unsigned iter_no, const std::string & suffix) {
  ostringstream tmp;
  tmp.str("");

  tmp << dirprefix.data() << "/" << path::basename(fname) << '.' << std::setfill('0') << std::setw(2) << iter_no << "." << suffix.data();
  return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, const string & suffix, int suffix_num ) {
  ostringstream tmp;
  tmp.str(""); tmp << dirprefix.data() << "/" << suffix.data() << "." << suffix_num;
  return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num ) {
  ostringstream tmp;
  tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data() << "." << suffix_num;
  return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num, const string & suffix2 ) {
  ostringstream tmp;
  tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data() << "." << suffix_num << "." << suffix2.data();
  return tmp.str();
}

void HammerTools::CorrectReadsBatch(std::vector<bool> &res,
                                    std::vector<Read> &reads, size_t buf_size,
                                    size_t &changedReads, size_t &changedNucleotides, size_t &uncorrectedNucleotides, size_t &totalNucleotides,
                                    const KMerData &data) {
  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  bool discard_singletons = cfg::get().bayes_discard_only_singletons;
  bool correct_threshold = cfg::get().correct_use_threshold;
  bool discard_bad = cfg::get().correct_discard_bad && !cfg::get().correct_notrim;

  ReadCorrector corrector(data);
# pragma omp parallel for shared(reads, res, data) num_threads(correct_nthreads)
  for (size_t i = 0; i < buf_size; ++i) {
    if (reads[i].size() >= K) {
      res[i] =
          corrector.CorrectOneRead(reads[i],
                                   correct_threshold, discard_singletons, discard_bad);
    } else
      res[i] = false;
  }

  changedReads += corrector.changed_reads();
  changedNucleotides += corrector.changed_nucleotides();
  uncorrectedNucleotides += corrector.uncorrected_nucleotides();
  totalNucleotides += corrector.total_nucleotides();
}

void HammerTools::CorrectReadFile(const KMerData &data,
                                  size_t &changedReads, size_t &changedNucleotides, size_t &uncorrectedNucleotides, size_t &totalNucleotides,
                                  const std::string &fname,
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
  while (!irs.eof()) {
    size_t buf_size = 0;
    for (; buf_size < read_buffer_size && !irs.eof(); ++buf_size) {
      irs >> reads[buf_size];
      reads[buf_size].trimNsAndBadQuality(trim_quality);
    }
    INFO("Prepared batch " << buffer_no << " of " << buf_size << " reads.");

    HammerTools::CorrectReadsBatch(res, reads, buf_size,
                                   changedReads, changedNucleotides, uncorrectedNucleotides, totalNucleotides,
                                   data);

    INFO("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) {
      reads[i].print(*(res[i] ? outf_good : outf_bad), qvoffset);
    }
    INFO("Written batch " << buffer_no);
    ++buffer_no;
  }
}

void HammerTools::CorrectPairedReadFiles(const KMerData &data,
                                         size_t &changedReads, size_t &changedNucleotides, size_t &uncorrectedNucleotides, size_t &totalNucleotides,
                                         const std::string &fnamel, const std::string &fnamer,
                                         ofstream * ofbadl, ofstream * ofcorl, ofstream * ofbadr, ofstream * ofcorr, ofstream * ofunp) {
  int qvoffset = cfg::get().input_qvoffset;
  int trim_quality = cfg::get().input_trim_quality;

  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  size_t read_buffer_size = correct_nthreads * cfg::get().correct_readbuffer;
  std::vector<Read> l(read_buffer_size);
  std::vector<Read> r(read_buffer_size);
  std::vector<bool> left_res(read_buffer_size, false);
  std::vector<bool> right_res(read_buffer_size, false);

  unsigned buffer_no = 0;

  ireadstream irsl(fnamel, qvoffset), irsr(fnamer, qvoffset);
  VERIFY(irsl.is_open()); VERIFY(irsr.is_open());

  while (!irsl.eof() && !irsr.eof()) {
    size_t buf_size = 0;
    for (; buf_size < read_buffer_size && !irsl.eof() && !irsr.eof(); ++buf_size) {
      irsl >> l[buf_size]; irsr >> r[buf_size];
      l[buf_size].trimNsAndBadQuality(trim_quality);
      r[buf_size].trimNsAndBadQuality(trim_quality);
    }
    INFO("Prepared batch " << buffer_no << " of " << buf_size << " reads.");

    HammerTools::CorrectReadsBatch(left_res, l, buf_size,
                                   changedReads, changedNucleotides, uncorrectedNucleotides, totalNucleotides,
                                   data);
    HammerTools::CorrectReadsBatch(right_res, r, buf_size,
                                   changedReads, changedNucleotides, uncorrectedNucleotides, totalNucleotides,
                                   data);

    INFO("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) {
      if (left_res[i] && right_res[i]) {
        l[i].print(*ofcorl, qvoffset);
        r[i].print(*ofcorr, qvoffset);
      } else {
        l[i].print(*(left_res[i] ? ofunp : ofbadl), qvoffset);
        r[i].print(*(right_res[i] ? ofunp : ofbadr), qvoffset);
      }
    }
    INFO("Written batch " << buffer_no);
    ++buffer_no;
  }
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

size_t HammerTools::CorrectAllReads() {
  // Now for the reconstruction step; we still have the reads in rv, correcting them in place.
  size_t changedReads = 0;
  size_t changedNucleotides = 0;
  size_t uncorrectedNucleotides = 0;
  size_t totalNucleotides = 0;

  int correct_nthreads = std::min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);

  INFO("Starting read correction in " << correct_nthreads << " threads.");

  const io::DataSet<> &dataset = cfg::get().dataset;
  io::DataSet<> outdataset;
  size_t idataset = 0;
  for (auto it = dataset.library_begin(), et = dataset.library_end(); it != et; ++it, ++idataset) {
    const auto& lib = *it;
    auto outlib = lib;
    outlib.clear();

    size_t ilib = 0;
    for (auto I = lib.paired_begin(), E = lib.paired_end(); I != E; ++I, ++ilib) {
      INFO("Correcting pair of reads: " << I->first << " and " << I->second);
      std::string unpaired = getLargestPrefix(I->first, I->second) + "_unpaired_" +
                             boost::lexical_cast<std::string>(idataset) + "_" +
                             boost::lexical_cast<std::string>(ilib) + path::extension(I->first);

      std::string outcorl = HammerTools::getReadsFilename(cfg::get().output_dir, I->first,  Globals::iteration_no, "cor.fastq");
      std::string outcorr = HammerTools::getReadsFilename(cfg::get().output_dir, I->second, Globals::iteration_no, "cor.fastq");
      std::string outcoru = HammerTools::getReadsFilename(cfg::get().output_dir, unpaired,  Globals::iteration_no, "cor.fastq");

      std::ofstream ofcorl(outcorl.c_str());
      std::ofstream ofbadl(HammerTools::getReadsFilename(cfg::get().output_dir, I->first,  Globals::iteration_no, "bad.fastq").c_str());
      std::ofstream ofcorr(outcorr.c_str());
      std::ofstream ofbadr(HammerTools::getReadsFilename(cfg::get().output_dir, I->second, Globals::iteration_no, "bad.fastq").c_str());
      std::ofstream ofunp (outcoru.c_str());

      HammerTools::CorrectPairedReadFiles(*Globals::kmer_data,
                                          changedReads, changedNucleotides, uncorrectedNucleotides, totalNucleotides,
                                          I->first, I->second,
                                          &ofbadl, &ofcorl, &ofbadr, &ofcorr, &ofunp);
      outlib.push_back_paired(outcorl, outcorr);
      outlib.push_back_single(outcoru);
    }

    for (auto I = dataset.single_begin(), E = dataset.single_end(); I != E; ++I) {
      INFO("Correcting single reads: " << *I);
      std::string outcor = HammerTools::getReadsFilename(cfg::get().output_dir, *I,  Globals::iteration_no, "cor.fastq");
      std::ofstream ofgood(outcor.c_str());
      std::ofstream ofbad(HammerTools::getReadsFilename(cfg::get().output_dir, *I,  Globals::iteration_no, "bad.fastq").c_str());

      HammerTools::CorrectReadFile(*Globals::kmer_data,
                                   changedReads, changedNucleotides, uncorrectedNucleotides, totalNucleotides,
                                   *I,
                                   &ofgood, &ofbad);
      outlib.push_back_single(outcor);
    }
    outdataset.push_back(outlib);
  }

  cfg::get_writable().dataset = outdataset;

  INFO("Correction done. Changed " << changedNucleotides << " bases in " << changedReads << " reads.");
  INFO("Failed to correct " << uncorrectedNucleotides << " bases out of " << totalNucleotides << ".");
  return changedReads;
}
