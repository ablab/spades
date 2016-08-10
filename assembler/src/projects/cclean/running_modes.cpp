//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "running_modes.hpp"

#include <string>
#include <unordered_map>
#include <algorithm>

#include "adapter_index.hpp"
#include "output.hpp"
#include "io/read_processor.hpp"
#include "pipeline/library.hpp"
#include "utils/logger/log_writers.hpp"
#include "job_wrappers.hpp"
#include "brute_force_clean.hpp"

AbstractCclean *Cleaner::getCleaner(std::ofstream *outf_alig_debug,
                                    std::ofstream *outf_bad_deb,
                                    const std::string &db, WorkModeType mode,
                                    unsigned mlen,
                                    const cclean::AdapterIndex &index,
                                    bool deb_info) {
  AbstractCclean *cleaner;  // Creating cleaner for reads
  if (mode == SINGLE_END || mode == SINGLE_END_Q)
    cleaner = new SimpleClean(*outf_alig_debug, *outf_bad_deb, db,
                              mode, mlen, index, deb_info);
  if (mode == BRUTE_SIMPLE || mode == BRUTE_WITH_Q)
    cleaner = new BruteForceClean(*outf_alig_debug, *outf_bad_deb, db,
                                  mode, mlen, index.GetSeqs(), deb_info);
  return cleaner;
}

void Cleaner::ProcessDataset() {
  // Options proceed
  const std::string db = cfg::get().database;
  const WorkModeType mode = getMode();

  cclean::AdapterIndex index;
  cclean::AdapterIndexBuilder().FillAdapterIndex(db, index);

  const io::DataSet<> &dataset = cfg::get().dataset;
  io::DataSet<> outdataset;
  // Proccessing dataset. Iterating through libraries
  for (auto it = dataset.library_begin(), et = dataset.library_end(); it != et; ++it) {
    const io::SequencingLibrary<> &lib = *it;
    io::SequencingLibrary<> outlib = lib;
    outlib.clear();
    // Iterating through paired reads in current library lib
    for (auto I = lib.paired_begin(), E = lib.paired_end(); I != E; ++I) {
      INFO("Correcting pair reads from " << I->first << " and " << I->second);

      const std::string &file_name_l = I->first;
      const std::string &file_name_r = I->second;
      const std::string outcorl = getReadsFilename(cfg::get().output_working_dir,
                                             file_name_l, "correct_l");
      const std::string outcorr = getReadsFilename(cfg::get().output_working_dir,
                                             file_name_r, "correct_r");
      const std::string unpaired = getPureFilename(file_name_l) + "_" +
                                   getPureFilename(file_name_r);
      const std::string outcoru = getReadsFilename(cfg::get().output_working_dir,
                                             unpaired, "correct_u");
      const std::string outbadl = getReadsFilename(cfg::get().output_working_dir,
                                                   file_name_l, "bad");
      const std::string outbadr = getReadsFilename(cfg::get().output_working_dir,
                                                   file_name_r, "bad");

      std::ofstream ofcorl(outcorl.c_str());
      std::ofstream ofbadl(outbadl.c_str());
      std::ofstream ofcorr(outcorr.c_str());
      std::ofstream ofbadr(outbadr.c_str());
      std::ofstream ofunp (outcoru.c_str());

      CorrectPairedReadFiles(index, file_name_l, file_name_r, &ofbadl, &ofcorl,
                             &ofbadr, &ofcorr, &ofunp, mode);
      outlib.push_back_paired(outcorl, outcorr);
      outlib.push_back_single(outcoru);
    }

    for (auto I = lib.single_begin(), E = lib.single_end(); I != E; ++I) {
      INFO("Correcting single reads from " << *I);

      const std::string reads_file_name = *I;
      const std::string outcor = getReadsFilename(cfg::get().output_working_dir,
                                                  reads_file_name, "correct");
      const std::string outbad = getReadsFilename(cfg::get().output_working_dir,
                                                  reads_file_name, "bad");

      std::ofstream ofgood(outcor.c_str());
      std::ofstream ofbad(outbad.c_str());

      CorrectReadFile(index, reads_file_name, &ofgood, &ofbad, mode);
      outlib.push_back_single(outcor);
    }
    outdataset.push_back(outlib);
  }

  cfg::get_writable().dataset = outdataset;
}

void Cleaner::CorrectReadFile(const cclean::AdapterIndex &index,
                              const std::string &fname, std::ofstream *outf_good,
                              std::ofstream *outf_bad, WorkModeType mode) {
  const unsigned nthreads = cfg::get().nthreads;
  const std::string db = cfg::get().database;
  const unsigned mlen = cfg::get().minimum_lenght;
  const size_t read_buffer_size = nthreads * cfg::get().buffer_size;
  std::vector<Read> reads(read_buffer_size);
  std::vector<bool> res(read_buffer_size, false);

  const bool deb_info = cfg::get().debug_information;
  std::string bad_out_debug = "";
  std::string aligned_out_debug = "";
  if (deb_info) {
    // Else ofstreams will be not used, so there is no sense to create empty files
    // So ofstreams will be created with empty strings
    bad_out_debug = getReadsFilename(cfg::get().output_working_dir,
                                     fname, "debug.bad");
    aligned_out_debug = getReadsFilename(cfg::get().output_working_dir,
                                       fname, "debug.alig");
  }
  std::ofstream ofbad_deb(bad_out_debug.c_str());
  std::ofstream ofalig_deb(aligned_out_debug.c_str());

  unsigned buffer_no = 0;
  unsigned count_bad = 0;
  unsigned count_total = 0;

  ireadstream irs(fname);
  VERIFY(irs.is_open());

  AbstractCclean *cleaner = getCleaner(&ofalig_deb, &ofbad_deb, db, mode, mlen,
                                       index, deb_info);

  while (!irs.eof()) {
    unsigned buf_size = 0;
    for (; buf_size < read_buffer_size && !irs.eof(); ++buf_size) {
      irs >> reads[buf_size];
    }
    if(deb_info) INFO("Prepared batch " << buffer_no << " of "
                      << buf_size << " reads.");
    count_bad += CorrectReadsBatch(cleaner, &res, &reads, buf_size, nthreads);
    count_total += buf_size;
    if (deb_info) INFO("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) { // Here output reads in files
      reads[i].print(*(res[i] ? outf_good : outf_bad), Read::PHRED_OFFSET);
    }
    if(deb_info) INFO("Written batch " << buffer_no);
    ++buffer_no;
  }

  delete cleaner;
  // Process info about results
  const double percent_val = static_cast<double>(count_total) / 100.0;
  std::ostringstream percent_bad;
  percent_bad << std::fixed << std::setprecision(2) <<
                   (static_cast<double>(count_bad) / percent_val);
  INFO("Total proceed " + std::to_string(count_total) + ", " +
       std::to_string(count_bad) + " reads (" + percent_bad.str() +
       " percents of total) is bad.");
}

void Cleaner::CorrectPairedReadFiles(const cclean::AdapterIndex &index,
                                     const std::string &fnamel,
                                     const std::string &fnamer, std::ofstream *ofbadl,
                                     std::ofstream *ofcorl, std::ofstream *ofbadr,
                                     std::ofstream *ofcorr, std::ofstream *ofunp,
                                     WorkModeType mode) {
  const unsigned nthreads = cfg::get().nthreads;
  const std::string db = cfg::get().database;
  const unsigned mlen = cfg::get().minimum_lenght;
  const size_t read_buffer_size = nthreads * cfg::get().buffer_size;

  std::vector<Read> left_reads(read_buffer_size);
  std::vector<Read> right_reads(read_buffer_size);
  std::vector<bool> left_res(read_buffer_size, false);
  std::vector<bool> right_res(read_buffer_size, false);

  ireadstream irsl(fnamel);
  ireadstream irsr(fnamer);
  VERIFY(irsl.is_open());
  VERIFY(irsr.is_open());

  const bool deb_info = cfg::get().debug_information;
  std::string bad_out_deb_l = "";
  std::string aligned_out_deb_l = "";
  std::string bad_out_deb_r = "";
  std::string aligned_out_deb_r = "";
  if (deb_info) {
    // Else ofstreams will be not used, so there is no sense to create empty files
    // So ofstreams will be created with empty strings
    bad_out_deb_l = getReadsFilename(cfg::get().output_working_dir,
                                     fnamel, "debug.bad");
    aligned_out_deb_l = getReadsFilename(cfg::get().output_working_dir,
                                       fnamel, "debug.alig");
    bad_out_deb_r = getReadsFilename(cfg::get().output_working_dir,
                                     fnamer, "debug.bad");
    aligned_out_deb_r = getReadsFilename(cfg::get().output_working_dir,
                                       fnamer, "debug.alig");
  }
  std::ofstream ofbad_deb_l(bad_out_deb_l.c_str());
  std::ofstream ofalig_deb_l(aligned_out_deb_l.c_str());
  std::ofstream ofbad_deb_r(bad_out_deb_r.c_str());
  std::ofstream ofalig_deb_r(aligned_out_deb_r.c_str());

  AbstractCclean *cleaner_l = getCleaner(&ofalig_deb_l, &ofbad_deb_l, db, mode,
                                         mlen, index, deb_info);
  AbstractCclean *cleaner_r = getCleaner(&ofalig_deb_r, &ofbad_deb_r, db, mode,
                                         mlen, index, deb_info);
  unsigned buffer_no = 0;
  unsigned count_bad_l = 0;
  unsigned count_bad_r = 0;
  unsigned count_total = 0;

  while (!irsl.eof() && !irsr.eof()) {
    unsigned buf_size = 0;
    for (; buf_size < read_buffer_size && !irsl.eof() &&
         !irsr.eof(); ++buf_size) {
      irsl >> left_reads[buf_size];
      irsr >> right_reads[buf_size];
    }
    if(deb_info) INFO("Prepared batch " << buffer_no << " of " << buf_size
                       << " reads.");

    count_bad_l += CorrectReadsBatch(cleaner_l, &left_res, &left_reads,
                                     buf_size, nthreads);
    count_bad_r += CorrectReadsBatch(cleaner_r, &right_res, &right_reads,
                                     buf_size, nthreads);
    count_total += buf_size;

    if(deb_info) INFO("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) {
      if (left_res[i] && right_res[i]) {
        left_reads[i].print(*ofcorl, Read::PHRED_OFFSET);
        right_reads[i].print(*ofcorr, Read::PHRED_OFFSET);
      }
      else {
        left_reads[i].print(*(left_res[i] ? ofunp : ofbadl),
                            Read::PHRED_OFFSET);
        right_reads[i].print(*(right_res[i] ? ofunp : ofbadr),
                             Read::PHRED_OFFSET);
      }
    }
    if(deb_info) INFO("Written batch " << buffer_no);
    ++buffer_no;
  }

  delete cleaner_l;
  delete cleaner_r;

  // Process info abouts results
  const double percent_val = static_cast<double>(count_total) / 100.0;
  std::ostringstream percent_bad_l;
  std::ostringstream percent_bad_r;
  percent_bad_l << std::fixed << std::setprecision(2) <<
                   (static_cast<double>(count_bad_l) / percent_val);
  percent_bad_r << std::fixed << std::setprecision(2) <<
                   (static_cast<double>(count_bad_r) / percent_val);
  INFO("Total proceed " + std::to_string(count_total) + ", " +
       std::to_string(count_bad_l) + " left reads (" +
       percent_bad_l.str() + " percents of total) is bad" + ", " +
       std::to_string(count_bad_r) + " right reads (" +
       percent_bad_r.str() + " percents of total) is bad.");
}
