//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/log_writers.hpp"

#include "io/reads/file_reader.hpp"
#include "io/sam/bam_reader.hpp"
#include "io/reads/paired_readers.hpp"
#include "io/reads/osequencestream.hpp"
#include "io/reads/read_processor.hpp"

#include "common/adt/concurrent_dsu.hpp"

#include "utils/segfault_handler.hpp"
#include "utils/memory_limit.hpp"

#include "HSeq.hpp"
#include "kmer_data.hpp"
#include "hamcluster.hpp"
#include "subcluster.hpp"
#include "err_helper_table.hpp"
#include "read_corrector.hpp"
#include "expander.hpp"
#include "config_struct.hpp"

#include "utils/openmp_wrapper.h"

#include "version.hpp"

#include <fstream>
#include <iomanip>

#include <bamtools/api/BamReader.h>
#include <bamtools/api/SamHeader.h>

void create_console_logger() {
  using namespace logging;

  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

struct UfCmp {
  bool operator()(const std::vector<unsigned long> &lhs,
                  const std::vector<unsigned long> &rhs) {
    return lhs.size() > rhs.size();
  }
};

// This is weird workaround for bug in gcc 4.4.7
static bool stage(hammer_config::HammerStage start, hammer_config::HammerStage current) {
  switch (start) {
    case hammer_config::HammerStage::KMerCounting:
      return true;
    case hammer_config::HammerStage::HammingClustering:
      return current != hammer_config::HammerStage::KMerCounting;
    case hammer_config::HammerStage::SubClustering:
      return (current != hammer_config::HammerStage::KMerCounting &&
              current != hammer_config::HammerStage::HammingClustering);
    case hammer_config::HammerStage::ReadCorrection:
      return current == hammer_config::HammerStage::ReadCorrection;
  }
  assert(0);
}

int main(int argc, char** argv) {
  segfault_handler sh;

  srand(42);
  srandom(42);

  try {
    create_console_logger();

    std::string config_file = "hammer-it.cfg";
    if (argc > 1) config_file = argv[1];
    INFO("Starting IonHammer, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);
    INFO("Loading config from " << config_file.c_str());
    cfg::create_instance(config_file);

    // hard memory limit
    const size_t GB = 1 << 30;
    limit_memory(cfg::get().hard_memory_limit * GB);

    KMerData kmer_data;
    if (stage(cfg::get().start_stage, hammer_config::HammerStage::KMerCounting)) {
      // FIXME: Actually it's num_files here
      KMerDataCounter(32).FillKMerData(kmer_data);
      if (cfg::get().debug_mode) {
        INFO("Debug mode on. Saving K-mer index.");
        std::ofstream ofs(path::append_path(cfg::get().working_dir, "count.kmdata"), std::ios::binary);
        kmer_data.binary_write(ofs);
      }
    } else {
      INFO("Loading K-mer index.");
      std::ifstream ifs(path::append_path(cfg::get().working_dir, "count.kmdata"), std::ios::binary);
      VERIFY(ifs.good());
      kmer_data.binary_read(ifs);
      INFO("Total " << kmer_data.size() << " entries were loader");
    }

    std::vector<std::vector<size_t> > classes;
    if (stage(cfg::get().start_stage, hammer_config::HammerStage::HammingClustering)) {
      ConcurrentDSU uf(kmer_data.size());
      KMerHamClusterer clusterer(cfg::get().tau);
      INFO("Clustering Hamming graph.");
      clusterer.cluster(path::append_path(cfg::get().working_dir, "kmers.hamcls"), kmer_data, uf);
      uf.get_sets(classes);
      size_t num_classes = classes.size();
      INFO("Clustering done. Total clusters: " << num_classes);

      if (cfg::get().debug_mode) {
        INFO("Debug mode on. Writing down clusters.");
        std::ofstream ofs(path::append_path(cfg::get().working_dir, "hamming.cls"), std::ios::binary);

        ofs.write((char*)&num_classes, sizeof(num_classes));
        for (size_t i=0; i < classes.size(); ++i) {
          size_t sz = classes[i].size();
          ofs.write((char*)&sz, sizeof(sz));
          ofs.write((char*)&classes[i][0], sz * sizeof(classes[i][0]));
        }
      }
    } else {
      INFO("Loading clusters.");
      std::ifstream ifs(path::append_path(cfg::get().working_dir, "hamming.cls"), std::ios::binary);
      VERIFY(ifs.good());

      size_t num_classes = 0;
      ifs.read((char*)&num_classes, sizeof(num_classes));
      classes.resize(num_classes);

      for (size_t i = 0; i < num_classes; ++i) {
        size_t sz = 0;
        ifs.read((char*)&sz, sizeof(sz));
        classes[i].resize(sz);
        ifs.read((char*)&classes[i][0], sz * sizeof(classes[i][0]));
      }
    }

    size_t singletons = 0;
    for (size_t i = 0; i < classes.size(); ++i)
      if (classes[i].size() == 1)
        singletons += 1;
    INFO("Singleton clusters: " << singletons);

    if (stage(cfg::get().start_stage, hammer_config::HammerStage::SubClustering)) {
      size_t nonread = 0;
#if 1
      INFO("Subclustering.");
#     pragma omp parallel for shared(nonread, classes, kmer_data)
      for (size_t i = 0; i < classes.size(); ++i) {
        auto& cluster = classes[i];

#       pragma omp atomic
        nonread += subcluster(kmer_data, cluster);
      }
#else
      INFO("Assigning centers");
#     pragma omp parallel for shared(nonread, classes, kmer_data)
      for (size_t i = 0; i < classes.size(); ++i) {
        const auto& cluster = classes[i];
#       pragma omp atomic
        nonread += assign(kmer_data, cluster);
      }
#endif
      INFO("Total " << nonread << " nonread kmers were generated");

      if (cfg::get().debug_mode) {
        INFO("Debug mode on. Saving K-mer index.");
        std::ofstream ofs(path::append_path(cfg::get().working_dir, "cluster.kmdata"), std::ios::binary);
        kmer_data.binary_write(ofs);
      }
    } else {
      INFO("Loading K-mer index.");
      std::ifstream ifs(path::append_path(cfg::get().working_dir, "cluster.kmdata"), std::ios::binary);
      VERIFY(ifs.good());
      kmer_data.binary_read(ifs);
      INFO("Total " << kmer_data.size() << " entries were loader");
    }

#if 0
    INFO("Starting solid k-mers expansion in " << cfg::get().max_nthreads << " threads.");
    while (true) {
        Expander expander(kmer_data);
        const io::DataSet<> &dataset = cfg::get().dataset;
        for (auto I = dataset.reads_begin(), E = dataset.reads_end(); I != E; ++I) {
            io::FileReadStream irs(*I, io::PhredOffset);
            hammer::ReadProcessor rp(cfg::get().max_nthreads);
            rp.Run(irs, expander);
            VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
        }
        INFO("" << expander.changed() << " solid k-mers were generated");
        if (expander.changed() == 0)
            break;
    }
#endif

#if 0
    std::ofstream fasta_ofs("centers.fasta");
    fasta_ofs << std::fixed << std::setprecision(6) << std::setfill('0');
    std::sort(classes.begin(), classes.end(),  UfCmp());
    for (size_t i = 0; i < classes.size(); ++i) {
      auto& cluster = classes[i];
      std::sort(cluster.begin(), cluster.end(), CountCmp(kmer_data));
      hammer::HKMer c = center(kmer_data, cluster);
      size_t idx = kmer_data.seq_idx(c);
      if (kmer_data[idx].kmer == c) {
        fasta_ofs << '>' << std::setw(6) << i
                  << "-cov_" << std::setw(0) << kmer_data[idx].count
                  << "-qual_" << 1.0 - kmer_data[idx].qual;

        if (cluster.size() == 1)
          fasta_ofs << "_singleton";
        fasta_ofs << '\n' << c << '\n';
      }
    }
#endif

    INFO("Correcting reads.");
    using namespace hammer::correction;
    SingleReadCorrector::NoDebug debug_pred;
    SingleReadCorrector::SelectAll select_pred;
    const auto& dataset = cfg::get().dataset;
    io::DataSet<> outdataset;
    size_t ilib = 0;
    for (auto it = dataset.library_begin(), et = dataset.library_end(); it != et; ++it, ++ilib) {
      const auto& lib = *it;
      auto outlib = lib;
      outlib.clear();

      size_t iread = 0;
      // First, correct all the paired FASTQ files
      for (auto I = lib.paired_begin(), E = lib.paired_end(); I != E; ++I, ++iread) {
          if (path::extension(I->first) == ".bam" || path::extension(I->second) == ".bam")
              continue;

          INFO("Correcting pair of reads: " << I->first << " and " << I->second);

          std::string usuffix = std::to_string(ilib) + "_" +
                                std::to_string(iread) + ".cor.fasta";

          std::string outcorl = path::append_path(cfg::get().output_dir, path::basename(I->first) + usuffix);
          std::string outcorr = path::append_path(cfg::get().output_dir, path::basename(I->second) + usuffix);

          io::PairedOutputSequenceStream ors(outcorl, outcorr);

          io::SeparatePairedReadStream irs(I->first, I->second, 0, false, false);
          PairedReadCorrector read_corrector(kmer_data, debug_pred, select_pred);
          hammer::ReadProcessor(cfg::get().max_nthreads).Run(irs, read_corrector, ors);

          outlib.push_back_paired(outcorl, outcorr);
      }

      // Second, correct all the single FASTQ files
      for (auto I = lib.single_begin(), E = lib.single_end(); I != E; ++I, ++iread) {
          if (path::extension(*I) == ".bam")
              continue;

          INFO("Correcting " << *I);

          std::string usuffix = std::to_string(ilib) + "_" +
                                std::to_string(iread) + ".cor.fasta";

          std::string outcor = path::append_path(cfg::get().output_dir, path::basename(*I) + usuffix);
          io::osequencestream ors(outcor);

          io::FileReadStream irs(*I, io::PhredOffset);
          SingleReadCorrector read_corrector(kmer_data, debug_pred, select_pred);
          hammer::ReadProcessor(cfg::get().max_nthreads).Run(irs, read_corrector, ors);

          outlib.push_back_single(outcor);
      }

      // Finally, correct all the BAM stuff in a row
      for (auto I = lib.reads_begin(), E = lib.reads_end(); I != E; ++I, ++iread) {
        if (path::extension(*I) != ".bam")
              continue;

        INFO("Correcting " << *I);

        std::string usuffix = std::to_string(ilib) + "_" +
                              std::to_string(iread) + ".cor.fasta";

        std::string outcor = path::append_path(cfg::get().output_dir, path::basename(*I) + usuffix);
        io::osequencestream ors(outcor);

        BamTools::BamReader bam_reader;
        bam_reader.Open(*I);
        auto header = bam_reader.GetHeader();
        bam_reader.Close();

        SingleReadCorrector read_corrector(kmer_data, &header, debug_pred, select_pred);
        io::UnmappedBamStream irs(*I);
        hammer::ReadProcessor(cfg::get().max_nthreads).Run(irs, read_corrector, ors);

        outlib.push_back_single(outcor);
      }

      outdataset.push_back(outlib);
    }
    cfg::get_writable().dataset = outdataset;

    std::string fname = path::append_path(cfg::get().output_dir, "corrected.yaml");
    INFO("Saving corrected dataset description to " << fname);
    cfg::get_writable().dataset.save(fname);

#if 0
    std::sort(classes.begin(), classes.end(),  UfCmp());
    for (size_t i = 0; i < classes.size(); ++i) {
      auto& cluster = classes[i];
      std::sort(cluster.begin(), cluster.end(), CountCmp(kmer_data));
      dump(kmer_data, cluster);
    }
#endif
  } catch (std::bad_alloc const& e) {
    std::cerr << "Not enough memory to run IonHammer. " << e.what() << std::endl;
    return EINTR;
  } catch (std::exception const& e) {
    std::cerr << "Exception caught " << e.what() << std::endl;
    return EINTR;
  } catch (...) {
    std::cerr << "Unknown exception caught " << std::endl;
    return EINTR;
  }

  return 0;
}
