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

#include "adt/concurrent_dsu.hpp"

#include "utils/segfault_handler.hpp"
#include "utils/memory_limit.hpp"

#include "HSeq.hpp"
#include "config_struct.hpp"
#include "err_helper_table.hpp"
#include "io_read_corrector.hpp"
#include "kmer_data.hpp"
#include "penalty_estimator.hpp"
#include "read_corrector_new.hpp"
#include "subcluster.hpp"

#include "utils/parallel/openmp_wrapper.h"

#include "hamcluster_1.h"
#include "quality_metrics.h"
#include "version.hpp"

#include <fstream>
#include <iomanip>

#include <bamtools/api/BamReader.h>
#include <bamtools/api/SamHeader.h>

#include "gamma_poisson_model.hpp"
#include "normal_quality_model.hpp"

void create_console_logger() {
  using namespace logging;

  logger* lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

struct UfCmp {
  bool operator()(const std::vector<unsigned long>& lhs,
                  const std::vector<unsigned long>& rhs) {
    return lhs.size() > rhs.size();
  }
};

using namespace n_gamma_poisson_model;

namespace hammer {
namespace correction {

using namespace n_gamma_poisson_model;
using namespace n_normal_model;

class TKMerDataEstimator {
  KMerData& Data;
  const uint NumFiles;
  const hammer_config::hammer_config& Config;
  std::vector<std::vector<size_t> > Classes;
  NormalClusterModel ClusterModel;

  // This is weird workaround for bug in gcc 4.4.7
  static bool stage(hammer_config::HammerStage start,
                    hammer_config::HammerStage current) {
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

  void SaveKMerData(const std::string &filename = "count.kmdata") {
    INFO("Debug mode on. Saving K-mer index.");
    std::ofstream ofs(fs::append_path(cfg::get().working_dir, filename), std::ios::binary);
    Data.binary_write(ofs);
  }

  void SaveClusters() {
    INFO("Debug mode on. Writing down clusters.");
    std::ofstream ofs(fs::append_path(Config.working_dir, "hamming.cls"),
                      std::ios::binary);
    const size_t num_classes = Classes.size();
    ofs.write((char*)&num_classes, sizeof(num_classes));
    for (size_t i = 0; i < Classes.size(); ++i) {
      size_t sz = Classes[i].size();
      ofs.write((char*)&sz, sizeof(sz));
      ofs.write((char*)&Classes[i][0], sz * sizeof(Classes[i][0]));
    }
  }

  void LoadKMerData(std::string filename) {
    INFO("Loading K-mer index.");
    std::ifstream ifs(fs::append_path(Config.working_dir, filename),
                      std::ios::binary);
    VERIFY(ifs.good());
    Data.binary_read(ifs);
    INFO("Total " << Data.size() << " entries were loader");
  }

  void CountKMers() { KMerDataCounter(NumFiles).FillKMerData(Data); }

  void ClusterHammingGraph() {
    INFO("Clustering Hamming graph.");
    {
      const auto num_threads = cfg::get().max_nthreads;
      TOneErrorClustering oneErrorClustering(Data, num_threads);
      oneErrorClustering.FillClasses(Classes);
    }
    const size_t num_classes = Classes.size();
    INFO("Clustering done. Total clusters: " << num_classes);
  }

  void LoadClusters() {
    INFO("Loading clusters.");
    std::ifstream ifs(fs::append_path(Config.working_dir, "hamming.cls"),
                      std::ios::binary);
    VERIFY(ifs.good());

    size_t num_classes = 0;
    ifs.read((char*)&num_classes, sizeof(num_classes));
    Classes.resize(num_classes);

    for (size_t i = 0; i < num_classes; ++i) {
      size_t sz = 0;
      ifs.read((char*)&sz, sizeof(sz));
      Classes[i].resize(sz);
      ifs.read((char*)&Classes[i][0], sz * sizeof(Classes[i][0]));
    }
  }

  void EstimateGenomicCenters() {
    const auto num_threads = cfg::get().max_nthreads;
    QualityTransform trans;
    n_normal_model::ModelEstimator priorEstimator(Data, cfg::get().max_nthreads,
                                                 50, false);

    ClusterModel = priorEstimator.Estimate(Classes);

    INFO("Subclustering.");
    TGenomicHKMersEstimator genomicHKMersEstimator(Data, ClusterModel, cfg::get().center_type);

#pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < Classes.size(); ++i) {
      auto& cluster = Classes[i];
      genomicHKMersEstimator.ProceedCluster(cluster);
    }
  }

  void CalcGenomicEstimationQuality(ClusteringQuality& quality) {
    const auto num_threads = cfg::get().max_nthreads;
    (void)num_threads;
#pragma omp parallel for num_threads(num_threads)
    for (size_t idx = 0; idx < Data.size(); ++idx) {
      if (Data[idx].count > 3) {
        quality.AddKMer(idx);
      }
    }
  }

 public:
  TKMerDataEstimator(KMerData& kmerData,
                     const hammer_config::hammer_config& config,
                     const uint numFiles = 32)
      : Data(kmerData), NumFiles(numFiles), Config(config) {}

  void Estimate() {
    if (stage(Config.start_stage, hammer_config::HammerStage::KMerCounting)) {
      CountKMers();
      if (Config.debug_mode) {
        SaveKMerData("count.kmdata");
      }
    } else {
      LoadKMerData("count.kmdata");
    }

    if (stage(Config.start_stage,
              hammer_config::HammerStage::HammingClustering)) {
      ClusterHammingGraph();
      if (Config.debug_mode) {
        SaveClusters();
      }
    } else {
      LoadClusters();
    }

    std::unique_ptr<TGenomReferenceOracle> oracle;
    std::unique_ptr<ClusteringQuality> clusteringQuality;
    std::string oraclePath = cfg::get().oracle_path;

    if (oraclePath.length()) {
      oracle.reset(new TGenomReferenceOracle(oraclePath));
      clusteringQuality.reset(new ClusteringQuality(*oracle, Data));
      for (size_t i = 0; i < Classes.size(); ++i) {
        clusteringQuality->AddCluster(Classes[i]);
      }
    }

    if (stage(Config.start_stage, hammer_config::HammerStage::SubClustering)) {
      EstimateGenomicCenters();

      if (clusteringQuality) {
        CalcGenomicEstimationQuality(*clusteringQuality);
        clusteringQuality->Info();
      }

      if (Config.debug_mode) {
        SaveKMerData("cluster.kmdata");
      }
    } else {
      LoadKMerData("cluster.kmdata");
    }
  }

  NormalClusterModel GetClusterModel() const { return ClusterModel; }

  void SaveCenters() {
    std::ofstream fasta_ofs("centers.fasta");
    fasta_ofs << std::fixed << std::setprecision(6) << std::setfill('0');
    std::sort(Classes.begin(), Classes.end(), UfCmp());
    for (size_t i = 0; i < Classes.size(); ++i) {
      auto& cluster = Classes[i];
      std::sort(cluster.begin(), cluster.end(), CountCmp(Data));
      hammer::HKMer c = TGenomicHKMersEstimator::Center(Data, cluster);
      size_t idx = Data.seq_idx(c);
      if (Data[idx].kmer == c) {
        fasta_ofs << '>' << std::setw(6) << i << "-cov_" << std::setw(0)
                  << Data[idx].count << "-qual_" << std::setw(14)
                  << 1.0 - Data[idx].qual;

        if (cluster.size() == 1) {
          fasta_ofs << "_singleton";
        }
        fasta_ofs << '\n' << c << '\n';
      }
    }
  }

#if 0
  void SolidKMerExpansion() {
  INFO("Starting solid k-mers expansion in " << Config.max_nthreads << " threads.");
        while (true) {
            Expander expander(Data);
            const io::DataSet<> &dataset = Config.dataset;
            for (auto I = dataset.reads_begin(), E = dataset.reads_end(); I != E; ++I) {
                io::FileReadStream irs(*I, io::PhredOffset);
                hammer::ReadProcessor rp(Config.max_nthreads);
                rp.Run(irs, expander);
                VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
            }
            INFO("" << expander.changed() << " solid k-mers were generated");
            if (expander.changed() == 0)
                break;
        }
  }
#endif
};

};  // namespace correction
};  // namespace hammer

int main(int argc, char** argv) {
  using namespace hammer::correction;
  using TCorrector = ReadCorrector<GammaPoissonLikelihoodCalcer>;
  using SingleReadsCorrector = SingleReadCorrector<TCorrector>;
  using PairedReadsCorrector = PairedReadCorrector<TCorrector>;

  utils::segfault_handler sh;
  srand(42);
  srandom(42);

  try {
    create_console_logger();
    std::string config_file = "hammer-it.cfg";
    if (argc > 1) config_file = argv[1];
    INFO("Starting IonHammer, built from " SPADES_GIT_REFSPEC
         ", git revision " SPADES_GIT_SHA1);
    INFO("Loading config from " << config_file.c_str());
    cfg::create_instance(config_file);
    INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << cfg::get().max_nthreads);

    // hard memory limit
    const size_t GB = 1 << 30;
    utils::limit_memory(cfg::get().hard_memory_limit * GB);

    KMerData kmerData;
    NormalClusterModel clusterModel;

    {
      TKMerDataEstimator estimator(kmerData, cfg::get());
      estimator.Estimate();
      clusterModel = estimator.GetClusterModel();
    }

    GammaPoissonLikelihoodCalcer::Factory calcerFactory(kmerData);

    INFO("Correcting reads.");
    using namespace hammer::correction;
    typename SingleReadsCorrector::NoDebug debug_pred;
    typename SingleReadsCorrector::SelectAll select_pred;
    const auto& dataset = cfg::get().dataset;
    io::DataSet<> outdataset;
    size_t ilib = 0;
    for (auto it = dataset.library_begin(), et = dataset.library_end();
         it != et; ++it, ++ilib) {
      const auto& lib = *it;
      auto outlib = lib;
      outlib.clear();

      size_t iread = 0;
      // First, correct all the paired FASTQ files
      for (auto I = lib.paired_begin(), E = lib.paired_end(); I != E;
           ++I, ++iread) {
        if (fs::extension(I->first) == ".bam" ||
            fs::extension(I->second) == ".bam") {
          continue;
        }

        INFO("Correcting pair of reads: " << I->first << " and " << I->second);

        std::string usuffix =
            std::to_string(ilib) + "_" + std::to_string(iread) + ".cor.fasta";

        std::string outcorl = fs::append_path(
            cfg::get().output_dir, fs::basename(I->first) + usuffix);
        std::string outcorr = fs::append_path(
            cfg::get().output_dir, fs::basename(I->second) + usuffix);

        io::OFastaPairedStream ors(outcorl, outcorr);

        io::SeparatePairedReadStream irs(I->first, I->second, 0);
        PairedReadsCorrector read_corrector(kmerData, calcerFactory, debug_pred,
                                            select_pred);
        hammer::ReadProcessor(cfg::get().max_nthreads)
            .Run(irs, read_corrector, ors);

        outlib.push_back_paired(outcorl, outcorr);
      }

      // Second, correct all the single FASTQ files
      for (auto I = lib.single_begin(), E = lib.single_end(); I != E;
           ++I, ++iread) {
        if (fs::extension(*I) == ".bam") {
          continue;
        }

        INFO("Correcting " << *I);

        std::string usuffix =
            std::to_string(ilib) + "_" + std::to_string(iread) + ".cor.fasta";

        std::string outcor = fs::append_path(cfg::get().output_dir,
                                               fs::basename(*I) + usuffix);
        io::OFastaReadStream ors(outcor);

        io::FileReadStream irs(*I, io::PhredOffset);
        SingleReadsCorrector read_corrector(kmerData, calcerFactory, debug_pred,
                                            select_pred);
        hammer::ReadProcessor(cfg::get().max_nthreads)
            .Run(irs, read_corrector, ors);

        outlib.push_back_single(outcor);
      }

      // Finally, correct all the BAM stuff in a row
      for (auto I = lib.reads_begin(), E = lib.reads_end(); I != E;
           ++I, ++iread) {
        if (fs::extension(*I) != ".bam") {
          continue;
        }

        INFO("Correcting " << *I);

        std::string usuffix =
            std::to_string(ilib) + "_" + std::to_string(iread) + ".cor.fasta";

        std::string outcor = fs::append_path(cfg::get().output_dir,
                                               fs::basename(*I) + usuffix);
        io::OFastaReadStream ors(outcor);

        BamTools::BamReader bam_reader;
        bam_reader.Open(*I);
        auto header = bam_reader.GetHeader();
        bam_reader.Close();

        SingleReadsCorrector read_corrector(kmerData, calcerFactory, &header,
                                            debug_pred, select_pred);
        io::UnmappedBamStream irs(*I);
        hammer::ReadProcessor(cfg::get().max_nthreads)
            .Run(irs, read_corrector, ors);

        outlib.push_back_single(outcor);
      }

      outdataset.push_back(outlib);
    }
    cfg::get_writable().dataset = outdataset;

    std::string fname = fs::append_path(cfg::get().output_dir, "corrected.yaml");
    INFO("Saving corrected dataset description to " << fname);
    cfg::get_writable().dataset.save(fname);
  } catch (std::bad_alloc const& e) {
    std::cerr << "Not enough memory to run IonHammer. " << e.what()
              << std::endl;
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
