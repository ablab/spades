#include "logger/log_writers.hpp"

#include "io/reader.hpp"
#include "io/ofastastream.hpp"
#include "io/read_processor.hpp"

#include "adt/concurrent_dsu.hpp"

#include "segfault_handler.hpp"
#include "memory_limit.hpp"

#include "HSeq.hpp"
#include "kmer_data.hpp"
#include "hamcluster.hpp"
#include "valid_hkmer_generator.hpp"
#include "err_helper_table.hpp"
#include "consensus.hpp"
#include "read_corrector.hpp"
#include "config_struct.hpp"

#include "openmp_wrapper.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <yaml-cpp/yaml.h>

void create_console_logger() {
  using namespace logging;

  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

struct UfCmp {
  bool operator()(const std::vector<unsigned> &lhs, const std::vector<unsigned> &rhs) {
    //return (lhs[0] < rhs[0]);
    return lhs.size() > rhs.size();
  }
};

struct CountCmp {
  const KMerData &kmer_data_;

  CountCmp(const KMerData &kmer_data)
      : kmer_data_(kmer_data) {}

  bool operator()(unsigned lhs, unsigned rhs) {
    return kmer_data_[lhs].count > kmer_data_[rhs].count;
  }
};


hammer::HKMer center(const KMerData &data, const std::vector<unsigned>& kmers) {
  hammer::HKMer res;
  namespace numeric = boost::numeric::ublas;

  for (unsigned i = 0; i < hammer::K; ++i) {
    numeric::matrix<double> scores(4, 64, 0);
    for (size_t j = 0; j < kmers.size(); ++j) {
      const hammer::KMerStat &k = data[kmers[j]];
      // FIXME: switch to MLE when we'll have use per-run quality values
#if 1
      scores(k.kmer[i].nucl, k.kmer[i].len) += k.count * (1 - k.qual);
#else
      for (unsigned n = 0; n < 4; ++n)
        for (unsigned l = 1; l < 64; ++l)
          scores(n, l) += k.count * (n == k.kmer[i].nucl && l == k.kmer[i].len ?
                                     log(1 - k.qual) : log(k.qual) - log(4*63 - 1));
#endif
    }

    res[i] = hammer::iontorrent::consensus(scores).first;
  }

  return res;
}

bool assign(KMerData &kmer_data, const std::vector<unsigned> &cluster) {
  hammer::HKMer c = center(kmer_data, cluster);
  bool nonread = false;

  size_t idx = kmer_data.seq_idx(c);
  if (kmer_data[idx].kmer != c) {
#   pragma omp critical
    {
      idx = kmer_data.push_back(hammer::KMerStat(0, c, 1.0));
    }
    nonread = true;
  }

  for (size_t j = 0; j < cluster.size(); ++j)
    kmer_data[cluster[j]].changeto = idx;

  return nonread;
}

void dump(const KMerData &kmer_data, const std::vector<unsigned> &cluster) {
  std::cerr << "{ \n";
  for (size_t j = 0; j < cluster.size(); ++j)
    std::cerr << kmer_data[cluster[j]].kmer << ": (" << kmer_data[cluster[j]].count << ", " << 1 - kmer_data[cluster[j]].qual << "), \n";
  hammer::HKMer c = center(kmer_data, cluster);
  size_t idx = kmer_data.seq_idx(c);
  if (kmer_data[idx].kmer == c)
    std::cerr << "center: ok " << c << '\n';
  else
    std::cerr << "center: not " << kmer_data[idx].kmer << ":" << c << '\n';
  std::cerr << "}" << std::endl;
}

size_t subcluster(KMerData &kmer_data, std::vector<unsigned> &cluster) {
  size_t nonread = 0;

  // First, sort the kmer indicies wrt count
  std::sort(cluster.begin(), cluster.end(), CountCmp(kmer_data));

  // The number of subclusters for now is really dumb: we assume that the quality should be 1.
  size_t k = 0;
  for (size_t i = 0; i < cluster.size(); ++i)
    k += kmer_data[cluster[i]].qual < sqrt(std::numeric_limits<double>::epsilon());

  if (k <= 1)
    return assign(kmer_data, cluster);

#if 0
  dump(kmer_data, cluster);
#endif

  // Find the closest center
  std::vector<std::vector<unsigned> > idx(k, std::vector<unsigned>());
  for (size_t i = 0; i < k; ++i)
    idx[i].push_back(cluster[i]);
  for (size_t i = k; i < cluster.size(); ++i) {
    unsigned dist = std::numeric_limits<unsigned>::max();
    size_t cidx = k;
    hammer::HKMer kmerx = kmer_data[cluster[i]].kmer;
    for (size_t j = 0; j < k; ++j) {
      hammer::HKMer kmery = kmer_data[cluster[j]].kmer;
      unsigned cdist = hammer::distanceHKMer(kmerx.begin(), kmerx.end(),
                                             kmery.begin(), kmery.end());
      if (cdist < dist) {
        cidx = j;
        dist = cdist;
      }
    }
    VERIFY(cidx < k);
    idx[cidx].push_back(cluster[i]);
  }

  for (auto it = idx.begin(), et = idx.end(); it != et; ++it) {
    const std::vector<unsigned> &subcluster = *it;

    if (assign(kmer_data, subcluster)) {
      nonread += 1;
#if 0
      dump(kmer_data, cluster);
      dump(kmer_data, subcluster);
#endif
    }
  }

  return nonread;
}


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
      return current != hammer_config::HammerStage::ReadCorrection;
  }
}


int main(int argc, char** argv) {
  segfault_handler sh;

  srand(42);
  srandom(42);

  try {
    create_console_logger();

    std::string config_file = "hammer-it.cfg";
    if (argc > 1) config_file = argv[1];
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
        std::ofstream ofs("count.kmdata", std::ios::binary);
        kmer_data.binary_write(ofs);
      }
    } else {
      INFO("Loading K-mer index.");
      std::ifstream ifs("count.kmdata", std::ios::binary);
      VERIFY(ifs.good());
      kmer_data.binary_read(ifs);
    }

    std::vector<std::vector<unsigned> > classes;
    if (stage(cfg::get().start_stage, hammer_config::HammerStage::HammingClustering)) {
      ConcurrentDSU uf(kmer_data.size());
      KMerHamClusterer clusterer(cfg::get().tau);
      INFO("Clustering Hamming graph.");
      clusterer.cluster("kmers.hamcls", kmer_data, uf);
      uf.get_sets(classes);
      size_t num_classes = classes.size();
      INFO("Clustering done. Total clusters: " << num_classes);

      if (cfg::get().debug_mode) {
        INFO("Debug mode on. Writing down clusters.");
        std::ofstream ofs("hamming.cls", std::ios::binary);

        ofs.write((char*)&num_classes, sizeof(num_classes));
        for (size_t i=0; i < classes.size(); ++i) {
          size_t sz = classes[i].size();
          ofs.write((char*)&sz, sizeof(sz));
          ofs.write((char*)&classes[i][0], sz * sizeof(classes[i][0]));
        }
      }
    } else {
      INFO("Loading clusters.");
      std::ifstream ifs("hamming.cls", std::ios::binary);
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
#     pragma omp parallel for shared(nonread, classes, kmer_data)
      INFO("Assigning centers");
      for (size_t i = 0; i < classes.size(); ++i) {
        const auto& cluster = classes[i];
#       pragma omp atomic
        nonread += assign(kmer_data, cluster);
      }
#endif
      INFO("Total " << nonread << " nonread kmers were generated");

      if (cfg::get().debug_mode) {
        INFO("Debug mode on. Saving K-mer index.");
        std::ofstream ofs("cluster.kmdata", std::ios::binary);
        kmer_data.binary_write(ofs);
      }
    } else {
      INFO("Loading K-mer index.");
      std::ifstream ifs("cluster.kmdata", std::ios::binary);
      VERIFY(ifs.good());
      kmer_data.binary_read(ifs);
    }

    INFO("Correcting reads.");
    const io::DataSet &dataset = cfg::get().dataset;
    for (auto it = dataset.reads_begin(), et = dataset.reads_end(); it != et; ++it) {
      INFO("Correcting " << *it);
      io::Reader irs(*it, io::PhredOffset);
      io::ofastastream ors(path::append_path(path::parent_path(*it), path::basename(*it) + ".fasta")); // FIXME: Proper filename

      using namespace hammer::correction;
      SingleReadCorrector read_corrector(kmer_data);
      hammer::ReadProcessor(cfg::get().max_nthreads).Run(irs, read_corrector, ors);
    }

#if 0
    std::sort(classes.begin(), classes.end(),  UfCmp());
    for (size_t i = 0; i < classes.size(); ++i) {
      unsigned modes = 0;
      auto& cluster = classes[i];
      const unsigned kCountThreshold = 50;
      for (size_t j = 0; j < cluster.size(); ++j) {
        if (kmer_data[cluster[j]].count > kCountThreshold) {
          ++modes;
          if (modes >= 2)
            break;
        }
      }

      if (modes < 2) continue; // skip uninteresting clusters

      std::sort(cluster.begin(), cluster.end(), CountCmp(kmer_data));

      std::cerr << i << ": { \n";
      for (size_t j = 0; j < cluster.size(); ++j)
        std::cerr << kmer_data[cluster[j]].kmer << ": (" << kmer_data[cluster[j]].count << ", " << 1 - kmer_data[cluster[j]].qual << "), \n";
      hammer::HKMer c = center(kmer_data, cluster);
      size_t idx = kmer_data.seq_idx(c);
      if (kmer_data[idx].kmer == c)
        std::cerr << "center: ok " << c << '\n';
      else
        std::cerr << "center: not " << kmer_data[idx].kmer << ":" << c << '\n';
      std::cerr << "}" << std::endl;
    }
#endif
  } catch (std::bad_alloc const& e) {
    std::cerr << "Not enough memory to run BayesHammer. " << e.what() << std::endl;
    return EINTR;
  } catch (const YAML::Exception &e) {
    std::cerr << "Error reading config file: " << e.what() << std::endl;
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
