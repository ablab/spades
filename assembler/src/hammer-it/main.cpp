#include "logger/log_writers.hpp"

#include "read/read.hpp"
#include "read/ireadstream.hpp"

#include "adt/concurrent_dsu.hpp"

#include "HSeq.hpp"
#include "kmer_data.hpp"
#include "hamcluster.hpp"
#include "valid_hkmer_generator.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

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

hammer::HKMer center(const KMerData &data, std::vector<unsigned> kmers) {
  hammer::HKMer res;
  double inf = -std::numeric_limits<double>::infinity();

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

    unsigned nucl = 0, len = 1; double max = inf;
    for (unsigned j = 0; j < 4; ++j)
      for (unsigned k = 1; k < 64; ++k)
        if (scores(j, k) > max) {
          nucl = j;
          len = k;
          max = scores(j, k);
        }
    res[i].nucl = nucl;
    res[i].len = len;
  }

  return res;
}


int main(void) {
  srand(42);
  srandom(42);

  create_console_logger();

  KMerData kmer_data;
  KMerDataCounter(1).FillKMerData(kmer_data);

  ConcurrentDSU uf(kmer_data.size());
  KMerHamClusterer clusterer(1);
  INFO("Clustering Hamming graph.");
  clusterer.cluster("kmers.hamcls", kmer_data, uf);
  std::vector<std::vector<unsigned> > classes;
  uf.get_sets(classes);
  size_t num_classes = classes.size();

#if 1
  std::sort(classes.begin(), classes.end(),  UfCmp());
  for (size_t i = 0; i < classes.size(); ++i) {
    std::cerr << i << ": { \n";
    for (size_t j = 0; j < classes[i].size(); ++j)
      std::cerr << kmer_data[classes[i][j]].kmer << ": (" <<   kmer_data[classes[i][j]].count << ", " << 1 - kmer_data[classes[i][j]].qual << "), \n";
    hammer::HKMer c = center(kmer_data, classes[i]);
    size_t idx = kmer_data.seq_idx(c);
    if (kmer_data[idx].kmer == c)
      std::cerr << "center: ok " << c << '\n';
    else
      std::cerr << "center: not " << kmer_data[idx].kmer << ":" << c << '\n';
    std::cerr << "}" << std::endl;
  }
#endif
  INFO("Clustering done. Total clusters: " << num_classes);

  return 0;
}
