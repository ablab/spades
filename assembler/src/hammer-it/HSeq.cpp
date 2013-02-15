#include "HSeq.hpp"

namespace hammer {
namespace iontorrent {

std::vector<HomopolymerRun> toHomopolymerRuns(const std::string &seq) {
  std::vector<HomopolymerRun> runs;
  if (seq.empty())
    return runs;
  runs.reserve(seq.size());

  char nucl = seq[0];
  uint8_t len = 1;
  for (size_t i = 1; i < seq.size(); ++i) {
    if (seq[i] != nucl) {
      runs.push_back(HomopolymerRun(dignucl(nucl), len));
      len = 1;
      nucl = seq[i];
    } else {
      ++len;
    }
  }
  if (len > 0) {
    runs.push_back(HomopolymerRun(dignucl(nucl), len));
  }
  return runs;
}

};
};
