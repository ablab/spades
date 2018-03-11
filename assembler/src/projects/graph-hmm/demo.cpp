#undef NDEBUG
#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include "graph.hpp"

extern "C" {
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"

#include "hmmer.h"
}

#include <limits>
#include "hmmfile.hpp"

// to compile: gcc this_prog.c -lz
#include <stdio.h>
#include <zlib.h>

#include "utils.hpp"

#include <kseq/kseq.h>
KSEQ_INIT(gzFile, gzread)

#include "demo.hpp"
#include "hmmpath.hpp"

std::string rev_comp(const std::string &s) {
  std::string result;
  std::unordered_map<char, char> rc = {{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}};
  for (auto it = s.rbegin(), end = s.rend(); it != end; ++it) {
    result += rc[*it];
  }
  return result;
}

std::vector<std::string> read_fasta_edges(const std::string &filename, bool add_rc) {
  std::vector<std::string> edges;

  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(filename.c_str(), "r");
  assert(fp);
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    // printf("name: %s\n", seq->name.s);
    // if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
    // printf("seq: %s\n", seq->seq.s);
    // if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
    std::string read = seq->seq.s;

    // U -> T
    for (char &ch : read) {
      if (ch == 'U') {
        ch = 'T';
      }
    }

    edges.push_back(read);
    if (add_rc) {
      edges.push_back(rev_comp(read));
    }
  }

  kseq_destroy(seq);
  gzclose(fp);

  return edges;
}

PathSet<ReversalGraphCursor<Graph::GraphCursor>> find_best_path_rev(const hmm::Fees &fees,
                                                                    const std::vector<ReversalGraphCursor<Graph::GraphCursor>> &initial) {
  return impl::find_best_path(fees, initial);
}

PathSet<ReversalGraphCursor<DBGraph::GraphCursor>> find_best_path_rev(const hmm::Fees &fees,
                                                                      const std::vector<ReversalGraphCursor<DBGraph::GraphCursor>> &initial) {
  return impl::find_best_path(fees, initial);
}

PathSet<DBGraph::GraphCursor> find_best_path(const hmm::Fees &fees,
                                             const std::vector<DBGraph::GraphCursor> &initial) {
  return impl::find_best_path(fees, initial);
}

PathSet<Graph::GraphCursor> find_best_path(const hmm::Fees &fees,
                                           const std::vector<Graph::GraphCursor> &initial) {
  return impl::find_best_path(fees, initial);
}
