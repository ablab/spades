#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>
#include <unordered_map>
#include "common/read/ireadstream.hpp"
#include "common/read/read.hpp"
#include "common/sequence/seq.hpp"
#include "hammer/valid_kmer_generator.hpp"
#include "hammer/quake_enhanced/quake.hpp"
#define SUPPRESS_UNUSED(X) ((void) (X))
using quake_enhanced::Quake;

using std::string;
using std::vector;
using std::unordered_map;
using io::Reader;
using io::SingleRead;

struct KMerInfo {
  int count;
  double q_count;
  double freq;
};

typedef Seq<kK> KMer;
typedef unordered_map<KMer, KMerInfo, KMer::hash> UnorderedMap;

/**
 * This function reads reads from the stream and splits them into
 * k-mers. Then k-mers are written to several file almost
 * uniformly. It is guaranteed that the same k-mers are written to the
 * same files.
 * @param ifs Steam to read reads from.
 * @param ofiles Files to write the result k-mers. They are written
 * one per line.
 */
void Quake::SplitToFiles(ireadstream ifs, const vector<FILE*> &ofiles,
                  uint8_t error_threshold) {
  uint32_t file_number = ofiles.size();
  while (!ifs.eof()) {
    Read r;
    ifs >> r;
    KMer::hash hash_function;
    for (ValidKMerGenerator<kK> gen(r, error_threshold); gen.HasMore(); gen.Next()) {
      KMer kmer = gen.kmer();
      if (KMer::less2()(!kmer, kmer)) {
        kmer = !kmer;
      }
      FILE *cur_file = ofiles[hash_function(kmer()) % file_number];
      KMer::BinWrite(cur_file, kmer);
      double q_count = gen.correct_probability();
      fwrite(&q_count, sizeof(q_count), 1, cur_file);
    }
  }
}

/**
 * This function reads k-mer and calculates number of occurrences for
 * each of them.
 * @param ifile File with k-mer to process. One per line.
 * @param ofile Output file. For each unique k-mer there will be a
 * line with k-mer itself and number of its occurrences.
 */
void Quake::EvalFile(FILE *ifile, FILE *ofile) {
  UnorderedMap stat_map;
  char buffer[kK + 1];
  buffer[kK] = 0;
  KMer kmer;
  while (KMer::BinRead(ifile, &kmer)) {
    KMerInfo &info = stat_map[kmer];
    double q_count = -1;
    bool readed =
        fread(&q_count, sizeof(q_count),
              1, ifile);
    assert(readed == 1);
    SUPPRESS_UNUSED(readed);
    double freq = 0;
    // ToDo 0.5 threshold ==>> command line option
    if (q_count > 0.5) {
      freq = 1 / q_count;
    }
    info.q_count += q_count;
    info.count += 1;
    info.freq += freq;
  }
  for (UnorderedMap::iterator it = stat_map.begin();
       it != stat_map.end(); ++it) {
    const KMerInfo &info = it->second;
    AddToHist(info.freq);
    fprintf(ofile, "%s %d %f %f\n", it->first.str().c_str(),
            info.count, info.q_count, info.freq);
  }
}

void Quake::Count(string ifile_name, string ofile_name, 
                  string hash_file_prefix, uint32_t hash_file_number, 
                  uint8_t quality_offset, uint8_t quality_threshold) {
  vector<FILE*> ofiles(hash_file_number);
  for (uint32_t i = 0; i < hash_file_number; ++i) {
    char filename[50];
    snprintf(filename, sizeof(filename), "%s%u.part", 
             hash_file_prefix.c_str(), i);
    ofiles[i] = fopen(filename, "wb");
    assert(ofiles[i] != NULL && "Too many files to open");
  }
  SplitToFiles(ireadstream(ifile_name, quality_offset),
               ofiles, quality_threshold);
  for (uint32_t i = 0; i < hash_file_number; ++i) {
    fclose(ofiles[i]);
  }
  FILE *ofile = fopen(ofile_name.c_str(), "w");
  assert(ofile != NULL && "Too many files to open");
  for (uint32_t i = 0; i < hash_file_number; ++i) {
    char ifile_name[50];
    snprintf(ifile_name, sizeof(ifile_name), "%s%u.part", 
             hash_file_prefix.c_str(), i);
    FILE *ifile = fopen(ifile_name, "rb");
    EvalFile(ifile, ofile);
    fclose(ifile);
    remove(ifile_name);
  }
  fclose(ofile);
  cur_state_ = kRealHistPrepared;
}
