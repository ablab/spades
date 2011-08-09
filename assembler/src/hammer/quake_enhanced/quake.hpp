#include "stdint.h"
#include <string>
#include <vector>
#include "common/io/reader.hpp"
#include "common/read/ireadstream.hpp"

// ToDo => to settings.hpp
const uint32_t kK = 3;

namespace quake_enhanced {
class Quake {
 public:
  Quake() {
    cur_state_ = kInitial;
  }
  void Count(std::string ifile, std::string ofile,
             std::string hash_file_prefix, uint32_t hash_file_number, 
             uint8_t quality_offset, uint8_t quality_threshold);
  void PrepareHists(std::string hist_file, std::string trusted_hist_file,
                    std::string bad_hist_file, uint32_t top_threshold,
                    double average_min);
 private:
  enum QuakeState {kInitial, kCountDone, kRealHistPrepared, kRealHistPrinted,
                   kTrustedHistPrepared, kLimitsCounted, kTrustedFiltered};
  QuakeState cur_state_;
  std::string kmer_count_file;
  std::vector<uint32_t> trusted_hist;
  std::vector<uint32_t> real_hist;

  // Count
  /**
   * This function reads reads from the stream and splits them into
   * k-mers. Then k-mers are written to several file almost
   * uniformly. It is guaranteed that the same k-mers are written to the
   * same files.
   * @param ifs Steam to read reads from.
   * @param ofiles Files to write the result k-mers. They are written
   * one per line.
   */
  void SplitToFiles(ireadstream ifs, const vector<FILE*> &ofiles,
                    uint8_t error_threshold);
  /**
   * This function reads k-mer and calculates number of occurrences for
   * each of them.
   * @param ifile File with k-mer to process. One per line.
   * @param ofile Output file. For each unique k-mer there will be a
   * line with k-mer itself and number of its occurrences.
   */
  void EvalFile(FILE *ifile, FILE *ofile);

  // PrepareHists
  void AddToHist(double freq);
  void PrepareRealHist();
  void PrintRealHist(std::string hist_file);
  void PrepareTrustedHist(std::string trusted_hist_file, 
                          std::string bad_hist_file, uint32_t top_threshold, 
                          double average_min);
};
}
