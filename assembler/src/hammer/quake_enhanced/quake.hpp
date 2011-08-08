#include "stdint.h"
#include <string>
#include "common/io/reader.hpp"
#include "common/read/ireadstream.hpp"

namespace quake_enhanced {
class Quake {
 public:
  Quake() {
    cur_state_ = kInitial;
  }
  void Count(std::string ifile, std::string ofile,
             std::string hash_file_prefix, uint32_t hash_file_number, 
             uint8_t quality_offset, uint8_t quality_threshold);
 private:
  enum QuakeState {kInitial, kCountDone};
  QuakeState cur_state_;
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
};
}
