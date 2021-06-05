#include "sequence_clusterer.hpp"
#include "common/minimap_output_reader.hpp"
#include "common/aligner_output_postprocessing.hpp"

namespace sequence_clusterer {

struct gcfg {
    std::string alignment_file;
    std::string output_file;
} cfg;


clipp::group GetCLI() {
  using namespace clipp;

  auto cli = (
      cfg.alignment_file << value("alignment file"),
      cfg.output_file << value("output file")
  );

  return cli;
}

#define WISHED_COLUMNS MColumns, MColumns::match, MColumns::strand,\
          MColumns::Q_name, MColumns::Q_size, MColumns::Q_start, MColumns::Q_end,\
          MColumns::T_name, MColumns::T_size, MColumns::T_start, MColumns::T_end

int main(int argc, char * argv[]) {
    START_BANNER("SPAdes standalone sequence clusterer");

    ReadMinimapOutput();

    INFO("SPAdes standalone sequence clusterer finished");
    return 0;
}

} // namespace sequence_clusterer
