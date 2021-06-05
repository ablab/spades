#include "sequence_clusterer.hpp"
#include "common/aligner_output_reader.hpp"
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

int main(int argc, char * argv[]) {
    START_BANNER("SPAdes standalone sequence clusterer");

    INFO("SPAdes standalone sequence clusterer finished");
    return 0;
}

} // namespace sequence_clusterer
