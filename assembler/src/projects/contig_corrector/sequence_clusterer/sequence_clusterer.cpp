#include "sequence_clusterer.hpp"
#include "helpers/minimap_output_reader.hpp"
#include "helpers/aligner_output_postprocessing.hpp"

namespace sequence_clusterer {

using namespace helpers;

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

template<class Columns, Columns ... columns>
bool Filter(Record<Columns, columns ...> const & element) {
    constexpr auto EDGE_LENGTH_ERROR_COEFF = 0.01;
    constexpr auto LENGTHS_ERROR_COEFF = 0.05;
    if (element.template Get<Columns::strand>() != '+')
        return false;
    auto contig_start = element.template Get<Columns::Q_start>();
    auto contig_end = element.template Get<Columns::Q_end>();
    auto contig_len = element.template Get<Columns::Q_size>();

    auto edge_start = element.template Get<Columns::T_start>();
    auto edge_end = element.template Get<Columns::T_end>();
    auto edge_len = element.template Get<Columns::T_size>();
    auto contig_delta = contig_end - contig_start;
    auto edge_delta = edge_end - edge_start;

    auto edge_len_difference = std::abs(edge_delta - edge_len);
    auto contig_len_difference = std::abs(contig_delta - contig_len);
    auto lens_difference = std::abs(contig_delta - edge_delta);

    auto const & edge_title = element.template Get<Columns::T_name>();
    return GetIDY(element) > 0.90 &&
            ((double)edge_len_difference <= (double)edge_len * EDGE_LENGTH_ERROR_COEFF || (double)contig_len_difference <= (double)contig_len * EDGE_LENGTH_ERROR_COEFF) &&
            (double) lens_difference < (double) edge_delta * LENGTHS_ERROR_COEFF &&
            GetCov(edge_title) > 2;
}

int main() {
    START_BANNER("SPAdes standalone sequence clusterer");

    // ReadMinimapOutput();

    INFO("SPAdes standalone sequence clusterer finished");
    return 0;
}

} // namespace sequence_clusterer
