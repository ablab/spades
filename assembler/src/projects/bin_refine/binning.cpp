//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning.hpp"

#include "blaze/math/dense/DynamicVector.h"
#include "io/utils/id_mapper.hpp"

#include "io/reads/io_helper.hpp"
#include "utils/filesystem/path_helper.hpp"

#include "csv/csv.h"
#include "utils/stl_utils.hpp"

using namespace debruijn_graph;
using namespace bin_stats;

const std::string Binning::UNBINNED_ID = "0";

EdgeLabels::EdgeLabels(const EdgeId e, const Binning& bin_stats)
        : e(e) {
    auto bins = bin_stats.edges_binning().find(e);
    is_binned = bins != bin_stats.edges_binning().end();
    labels_probabilities.resize(bin_stats.bins().size());

    is_repetitive = false;
    if (is_binned) {
        size_t sz = bins->second.size();
        //size_t sz = bin_stats.multiplicities().at(e);
        for (bin_stats::Binning::BinId bin : bins->second)
            labels_probabilities.set(bin, 1.0 / static_cast<double>(sz));
        is_repetitive = bin_stats.multiplicities().at(e) > 1;
    }
}

void Binning::ScaffoldsToEdges() {
  edges_binning_.clear();
  unbinned_edges_.clear();

  for (const auto &scaffold_entry : scaffolds_binning_) {
      ScaffoldId scaffold = scaffold_entry.first;
      BinId bin_id = scaffold_entry.second;

      for (EdgeId e : scaffolds_paths_.at(scaffold)) {
          if (bin_id != UNBINNED) {
              edges_binning_[e].insert(bin_id);
              edges_binning_[graph_.conjugate(e)].insert(bin_id);
          }
          edges_multiplicity_[e] += 1;
          edges_multiplicity_[graph_.conjugate(e)] += 1;
      }
  }

  // Extract unbinned edges
  for (EdgeId edge : graph_.edges()) {
    if (edges_binning_.count(edge))
      continue;

    unbinned_edges_.insert(edge);
    unbinned_edges_.insert(graph_.conjugate(edge));
  }


  // Determine bin coverage
  size_t nbins = bins_.size();
  blaze::DynamicVector<double> bin_covs(nbins, 0), bin_lens(nbins, 0), m2_covs(nbins, 0);
  for (EdgeId e : graph_.canonical_edges()) {
      auto bin = edges_binning_.find(e);
      if (bin == edges_binning_.end())
          continue;

      // Skip repeats
      if (edges_multiplicity_.at(e) > 1)
          continue;

      VERIFY(bin->second.size() == 1);
      BinId binid = *bin->second.begin();
      bin_covs[binid] += double(graph_.kmer_multiplicity(e));
      bin_lens[binid] += double(graph_.length(e));
      m2_covs[binid] += double(graph_.kmer_multiplicity(e)) * double(graph_.coverage(e));
  }

  for (size_t i = 0; i < nbins; ++i) {
      bin_covs[i] /= bin_lens[i];
      bin_stats_[i].mean_cov = bin_covs[i];
      bin_stats_[i].len = bin_lens[i];
  }

  for (size_t i = 0; i < nbins; ++i) {
      m2_covs[i] /= bin_lens[i];
      bin_stats_[i].m2_cov = m2_covs[i];
      bin_stats_[i].sd_cov = sqrt(m2_covs[i] - bin_covs[i] * bin_covs[i]);
  }
}

void Binning::LoadBinning(const std::string &binning_file,
                          bool cami) {
  fs::CheckFileExistenceFATAL(binning_file);

  scaffolds_binning_.clear();
  bins_.clear();
  std::ifstream binning_reader(binning_file);

  io::CSVReader<2,
                io::trim_chars<' '>, io::no_quote_escape<'\t'>,
                io::throw_on_overflow, io::single_and_empty_line_comment<'#'>>
          reader(binning_file);

  // Read the input in CAMI bioboxes format
  // Parse the header
  if (cami) {
      while (char *cline = reader.next_line()) {
          std::string line(cline);
          if (utils::starts_with(line, "#"))
              continue;

          if (utils::starts_with(utils::str_tolower(line), "@sampleid:")) // SampleID
              sample_id_ = line.substr(strlen("@SAMPLEID:"));
          else if (utils::starts_with(line, "@@")) { // Actual header
              reader.set_header("@@SEQUENCEID", "BINID");
              reader.parse_header_line(io::ignore_extra_column, cline);
              break;
          } else if (!utils::starts_with(line, "@")) {
              FATAL_ERROR("Invalid CAMI bioboxies input format!");
          }
      }
      if (sample_id_.empty())
          FATAL_ERROR("Invalid sample ID!");

      INFO("Sample ID: " << sample_id_);
  } else {
      reader.set_header("@@SEQUENCEID", "BINID");
  }

  std::string scaffold_name;
  BinLabel bin_label;
  BinId max_bin_id = 0;
  while (reader.read_row(scaffold_name, bin_label)) {
      BinId cbin_id;
      if (bin_label == UNBINNED_ID) // unbinned scaffold
          cbin_id = UNBINNED;
      else {
          auto entry = bins_.find(bin_label);
          if (entry == bins_.end()) { // new bin label
              cbin_id = max_bin_id++;
              bin_labels_.emplace(cbin_id, bin_label);
              bins_.emplace(bin_label, cbin_id);
          } else {
              cbin_id = entry->second;
          }
      }

      auto scaffold_entry = scaffolds_.find(scaffold_name);
      if (scaffold_entry == scaffolds_.end()) {
          INFO("Unknown scaffold: " << scaffold_name);
          continue;
      }

      scaffolds_binning_[scaffold_entry->second] = cbin_id;
  }

  ScaffoldsToEdges();
}

void Binning::WriteToBinningFile(const std::string& prefix, uint64_t output_options,
                                 const SoftBinsAssignment &soft_edge_labels, const BinningAssignmentStrategy& assignment_strategy,
                                 const io::IdMapper<std::string> &edge_mapper) {
    std::ofstream out_tsv(fs::append_path(prefix, "binning.tsv"));
    std::ofstream out_bins(fs::append_path(prefix, "bin_stats.tsv"));
    std::ofstream out_weights(fs::append_path(prefix, "bin_weights.tsv"));
    std::ofstream out_edges(fs::append_path(prefix, "edge_weights.tsv"));

    auto weight_sorter = [] (const auto &lhs, const auto &rhs) {
        if (math::eq(rhs.second, lhs.second))
            return lhs.first < rhs.first;

        return rhs.second < lhs.second;
    };

    std::vector<std::pair<ScaffoldId, size_t>> scaffolds_by_length;
    for (const auto &path_entry : scaffolds_paths_) {
        size_t length = graph_.k();
        for (EdgeId e : path_entry.second)
            length += graph_.length(e);

        scaffolds_by_length.emplace_back(path_entry.first, length);
    }

    std::sort(scaffolds_by_length.begin(), scaffolds_by_length.end(),
              [](const auto &lhs, const auto &rhs) {
                  if (lhs.second == rhs.second)
                      return lhs.first < lhs.first;

                  return lhs.second > rhs.second;
              });

    // Add CAMI bioboxes format header, if necessary
    if (output_options & OutputOptions::CAMI)
        out_tsv << "@Version:0.9.0\n"
                << "@SAMPLEID:" << sample_id_ << '\n'
                << "@@SEQUENCEID\tBINID\n";

    for (const auto &entry : scaffolds_by_length) {
        const std::string &scaffold_name = scaffolds_labels_.at(entry.first);
        const auto &bins_weights = scaffolds_bin_weights_.at(entry.first);
        std::vector<BinId> new_bin_id =
                assignment_strategy.ChooseMajorBins(bins_weights,
                                                    soft_edge_labels, *this);
        if (new_bin_id.empty()) {
            if (output_options & OutputOptions::EmitZeroBin)
                out_tsv << scaffold_name << '\t' << UNBINNED_ID << '\n';
        } else {
            if (output_options & OutputOptions::TallMulti) {
                for (BinId bin : new_bin_id)
                    out_tsv << scaffold_name << '\t' << bin_labels_.at(bin) << '\n';
            } else {
                out_tsv << scaffold_name;
                for (BinId bin : new_bin_id)
                    out_tsv << '\t' << bin_labels_.at(bin);
                out_tsv << '\n';
            }
        }

        out_weights << scaffold_name << '\t';
        out_weights << "nz: " << bins_weights.nonZeros();
        std::vector<std::pair<BinId, double>> weights;
        for (const auto &entry : bins_weights)
            weights.emplace_back(entry.index(), entry.value());
        std::sort(weights.begin(), weights.end(), weight_sorter);
        for (const auto &entry : weights)
            out_weights << '\t' << bin_labels_.at(entry.first) << ":" << entry.second;

        out_weights << '\n';
    }

    out_edges.precision(3);
    out_edges << "# edge_id\tinternal_edge_id\tlength\tcoverage\tbinned\twas_binned\trepetitive\tedge probs\n";
    for (EdgeId e : graph_.canonical_edges()) {
        const EdgeLabels& edge_labels = soft_edge_labels.at(e);
        out_edges << edge_mapper[graph_.int_id(e)] << '\t' << e << '\t' << graph_.length(e) + graph_.k() << '\t' << graph_.coverage(e) << '\t'
                  << !unbinned_edges_.count(e) << '\t' << edge_labels.is_binned << '\t' << edge_labels.is_repetitive << '\t'
                  << "nz: " << edge_labels.labels_probabilities.nonZeros();
        std::vector<std::pair<BinId, double>> weights;
        for (const auto &entry : edge_labels.labels_probabilities)
            weights.emplace_back(entry.index(), entry.value());
        std::sort(weights.begin(), weights.end(), weight_sorter);
        for (const auto &entry : weights)
            out_edges << '\t' << bin_labels_.at(entry.first) << ":" << entry.second;
        out_edges << '\n';
    }

    out_bins << "bin\tcoverage mean\tcov sd\tbin length\tcoverage m^2\n";
    for (size_t i = 0; i < bins_.size(); ++i) {
        out_bins << bin_labels_.at(i) << '\t' << bin_stats_[i].mean_cov << '\t' << bin_stats_[i].sd_cov << '\t'
                 << bin_stats_[i].len << '\t' << bin_stats_[i].m2_cov << '\n';
    }
}

void Binning::AssignBins(const SoftBinsAssignment& soft_edge_labels,
                         const BinningAssignmentStrategy& assignment_strategy) {
    assignment_strategy.AssignEdgeBins(soft_edge_labels, *this);
    for (const auto &path_entry : scaffolds_paths_) {
        std::vector<EdgeId> scaffold_path(path_entry.second.begin(), path_entry.second.end());
        scaffolds_bin_weights_[path_entry.first] =
                assignment_strategy.AssignScaffoldBins(scaffold_path,
                                                       soft_edge_labels, *this);
    }
}

SoftBinsAssignment LabelInitializer::InitLabels(const bin_stats::Binning &bin_stats) const {
    SoftBinsAssignment state(bin_stats.graph().max_eid());
    for (debruijn_graph::EdgeId e : g_.canonical_edges()) {
        EdgeLabels labels(e, bin_stats);
        state.emplace(e, labels);
        state.emplace(g_.conjugate(e), std::move(labels));
    }

    return state;
}
LabelInitializer::LabelInitializer(const Graph &g) : g_(g) {}

static double WJaccard(const blaze::CompressedMatrix<double, blaze::rowMajor> &m,
                       size_t i, size_t j) {
    auto lhs = blaze::row(m, i), rhs = blaze::row(m, j);
    double sl = blaze::sum(lhs), sr = blaze::sum(rhs), cap = blaze::sum(lhs * rhs);
    return cap / (sl + sr - cap);
}

static double PJaccard(const blaze::CompressedMatrix<double, blaze::rowMajor> &m,
                       size_t i, size_t j) {
    auto lb = m.begin(i), le = m.end(i), rb = m.begin(j), re = m.end(j);
    double res = 0;
    while (lb != le && rb != re) {
        size_t li = lb->index(), ri = rb->index();
        if (li == ri) {
            double x = lb->value(), y = rb->value();
            double sum = 0;

            if (x * y > 0) {
                auto ljb = m.begin(i), lje = m.end(i), rjb = m.begin(j), rje = m.end(j);
                while (ljb != lje || rjb != rje) {
                    if (ljb != lje && rjb != rje &&
                        ljb->index() == rjb->index()) {
                        sum += blaze::max(ljb->value() * y, rjb->value() * x);
                        ++ljb, ++rjb;
                    } else if (rjb == rje ||
                               (ljb != lje && ljb->index() < rjb->index())) {
                        sum += ljb->value() * y;
                        ++ljb;
                    } else {
                        sum += rjb->value() * x;
                        ++rjb;
                    }
                }

                res += x*y / sum;
            }

            ++lb, ++rb;
        } else if (li < ri)
            ++lb;
        else
            ++rb;
    }

    return res;
}

blaze::DynamicMatrix<double> Binning::BinDistance(const SoftBinsAssignment& soft_bins_assignment,
                                                  bool edges) {
    size_t nbins = bins_.size();
    size_t nrows = 0, nz = 0;

    if (edges) {
        for (EdgeId e : graph_.canonical_edges()) {
            const EdgeLabels& edge_labels = soft_bins_assignment.at(e);
            nrows += 1;
            nz += edge_labels.labels_probabilities.nonZeros();
        }
    } else {
        for (const auto &entry: scaffolds_bin_weights_) {
            //if (entry.second.nonZeros() <= 1)
            //    continue;
            nrows += 1;
            nz += entry.second.nonZeros();
        }
    }

    blaze::CompressedMatrix<double, blaze::rowMajor> bin_probs(nrows, nbins);
    bin_probs.reserve(nz);
    size_t rnz = 0;
    if (edges) {
        size_t row = 0;
        for (EdgeId e : graph_.canonical_edges()) {
            const EdgeLabels& edge_labels = soft_bins_assignment.at(e);
            for (const auto &entry : edge_labels.labels_probabilities) {
                //if (entry.value() < 1e-2)
                //    continue;
                bin_probs.append(row, entry.index(), entry.value());
                rnz += 1;
            }
            bin_probs.finalize(row);
            row += 1;
        }
    } else {
        size_t row = 0;
        for (const auto &entry: scaffolds_bin_weights_) {
            //if (entry.second.nonZeros() <= 1)
            //    continue;

            for (const auto &prob : entry.second) {
                //if (prob.value() < 1e-2)
                //    continue;

                bin_probs.append(row, prob.index(), prob.value());
                rnz += 1;
            }
            bin_probs.finalize(row);
            row += 1;
        }
    }

    bin_probs.transpose();

    // Now we're having single bin in a row
    blaze::DynamicMatrix<double> dist(nbins, nbins);
    for (size_t i = 0; i < nbins; ++i) {
        dist(i, i) = 1.0;
        VERBOSE_POWER_T2(i, 0, "Processed " << i << " bins");
        #pragma omp parallel for
        for (size_t j = i + 1; j < nbins; ++j) {
            dist(i, j) = PJaccard(bin_probs, i, j);
            dist(j, i) = dist(i, j);
        }
    }

    return dist;
}

namespace bin_stats {
std::ostream &operator<<(std::ostream &os, const Binning &stats) {
    const auto &graph = stats.graph();

    os << "Total bins: " << stats.bins().size() << std::endl;
    size_t sum_length = 0;
    size_t unbinned_length = 0;
    for (const EdgeId e : graph.edges()) {
      const size_t length = graph.length(e);
      sum_length += length;
      if (stats.unbinned_edges().count(e)) {
        unbinned_length += length;
      }
    }

    os << "Total edges: " << graph.e_size()
       << ", binned: " << stats.edges_binning().size()
       << ", unbinned: " << stats.unbinned_edges().size()
       << " (" << (static_cast<double>(stats.unbinned_edges().size()) / static_cast<double>(graph.e_size()) * 100) << "%)" << std::endl
       << "Sum edge length: " << sum_length << std::endl
       << "Unbinned edges sum length: " << unbinned_length
       << " (" << (static_cast<double>(unbinned_length) / static_cast<double>(sum_length) * 100) << "%)";

    return os;
}

std::ostream &operator<<(std::ostream &os, const EdgeLabels &labels) {
    os << labels.is_binned << '\t';
    os << "nz: " << labels.labels_probabilities.nonZeros();
    for (const auto &entry : labels.labels_probabilities)
        os << '\t' << entry.index() << ":" << entry.value();

    return os;
}

} // namespace bin_stats
