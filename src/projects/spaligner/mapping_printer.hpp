//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "alignment/pacbio/g_aligner.hpp"
#include "assembly_graph/core/graph.hpp"
#include "io/utils/edge_namer.hpp"
#include "io/utils/id_mapper.hpp"

#include <filesystem>
#include <fstream>

namespace sensitive_aligner {

class MappingPrinter {
 public:

  MappingPrinter(const debruijn_graph::ConjugateDeBruijnGraph &g,
                 const io::CanonicalEdgeHelper<debruijn_graph::Graph> &edge_namer,
                 const std::filesystem::path &output_dir)
    : g_(g), edge_namer_(edge_namer), output_dir_(output_dir)
  {}

  virtual void SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) = 0;

  virtual ~MappingPrinter () {};

 protected:

  std::string StrId(const EdgeId &e) const;

  const debruijn_graph::ConjugateDeBruijnGraph &g_;
  const io::CanonicalEdgeHelper<debruijn_graph::Graph> &edge_namer_;
  std::filesystem::path output_dir_;
  std::ofstream output_file_;
};

class MappingPrinterTSV: public MappingPrinter {
 public:
  MappingPrinterTSV(const debruijn_graph::ConjugateDeBruijnGraph &g,
                    const io::CanonicalEdgeHelper<debruijn_graph::Graph> &edge_namer,
                    const std::filesystem::path &output_dir)
    : MappingPrinter(g, edge_namer, output_dir) {
    create_directories(output_dir_);
    output_file_.open(output_dir_ / "alignment.tsv", std::ofstream::out);
  }

  void SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) override;

  ~MappingPrinterTSV() {
    output_file_.close();
  }
};


class MappingPrinterFasta: public MappingPrinter {
 public:
  MappingPrinterFasta(const debruijn_graph::ConjugateDeBruijnGraph &g,
                    const io::CanonicalEdgeHelper<debruijn_graph::Graph> &edge_namer,
                    const std::filesystem::path &output_dir)
    : MappingPrinter(g, edge_namer, output_dir) {
    create_directories(output_dir_);
    output_file_.open(output_dir_ / "alignment.fasta", std::ofstream::out);
  }

  void SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) override;

  ~MappingPrinterFasta() {
    output_file_.close();
  }
};

class MappingPrinterGPA : public MappingPrinter {
 public:
  MappingPrinterGPA(const debruijn_graph::ConjugateDeBruijnGraph &g,
                    const io::CanonicalEdgeHelper<debruijn_graph::Graph> &edge_namer,
                    const std::filesystem::path &output_dir)
    : MappingPrinter(g, edge_namer, output_dir) {
    create_directory(output_dir_);
    output_file_.open(output_dir_ / "alignment.gpa", std::ofstream::out);
    output_file_ << "H\n";
  }

  std::string Print(std::map<std::string, std::string> &line) const;


  std::string FormCigar(const std::string &read, const std::string &aligned) const;

  void FormEdgeCigar(const std::string &subread, const std::string &path_seq, const std::vector<size_t> &edgeblocks,
                    std::vector<std::string> &edgecigar, std::vector<Range> &edgeranges) const;

  void GeneratePath(const std::vector<debruijn_graph::EdgeId> &path,
               const PathRange &path_range,
               std::string &aligned, std::vector<size_t> &edgeblocks) const;

  std::string GenerateSubread(const Sequence &read, const PathRange &path_range) const;

  std::string FormGPAOutput(const io::SingleRead &read,
                            const std::vector<debruijn_graph::EdgeId> &path,
                            const std::vector<std::string> &edgecigar,
                            const std::vector<Range> &edgeranges,
                            int &nameIndex, const PathRange &path_range) const;

  void SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) override;

  ~MappingPrinterGPA() {
    output_file_.close();
  }

};

class MappingPrinterHub {
 public:
  MappingPrinterHub(const debruijn_graph::ConjugateDeBruijnGraph &g,
                    const io::CanonicalEdgeHelper<debruijn_graph::Graph> &edge_namer,
                    const std::filesystem::path &output_dir,
                    const std::string formats) {
    if (formats.find("tsv") != std::string::npos) {
      mapping_printers_.push_back(new MappingPrinterTSV(g, edge_namer, output_dir));
    }
    if (formats.find("gpa") != std::string::npos) {
      mapping_printers_.push_back(new MappingPrinterGPA(g, edge_namer, output_dir));
    }
    if (formats.find("fasta") != std::string::npos) {
      mapping_printers_.push_back(new MappingPrinterFasta(g, edge_namer, output_dir));
    }
  }

  void SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read)  {
    for (auto printer : mapping_printers_) {
      printer->SaveMapping(aligned_mappings, read);
    }
  }

  ~MappingPrinterHub() {
    for (auto printer : mapping_printers_) {
      delete printer;
    }
    mapping_printers_.clear();
  }

 private:
  std::vector<MappingPrinter*> mapping_printers_;
};

} // namespace sensitive_aligner
