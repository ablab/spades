#pragma once

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/pacbio/g_aligner.hpp"

namespace sensitive_aligner {

class MappingPrinter {
  public:

    MappingPrinter(const debruijn_graph::ConjugateDeBruijnGraph &g, const std::string &output_file_prefix)
        : g_(g), output_file_prefix_(output_file_prefix)
    {}

    virtual void SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) = 0;

    virtual ~MappingPrinter () {};

  protected:
    const debruijn_graph::ConjugateDeBruijnGraph &g_;
    std::string output_file_prefix_;
    std::ofstream output_file_;
};

class MappingPrinterTSV: public MappingPrinter {
  public:
    MappingPrinterTSV(const debruijn_graph::ConjugateDeBruijnGraph &g, const std::string &output_file_prefix)
        : MappingPrinter(g, output_file_prefix) {
        output_file_.open(output_file_prefix_ + ".tsv", std::ofstream::out);
        //tsv_file << "read\tstart_pos\tend_pos\tread_length\tgraph_path\tedges_lengths\tmapping\ted\n";
    }

    virtual void SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read);

    ~MappingPrinterTSV() {
        output_file_.close();
    }
};

class MappingPrinterGPA : public MappingPrinter {
  public:
    MappingPrinterGPA(const debruijn_graph::ConjugateDeBruijnGraph &g, const std::string &output_file_prefix)
        : MappingPrinter(g, output_file_prefix) {
        output_file_.open(output_file_prefix_ + ".gpa", std::ofstream::out);
        output_file_ << "H\n";
    }

    std::string Print(std::map<std::string, std::string> &line);


    std::string getCigar(const std::string &read, const std::string &aligned);

    void getEdgeCigar(const std::string &subread, const std::string &path_seq, const std::vector<size_t> &edgeblocks,
                      std::vector<std::string> &edgecigar, std::vector<Range> &edgeranges);

    void getPath(const std::vector<debruijn_graph::EdgeId> &path,
                 const PathRange &path_range,
                 std::string &aligned, std::vector<size_t> &edgeblocks);

    std::string getSubread(const Sequence &read, const PathRange &path_range);

    std::string formGPAOutput(const io::SingleRead &read,
                         const std::vector<debruijn_graph::EdgeId> &path,
                         const std::vector<std::string> &edgecigar,
                         const std::vector<Range> &edgeranges,
                         int &nameIndex, const PathRange &path_range);

    virtual void SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read);

    ~MappingPrinterGPA() {
        output_file_.close();
    }

};

class MappingPrinterHub {
  public:
    MappingPrinterHub(const debruijn_graph::ConjugateDeBruijnGraph &g, const std::string &output_file_prefix, const std::string formats) {
        if (formats.find("tsv") != std::string::npos) {
            mapping_printers_.push_back(new MappingPrinterTSV(g, output_file_prefix));
        }
        if (formats.find("gpa") != std::string::npos) {
            mapping_printers_.push_back(new MappingPrinterGPA(g, output_file_prefix));
        }
    }

    void SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) {
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