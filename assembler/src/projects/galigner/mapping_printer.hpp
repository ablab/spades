#pragma once

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/pacbio/g_aligner.hpp"

namespace graph_aligner {

class MappingPrinter {
public:

    MappingPrinter(const debruijn_graph::ConjugateDeBruijnGraph &g, const std::string &output_file_prefix)
        : g_(g), output_file_prefix_(output_file_prefix)
    {}

    virtual void SaveMapping(const pacbio::OneReadMapping &aligned_mappings, const io::SingleRead &read) = 0;

    virtual ~MappingPrinter () {};

protected:
    const debruijn_graph::ConjugateDeBruijnGraph &g_;
    std::string output_file_prefix_;
    ofstream output_file_;
};

class MappingPrinterTSV: public MappingPrinter {
public:
    MappingPrinterTSV(const debruijn_graph::ConjugateDeBruijnGraph &g, const std::string &output_file_prefix)
        : MappingPrinter(g, output_file_prefix)
    {
        output_file_.open(output_file_prefix_ + ".tsv", std::ofstream::out);
        //tsv_file << "read\tstart_pos\tend_pos\tread_length\tgraph_path\tedges_lengths\tmapping\ted\n";
    }

    virtual void SaveMapping(const pacbio::OneReadMapping &aligned_mappings, const io::SingleRead &read);

    ~MappingPrinterTSV() {
        output_file_.close();
    }
};

class MappingPrinterGPA : public MappingPrinter {
public:
    MappingPrinterGPA(const debruijn_graph::ConjugateDeBruijnGraph &g, const std::string &output_file_prefix)
        : MappingPrinter(g, output_file_prefix)
    {
        output_file_.open(output_file_prefix_ + ".gpa", std::ofstream::out);
        output_file_ << "H\n";
    }

    std::string Print(map<string, string> &line);


    void CIGAR(std::string &read, std::string aligned, std::string &cigar, int &score);

    void DivideByEdgeCIGAR(string &read, string &aligned, std::vector<size_t> &edgeblocks, size_t start,
                        std::vector<string> &edgecigar,
                        std::vector<Range> &edge_initial_ranges, int &score);

    void MappedString(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappingpath, string &aligned, std::vector<size_t> &edgeblocks);

    void MappingOnRead(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappingpath, size_t &start, size_t &end);

    std::string SubRead(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappingpath, const io::SingleRead &read);

    virtual void SaveMapping(const pacbio::OneReadMapping &aligned_mappings, const io::SingleRead &read);

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

    void SaveMapping(const pacbio::OneReadMapping &aligned_mappings, const io::SingleRead &read) {
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
    vector<MappingPrinter*> mapping_printers_;
};

} // namespace graph_aligner