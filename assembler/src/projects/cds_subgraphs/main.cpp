//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "stop_condon_finder.hpp"
#include "subgraph_extraction.hpp"
#include "io/reads/file_reader.hpp"

#include "toolchain/edge_label_helper.hpp"
#include "toolchain/utils.hpp"
#include "toolchain/subgraph_utils.hpp"

#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include <clipp/clipp.h>
#include <string>
#include <set>
#include <numeric>
#include <sys/types.h>

using namespace debruijn_graph;

namespace cds_subgraphs {

struct PartialGeneInfo {
    size_t unitig_id;
    Range r;
    std::string gene_id;
    bool strand;
};

inline std::istream &operator>>(std::istream &is, PartialGeneInfo &pgi) {
    is >> pgi.unitig_id;
    is >> pgi.r.start_pos;
    is >> pgi.r.end_pos;
    is >> pgi.gene_id;
    std::string strand_symbol;
    is >> strand_symbol;
    if (strand_symbol == "+")
        pgi.strand = true;
    else if (strand_symbol == "-")
        pgi.strand = false;
    else CHECK_FATAL_ERROR(false, "Unsupported strand symbol");
    return is;
}

inline std::ostream &operator<<(std::ostream &os, const PartialGeneInfo &pgi) {
    os << "Unitig id: " << pgi.unitig_id << "; ";
    os << "Range: [" << pgi.r.start_pos << ", " << pgi.r.end_pos << "]; ";
    os << "Gene id:" << pgi.gene_id << "; ";
    os << "Strand: " << (pgi.strand ? "+" : "-");
    return os;
}

static Range ConjugateRange(const Graph &g, EdgeId e, const Range &r) {
    return Range(g.length(e) - r.end_pos, g.length(e) - r.start_pos);
}

using GeneInitSeq = std::multimap<std::string, std::string>;

static GeneInitSeq
PredictionsFromDescFile(const Graph &g,
                        const omnigraph::GraphElementFinder<Graph> &element_finder,
                        const std::string &desc_file) {
    fs::CheckFileExistenceFATAL(desc_file);
    std::ifstream descs(desc_file);
    std::string l;
    PartialGeneInfo info;

    size_t i = 1;
    GeneInitSeq starting_seqs;
    while (std::getline(descs, l)) {
        std::istringstream ss(l);
        ss >> info;

        INFO("Gene prediction #" << i << ": " << info);
        EdgeId e = element_finder.ReturnEdgeId(info.unitig_id);
        VERIFY(e != EdgeId());
        Range r = info.r;
        //exclusive right position
        r.end_pos += 1;
        //to 0-based coordinates
        r.shift(-1);
        VERIFY(r.size() > g.k());
        //to k+1-mer coordinates
//        r.end_pos -= g.k();

        if (!info.strand) {
            r = ConjugateRange(g, e, r);
            e = g.conjugate(e);
        }

        starting_seqs.insert(make_pair(info.gene_id, g.EdgeNucls(e).Subseq(r.start_pos, r.end_pos).str()));
        ++i;
    }
    return starting_seqs;
}

std::unordered_map<std::string, size_t>
CDSLengthsFromFile(const std::string &fn) {
    INFO("Parsing estimated CDS lengths from " << fn);
    fs::CheckFileExistenceFATAL(fn);
    std::unordered_map<std::string, size_t> answer;
    std::ifstream in(fn);
    std::string l;
    std::string gene_id;
    double est_len;
    while (std::getline(in, l)) {
        std::istringstream ss(l);
        ss >> gene_id;
        ss >> est_len;
        answer[gene_id] = RoundedProduct(3, est_len);
    }
    return answer;
}

static std::string GeneNameFromFasta(const std::string &header) {
    std::stringstream ss(header);
    std::string s;
    ss >> s;
    ss >> s;
    return s;
}

static GeneInitSeq PredictionsFromFastaFile(const std::string &fasta_fn) {
    fs::CheckFileExistenceFATAL(fasta_fn);
    io::FileReadStream gene_frags(fasta_fn);

    io::SingleRead gene_frag;

    size_t i = 1;
    GeneInitSeq starting_seqs;
    while (!gene_frags.eof()) {
        gene_frags >> gene_frag;
        INFO("Gene prediction #" << i << ": " << gene_frag.name());
        starting_seqs.insert(make_pair(GeneNameFromFasta(gene_frag.name()), gene_frag.GetSequenceString()));
        ++i;
    }
    return starting_seqs;
}

static void WriteComponent(const omnigraph::GraphComponent<Graph> &component, const std::string &prefix,
                           const std::set<GraphPos> &stop_codon_poss, const io::EdgeNamingF<Graph> &naming_f) {

    const auto &g = component.g();
    toolchain::WriteComponentWithDeadends(component, prefix, naming_f);

    INFO("Writing potential stop-codon positions to " << prefix << ".stops")
    std::ofstream stop_codon_os(prefix + ".stops");
    io::CanonicalEdgeHelper<Graph> canonical_helper(g, naming_f);
//1-based coordinate gives the start of the stop-codon
    for (GraphPos gpos : stop_codon_poss) {
        if (component.edges().count(gpos.first)) {
            stop_codon_os << canonical_helper.EdgeOrientationString(gpos.first, "\t")
                          << "\t" << (gpos.second + g.k() - 2) << "\n";
        } else {
            WARN("Earlier detected stop codon " << naming_f(g, gpos.first) << " "
                                                << gpos.second << " is outside the component");
        }
    }
}

void ExtractCDSSubgraphs(const GraphPack &gp,
                         const GeneInitSeq &starting_seqs,
                         const std::unordered_map<std::string, size_t> &cds_len_ests,
                         const io::EdgeNamingF<Graph> &edge_naming_f,
                         const std::string &out_folder) {
    const auto &graph = gp.get<Graph>();
    //fixme rename
    PartialGenePathProcessor partial_path_processor(graph, edge_naming_f);
    auto mapper = MapperInstance(gp);
    CDSSubgraphExtractor subgraph_extractor(graph, *mapper, partial_path_processor);
    INFO("Searching relevant subgraphs in parallel for all partial predictions");
    for (const auto &gene_id : utils::key_set(starting_seqs)) {
        INFO("Processing gene " << gene_id);
        INFO("Subgraphs extracted for all " << starting_seqs.count(gene_id)
                                            << " of its partial predictions will be united");
        std::set<EdgeId> edges;
        std::set<GraphPos> stop_codon_poss;
        for (const std::string &s : utils::get_all(starting_seqs, gene_id)) {
            auto gc = subgraph_extractor.ProcessPartialCDS(s,
                                                 utils::get(cds_len_ests, gene_id),
                                                 &stop_codon_poss);
            //TODO remove
            INFO("'Closing' gathered component");
            toolchain::ComponentExpander expander(graph);
            gc = expander.Expand(gc);
            utils::insert_all(edges, gc.edges());
        }
        toolchain::ComponentExpander expander(graph);
        auto component = expander.Expand(omnigraph::GraphComponent<Graph>::FromEdges(graph, edges.begin(),
                                                                          edges.end(), /*add conjugate*/true));

        if (component.e_size() > 0) {
            WriteComponent(component, out_folder + gene_id, stop_codon_poss, edge_naming_f);
        } else {
            INFO("Couldn't find a non-trivial component for gene " << gene_id);
        }
    }
}

void ParallelExtractCDSSubgraphs(const GraphPack &gp,
                                 const GeneInitSeq &starting_seqs,
                                 const std::unordered_map<std::string, size_t> &cds_len_ests,
                                 const io::EdgeNamingF<Graph> &edge_naming_f,
                                 const std::string &out_folder) {
    const auto &graph = gp.get<Graph>();
    //fixme rename
    PartialGenePathProcessor partial_path_processor(graph, edge_naming_f);
    auto mapper = MapperInstance(gp);
    CDSSubgraphExtractor subgraph_extractor(graph, *mapper, partial_path_processor);
    INFO("Searching relevant subgraphs in parallel for all partial predictions");
    std::vector<std::string> flattened_ids;
    std::vector<std::string> flattened_part_genes;
    for (const auto &id_part : starting_seqs) {
        flattened_ids.push_back(id_part.first);
        flattened_part_genes.push_back(id_part.second);
    }

    size_t n = flattened_ids.size();

    std::vector<std::set<EdgeId>> flattened_relevant_edges(n);
    std::vector<std::set<GraphPos>> flattened_stop_poss(n);

#pragma omp parallel for schedule(guided)
    for (size_t i = 0; i < n; ++i) {
        auto gc = subgraph_extractor.ProcessPartialCDS(flattened_part_genes[i],
                                                             utils::get(cds_len_ests, flattened_ids[i]),
                                                             &flattened_stop_poss[i]);
        //TODO remove
        INFO("'Closing' gathered component");
        toolchain::ComponentExpander expander(graph);
        gc = expander.Expand(gc);
        flattened_relevant_edges[i] = subgraph_extractor.ProcessPartialCDS(flattened_part_genes[i],
                                                                            utils::get(cds_len_ests, flattened_ids[i]),
                                                                            &flattened_stop_poss[i]).edges();
    }
    INFO("Done searching subgraphs");

    std::set<EdgeId> edges;
    std::set<GraphPos> stop_codon_poss;
    for (size_t i = 0; i < n; ++i) {
        const std::string &gene_id = flattened_ids[i];
        utils::insert_all(edges, flattened_relevant_edges[i]);
        utils::insert_all(stop_codon_poss, flattened_stop_poss[i]);

        if (i == n - 1 || flattened_ids[i + 1] != gene_id) {
            toolchain::ComponentExpander expander(graph);
            auto component = expander.Expand(omnigraph::GraphComponent<Graph>::FromEdges(graph, edges.begin(),
                                                                              edges.end(), /*add conjugate*/true));

            if (component.e_size() > 0) {
                WriteComponent(component, out_folder + gene_id, stop_codon_poss, edge_naming_f);
            } else {
                INFO("Couldn't find a non-trivial component for gene " << gene_id);
            }
            edges.clear();
            stop_codon_poss.clear();
        }
    }
}

}

struct gcfg {
    gcfg()
            : k(0),
              nthreads(omp_get_max_threads() / 2 + 1)
    {}

    unsigned k;
    std::string graph;
    std::string tmpdir;
    std::string outdir;
    std::string genes_desc;
    std::string genes_seq;
    std::string cds_len_fn;
    unsigned nthreads;
};

static void process_cmdline(int argc, char **argv, gcfg &cfg) {
    using namespace clipp;

    auto cli = (
            (required("-o", "--output-folder") & value("dir", cfg.outdir)) % "output folder to use for GFA files",
                    one_of((option("--part-desc") & value("file", cfg.genes_desc)) % "file with partial genes description (.gff)",
                           (option("--part-seq") & value("file", cfg.genes_seq)) % "file with partial genes sequences (.fasta)"),
                    (required("--graph") & value("graph", cfg.graph)) % "In GFA (ending with .gfa) or prefix to SPAdes graph pack",
                    (required("--cds-len-est") & value("file", cfg.cds_len_fn)) % "file with cds length estimamtes",
                    (required("-k") & integer("value", cfg.k)) % "k-mer length to use",
                    (option("-t", "--threads") & integer("value", cfg.nthreads)) % "# of threads to use (default: max_threads / 2)",
                    (option("--tmpdir") & value("dir", cfg.tmpdir)) % "scratch directory to use (default: <outdir>/tmp)"
    );

    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cout << make_man_page(cli, argv[0]);
        exit(1);
    }
}

int main(int argc, char** argv) {
    utils::segfault_handler sh;
    gcfg cfg;

    process_cmdline(argc, argv, cfg);

    toolchain::create_console_logger(/*logging::L_TRACE*/);
    START_BANNER("Extracting relevant subgraphs for (partial) predicted CDS");

    try {
        unsigned nthreads = cfg.nthreads;
        unsigned k = cfg.k;
        std::string out_folder = cfg.outdir + "/";
        fs::make_dirs(out_folder);

        INFO("K-mer length set to " << k);

        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);
        INFO("# of threads to use: " << nthreads);

        std::string tmpdir = cfg.tmpdir.empty() ? out_folder + "tmp" : cfg.tmpdir;
        fs::make_dirs(tmpdir);
        debruijn_graph::GraphPack gp(k, tmpdir, 0);
        
        const auto &graph = gp.get<Graph>();
        omnigraph::GraphElementFinder<Graph> element_finder(graph);
        INFO("Loading de Bruijn graph from " << cfg.graph);
        gp.get_mutable<KmerMapper<Graph>>().Attach(); // TODO unnecessary
        io::EdgeLabelHelper<Graph> label_helper(element_finder,
                toolchain::LoadGraph(gp, cfg.graph));

        gp.EnsureBasicMapping();

        VERIFY(cfg.genes_desc.empty() != cfg.genes_seq.empty());

        auto starting_seqs = cfg.genes_desc.empty() ?
                              cds_subgraphs::PredictionsFromFastaFile(cfg.genes_seq) :
                             cds_subgraphs::PredictionsFromDescFile(graph, element_finder, cfg.genes_desc);

        auto cds_len_ests = cds_subgraphs::CDSLengthsFromFile(cfg.cds_len_fn);

        static const bool parallel = false;
        if (parallel) {
            //Experimental parallel mode
            cds_subgraphs::ParallelExtractCDSSubgraphs(gp, starting_seqs, cds_len_ests,
                              label_helper.edge_naming_f(), out_folder);
        } else {
            cds_subgraphs::ExtractCDSSubgraphs(gp, starting_seqs, cds_len_ests,
                              label_helper.edge_naming_f(), out_folder);
        }

        INFO("Done");
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }
}
