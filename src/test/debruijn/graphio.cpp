//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graphio.hpp"

#include "alignment/edge_index.hpp"
#include "alignment/kmer_mapper.hpp"
#include "assembly_graph/core/construction_helper.hpp"
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "paired_info/index_point.hpp"
#include "pipeline/graph_pack.hpp"

namespace debruijn_graph {

inline void DeserializePoint(FILE* file, size_t& e1, size_t& e2, omnigraph::de::RawPoint &p) {
    float unused;
    size_t read_count = fscanf(file, "%zu %zu %f %f %f .\n", &e1, &e2,
                               (float *)&p.d, (float *)&p.weight, (float *)&unused);
    VERIFY(read_count == 5);

}

inline void DeserializePoint(FILE* file, size_t& e1, size_t& e2, omnigraph::de::Point &p) {
    size_t read_count = fscanf(file, "%zu %zu %f %f %f .\n", &e1, &e2,
                               (float *)&p.d, (float *)&p.weight, (float *)&p.var);
    VERIFY(read_count == 5);
}

class LegacyTextIO {

public:
    bool LoadGraph(const std::string &file_name, Graph &graph) {
        INFO("Trying to read conjugate de bruijn graph from " << file_name << ".grp");
        FILE* file = fopen((file_name + ".grp").c_str(), "r");
        if (file == NULL) {
            ERROR("Couldn't find file " << (file_name + ".grp"));
            return false;
        }
        FILE* sequence_file = fopen((file_name + ".sqn").c_str(), "r");
        if (sequence_file == NULL) {
            ERROR("Couldn't find file " << (file_name + ".sqn"));
            return false;
        }
        INFO("Reading conjugate de bruijn  graph from " << file_name << " started");
        typename Graph::HelperT helper(graph);
        size_t vertex_count;
        size_t edge_count;
        int flag = fscanf(file, "%zu %zu \n", &vertex_count, &edge_count);
        VERIFY(flag == 2);
        for (size_t i = 0; i < vertex_count; i++) {
            size_t vid, conj_vid;
            flag = fscanf(file, "Vertex %zu ~ %zu .\n", &vid, &conj_vid);
            TRACE("Vertex " << vid << " ~ " << conj_vid << " .");
            if (flag != 2) {
                ERROR("Malformed saves");
                return false;
            }

            if (!graph.contains(VertexId(vid))) {
                size_t needed = std::max(vid, conj_vid);
                if (graph.vreserved() <= needed)
                    graph.vreserve(needed * 2 + 1);
                VertexId new_id = graph.AddVertex(typename Graph::VertexData(graph.k()), vid, conj_vid);
                VERIFY(new_id = vid);
                VERIFY(graph.conjugate(new_id) = conj_vid);
            }
        }

        if (!edge_count) {
            fclose(file);
            fclose(sequence_file);
        }

        char first_char = (char) getc(sequence_file);
        VERIFY(!ferror(sequence_file));
        ungetc(first_char, sequence_file);
        bool fasta = (first_char == '>'); // if it's not fasta, then it's old .sqn

        if (!fasta) {
            size_t tmp_edge_count;
            flag = fscanf(sequence_file, "%zu", &tmp_edge_count);
            if (flag != 1) {
                ERROR("Malformed saves");
                return false;
            }
            if (edge_count != tmp_edge_count) {
                ERROR("Wrong edge count");
                return false;
            }
        }

        static const size_t longstring_size = 1000500; //Enough for tests
        char longstring[longstring_size];
        for (size_t i = 0; i < edge_count; i++) {
            size_t eid, start_id, fin_id, length, conj_eid;
            flag = fscanf(file, "Edge %zu : %zu -> %zu, l = %zu ~ %zu .\n",
                          &eid, &start_id, &fin_id, &length, &conj_eid);
            if (flag != 5) {
                ERROR("Malformed saves");
                return false;
            }
            VERIFY(length < longstring_size);
            size_t eid_tmp;
            if (fasta) {
                flag = fscanf(sequence_file, ">%zu\n%s\n", &eid_tmp, longstring);
            } else {
                flag = fscanf(sequence_file, "%zu %s .", &eid_tmp, longstring);
            }
            VERIFY(eid == eid_tmp);
            if (flag != 2) {
                ERROR("Malformed saves");
                return false;
            }
            TRACE("Edge " << eid << " : " << start_id << " -> "
                          << fin_id << " l = " << length << " ~ " << conj_eid);
            if (!graph.contains(EdgeId(eid))) {
                size_t needed = std::min(eid, conj_eid);
                if (graph.ereserved() <= needed)
                    graph.ereserve(needed * 2 + 1);
                Sequence tmp(longstring);
                EdgeId new_id = graph.AddEdge(start_id, fin_id, Graph::EdgeData(tmp), eid, conj_eid);
                VERIFY(new_id == eid);
                VERIFY(graph.conjugate(new_id) == conj_eid);
            }
        }

        fclose(file);
        fclose(sequence_file);
        return true;
    }

    bool LoadCoverage(const std::string& file_name, omnigraph::CoverageIndex<Graph> &cov) {
        INFO("Reading coverage from " << file_name);
        std::ifstream in(file_name + ".cvr");
        if (!in.good()) {
            ERROR("Failed to open: " << (file_name + ".cvr"));
            return false;
        }

        LoadEdgeAssociatedInfo(in, cov);
        return true;
    }

    bool LoadFlankingCoverage(const std::string &file_name, omnigraph::FlankingCoverage<Graph> &flanking_cov) {
        if (!std::filesystem::exists(file_name + ".flcvr")) {
            INFO("Flanking coverage saves are missing");
            return false;
        }
        INFO("Reading flanking coverage from " << file_name);
        std::ifstream in(file_name + ".flcvr");
        LoadEdgeAssociatedInfo(in, flanking_cov);
        return true;
    }

    bool LoadKmerMapper(const std::string &file_name, KmerMapper<Graph> &kmer_mapper) {
        kmer_mapper.clear();
        std::ifstream file;
        file.open((file_name + ".kmm").c_str(),
                  std::ios_base::binary | std::ios_base::in);
        if (!file.is_open()) {
            return false;
        }
        INFO("Reading kmer mapper, " << file_name <<" started");

        uint32_t k_;
        file.read((char *) &k_, sizeof(uint32_t));

        VERIFY_MSG(k_ == kmer_mapper.k(), "Cannot read kmer mapper, different Ks");
        kmer_mapper.BinRead(file);

        file.close();
        return true;
    }

    bool LoadPositions(const std::string &file_name, omnigraph::EdgesPositionHandler<Graph> &edge_pos) {
        FILE* file = fopen((file_name + ".pos").c_str(), "r");
        if (file == NULL) {
            INFO("No positions were saved");
            return false;
        }
        VERIFY(!edge_pos.IsAttached());
        edge_pos.Attach();
        INFO("Reading edges positions, " << file_name << " started");
        VERIFY(file != NULL);
        size_t pos_count;
        int read_count = fscanf(file, "%zu\n", &pos_count);
        VERIFY(read_count == 1);
        for (size_t i = 0; i < pos_count; i++) {
            size_t eid, pos_info_count;
            char contigId[500];
            char cur_str[500];
            read_count = fscanf(file, "%zu %zu\n", &eid, &pos_info_count);
            VERIFY(read_count == 2);
            for (size_t j = 0; j < pos_info_count; j++) {
                int start_pos, end_pos;
                int m_start_pos, m_end_pos;
                read_count = fscanf(file, "%[^\n]s", cur_str);
                read_count = fscanf(file, "\n");
                read_count = sscanf(cur_str, "%s [%d - %d] --> [%d - %d]", contigId,
                                    &start_pos, &end_pos, &m_start_pos, &m_end_pos);
                VERIFY(read_count == 5);
                edge_pos.AddEdgePosition(eid, contigId, start_pos - 1, end_pos, m_start_pos - 1, m_end_pos);
            }
        }
        fclose(file);
        return true;
    }

    bool LoadBasicGraph(const std::string &file_name, Graph &g) {
        if (!LoadGraph(file_name, g))
            return false;
        if (!LoadCoverage(file_name, g.coverage_index()))
            return false;
        return true;
    }

    bool LoadGraphPack(const std::string &file_name, graph_pack::GraphPack &gp) {
        auto &graph = gp.get_mutable<Graph>();
        auto &index = gp.get_mutable<EdgeIndex<Graph>>();
        auto &edge_pos = gp.get_mutable<omnigraph::EdgesPositionHandler<Graph>>();
        auto &flanking_cov = gp.get_mutable<omnigraph::FlankingCoverage<Graph>>();
        auto &kmer_mapper = gp.get_mutable<KmerMapper<Graph>>();

        if (!LoadBasicGraph(file_name, graph))
            return false;
        index.Attach();
        index.Refill();

        LoadPositions(file_name, edge_pos);
        //load kmer_mapper only if needed
        if (kmer_mapper.IsAttached())
            if (!LoadKmerMapper(file_name, kmer_mapper)) {
                WARN("Cannot load kmer_mapper, information on projected kmers will be missed");
            }
        if (!LoadFlankingCoverage(file_name, flanking_cov)) {
            FATAL_ERROR("Cannot load flanking coverage");
        }
        return true;
    }

private:
    template<class T>
    void LoadEdgeAssociatedInfo(std::istream &in, T &t) const {
        size_t cnt;
        in >> cnt;
        for (size_t i = 0 ; i < cnt; ++i) {
            size_t eid;
            in >> eid;
            t.Load(eid, in);
            std::string delim;
            in >> delim;
            VERIFY(delim == ".");
        }
    }
};

//Legacy wrappers
namespace graphio {

bool ScanBasicGraph(const std::string &file_name, debruijn_graph::Graph &g) {
    return LegacyTextIO().LoadBasicGraph(file_name, g);
}

bool ScanGraphPack(const std::string &file_name, graph_pack::GraphPack &gp) {
    return LegacyTextIO().LoadGraphPack(file_name, gp);
}

}

} // Namespace debruijn_graph
