//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/construction_helper.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
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
    void LoadGraph(const std::string &file_name, Graph &graph) {
        INFO("Trying to read conjugate de bruijn graph from " << file_name << ".grp");
        FILE* file = fopen((file_name + ".grp").c_str(), "r");
        VERIFY_MSG(file != NULL, "Couldn't find file " << (file_name + ".grp"));
        FILE* sequence_file = fopen((file_name + ".sqn").c_str(), "r");
        VERIFY_MSG(file != NULL, "Couldn't find file " << (file_name + ".sqn"));
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
            VERIFY(flag == 2);

            if (!graph.contains(VertexId(vid))) {
                size_t needed = std::max(vid, conj_vid);
                if (graph.vreserved() <= needed)
                    graph.vreserve(needed * 2 + 1);
                VertexId new_id = graph.AddVertex(typename Graph::VertexData(), vid, conj_vid);
                VERIFY(new_id = vid);
                VERIFY(graph.conjugate(new_id) = conj_vid);
            }
        }

        if (!edge_count) {
            fclose(file);
            fclose(sequence_file);
            return;
        }

        char first_char = (char) getc(sequence_file);
        VERIFY(!ferror(sequence_file));
        ungetc(first_char, sequence_file);
        bool fasta = (first_char == '>'); // if it's not fasta, then it's old .sqn

        if (!fasta) {
            size_t tmp_edge_count;
            flag = fscanf(sequence_file, "%zu", &tmp_edge_count);
            VERIFY(flag == 1);
            VERIFY(edge_count == tmp_edge_count);
        }

        static const size_t longstring_size = 1000500; //Enough for tests
        char longstring[longstring_size];
        for (size_t i = 0; i < edge_count; i++) {
            size_t eid, start_id, fin_id, length, conj_eid;
            flag = fscanf(file, "Edge %zu : %zu -> %zu, l = %zu ~ %zu .\n",
                          &eid, &start_id, &fin_id, &length, &conj_eid);
            VERIFY(flag == 5);
            VERIFY(length < longstring_size);
            size_t eid_tmp;
            if (fasta) {
                flag = fscanf(sequence_file, ">%zu\n%s\n", &eid_tmp, longstring);
            } else {
                flag = fscanf(sequence_file, "%zu %s .", &eid_tmp, longstring);
            }
            VERIFY(eid == eid_tmp);
            VERIFY(flag == 2);
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
    }

    void LoadCoverage(const std::string& file_name, omnigraph::CoverageIndex<Graph> &cov) {
        INFO("Reading coverage from " << file_name);
        std::ifstream in(file_name + ".cvr");
        LoadEdgeAssociatedInfo(in, cov);
    }

    bool LoadFlankingCoverage(const std::string &file_name, omnigraph::FlankingCoverage<Graph> &flanking_cov) {
        if (!fs::FileExists(file_name + ".flcvr")) {
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

    template<typename Index>
    void LoadPaired(const std::string &file_name, Index &paired_index, bool force_exists = true) {
        FILE* file = fopen((file_name + ".prd").c_str(), "r");
        INFO((file_name + ".prd"));
        if (force_exists) {
            VERIFY(file != NULL);
        } else if (file == NULL) {
            INFO("Paired info not found, skipping");
            return;
        }
        INFO("Reading paired info from " << file_name << " started");

        size_t paired_count;
        int read_count = fscanf(file, "%zu \n", &paired_count);
        VERIFY(read_count == 1);
        while (!feof(file)) {
            size_t e1, e2;

            typename Index::Point point;
            DeserializePoint(file, e1, e2, point);

            TRACE(e1 << " " << e2 << " " << point);
            //Need to prevent doubling of self-conjugate edge pairs
            //Their weight would be always even, so we don't lose precision
            auto ep = std::make_pair(e1, e2);
            if (ep == paired_index.ConjugatePair(ep))
                point.weight = math::round(point.weight / 2);
            paired_index.Add(e1, e2, point);
        }
        DEBUG("PII SIZE " << paired_index.size());
        fclose(file);
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

    void LoadBasicGraph(const std::string &file_name, Graph &g) {
        LoadGraph(file_name, g);
        LoadCoverage(file_name, g.coverage_index());
    }

    void LoadGraphPack(const std::string &file_name, GraphPack &gp) {
        auto &graph = gp.get_mutable<Graph>();
        auto &index = gp.get_mutable<EdgeIndex<Graph>>();
        auto &edge_pos = gp.get_mutable<omnigraph::EdgesPositionHandler<Graph>>();
        auto &flanking_cov = gp.get_mutable<omnigraph::FlankingCoverage<Graph>>();
        auto &kmer_mapper = gp.get_mutable<KmerMapper<Graph>>();

        LoadBasicGraph(file_name, graph);

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

void ScanBasicGraph(const std::string &file_name, Graph &g) {
    LegacyTextIO().LoadBasicGraph(file_name, g);
}

void ScanGraphPack(const std::string &file_name, GraphPack &gp) {
    LegacyTextIO().LoadGraphPack(file_name, gp);
}

}

} // namespace debruijn_graph
