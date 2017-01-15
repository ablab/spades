//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/standard_base.hpp"

#include "assembly_graph/handlers/id_track_handler.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "assembly_graph/components/graph_component.hpp"

#include "paired_info/paired_info.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"
#include "modules/alignment/long_read_storage.hpp"

#include "assembly_graph/core/order_and_law.hpp"

#include <cmath>
#include <set>
#include <map>
#include <algorithm>
#include <fstream>
#include <cstdio>

namespace debruijn_graph {

namespace graphio {

using namespace omnigraph;
using namespace omnigraph::de;
//todo think of inner namespace

template<class KmerMapper>
void SaveKmerMapper(const string& file_name,
                    const KmerMapper& mapper) {
    std::ofstream file;
    file.open((file_name + ".kmm").c_str(),
              std::ios_base::binary | std::ios_base::out);
    DEBUG("Saving kmer mapper, " << file_name <<" created");
    VERIFY(file.is_open());

    uint32_t k = (uint32_t) mapper.k();
    file.write((char *) &k, sizeof(uint32_t));
    mapper.BinWrite(file);

    file.close();
    DEBUG("kmer mapper saved ")
}

template<class KmerMapper>
bool LoadKmerMapper(const string& file_name,
                    KmerMapper& kmer_mapper) {
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

template<class EdgeIndex>
void SaveEdgeIndex(const std::string& file_name,
                   const EdgeIndex& index) {
    std::ofstream file;
    file.open((file_name + ".kmidx").c_str(),
              std::ios_base::binary | std::ios_base::out);
    DEBUG("Saving kmer index, " << file_name <<" created");
    VERIFY(file.is_open());

    uint32_t k_ = index.k();
    file.write((char *) &k_, sizeof(uint32_t));
    index.BinWrite(file);

    file.close();
    DEBUG("index saved ")
}

template<class EdgeIndex>
bool LoadEdgeIndex(const std::string& file_name,
                   EdgeIndex& index) {
    std::ifstream file;
    file.open((file_name + ".kmidx").c_str(),
              std::ios_base::binary | std::ios_base::in);
    INFO("Reading kmer index, " << file_name <<" started");
    if (!file.is_open())
        return false;

    uint32_t k_;
    file.read((char *) &k_, sizeof(uint32_t));
    VERIFY_MSG(k_ == index.k(), "Cannot read edge index, different Ks:");

    index.BinRead(file, file_name + ".kmidx");

    file.close();

    return true;
}

inline
void SaveMapCoverage(const std::string& path, const std::map<int, int>& data ) {
    std::ofstream outFile;
    outFile.open(path.c_str());

    INFO("Saving detailed coverage in file " << path <<" started");
    outFile << data.size() << "\n";
    for (auto dataIterator = data.begin(); dataIterator != data.end(); ++dataIterator){
        outFile << dataIterator->first << " " << dataIterator->second << " .\n";
    }
}

template<class KmerIndex>
void SaveDetailCoverage(const std::string& pathInCov, const std::string& pathOutCov, const KmerIndex& index ) {
    SaveMapCoverage(pathInCov, index.inCoverage);
    SaveMapCoverage(pathOutCov, index.outCoverage);
}

inline void SerializePoint(FILE* file, size_t e1, size_t e2, const RawPoint &p) {
    fprintf(file, "%zu %zu %.2f %.2f 0.00 .\n", e1, e2, (double)p.d, (double)p.weight);
}

inline void SerializePoint(FILE* file, size_t e1, size_t e2, const Point &p) {
    fprintf(file, "%zu %zu %.2f %.2f %.2f .\n", e1, e2, (double)p.d, (double)p.weight, (double)p.var);
}

inline void DeserializePoint(FILE* file, size_t& e1, size_t& e2, RawPoint &p) {
    float unused;
    size_t read_count = fscanf(file, "%zu %zu %f %f %f .\n", &e1, &e2,
        (float *)&p.d, (float *)&p.weight, (float *)&unused);
    VERIFY(read_count == 5);

}

inline void DeserializePoint(FILE* file, size_t& e1, size_t& e2, Point &p) {
    size_t read_count = fscanf(file, "%zu %zu %f %f %f .\n", &e1, &e2,
                               (float *)&p.d, (float *)&p.weight, (float *)&p.var);
    VERIFY(read_count == 5);
}


template<class Graph>
class DataPrinter {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    //todo reduce duplication
    template<class T>
    void SaveEdgeAssociatedInfo(std::function<T (EdgeId)> access_f, ostream& out) const {
        out << component_.e_size() << endl;
        for (auto iter = component_.e_begin(); iter != component_.e_end(); ++iter) {
            EdgeId e = *iter;
            //todo fixme currently matches old format .cvr format
            out << e.int_id()/* << endl*/;
            out << " " << access_f(e) << " ." << endl;
        }
    }

//    template<class C>
//    void SaveEdgeAssociatedInfo(const C& c, ostream& out) const {
//        SaveEdgeAssociatedInfo<decltype(C::operator[])>(boost::bind(&C::operator[], c, _1), out);
//    }

    template<class C>
    void SaveEdgeAssociatedInfo(const C& c, ostream& out) const {
        out << component_.e_size() << endl;
        for (auto iter = component_.e_begin(); iter != component_.e_end(); ++iter) {
            EdgeId e = *iter;
            //todo fixme currently matches old format .cvr format
            out << e.int_id()/* << endl*/;
            out << " ";
            c.Save(e, out);
            out << " ." << endl;
        }
    }

  public:

    void SaveGraph(const string& file_name) const {
        FILE* gid_file = fopen((file_name + ".gid").c_str(), "w");
        size_t max_id = this->component().g().GetGraphIdDistributor().GetMax();
        fprintf(gid_file, "%zu\n", max_id);
        fclose(gid_file);
        FILE* file = fopen((file_name + ".grp").c_str(), "w");
        DEBUG("Graph saving to " << file_name << " started");
        VERIFY_MSG(file != NULL,
                   "Couldn't open file " << (file_name + ".grp") << " on write");
        size_t vertex_count = component_.v_size();
        size_t edge_count = component_.e_size();
        fprintf(file, "%zu %zu \n", vertex_count, edge_count);
        for (auto iter = component_.v_begin(); iter != component_.v_end(); ++iter) {
            Save(file, *iter);
        }

        fprintf(file, "\n");

        for (auto iter = component_.e_begin(); iter != component_.e_end(); ++iter) {
            Save(file, *iter);
        }
        DEBUG("Graph saving to " << file_name << " finished");

        fclose(file);
    }

    void SaveEdgeSequences(const string& file_name) const {
        ofstream out(file_name + ".sqn");
        //todo switch to general function after its switching to fasta
        //SaveEdgeAssociatedInfo<Sequence>(boost::bind(&Graph::EdgeNucls, component_.g(), _1), out);
        DEBUG("Saving sequences, " << file_name <<" created");
        for (auto iter = component_.e_begin(); iter != component_.e_end(); ++iter) {
            EdgeId e = *iter;
            out << ">" << e.int_id() << endl;
            out << component_.g().EdgeNucls(e) << endl;
        }
    }

    void SaveCoverage(const string& file_name) const {
        ofstream out(file_name + ".cvr");
        DEBUG("Saving coverage, " << file_name <<" created");
        SaveEdgeAssociatedInfo(component_.g().coverage_index(), out);
    }

    void SaveFlankingCoverage(const string& file_name, const FlankingCoverage<Graph>& flanking_cov) const {
        ofstream out(file_name + ".flcvr");
        DEBUG("Saving flanking coverage, " << file_name <<" created");
        SaveEdgeAssociatedInfo(flanking_cov, out);
    }

    template<class Index>
    void SavePaired(const string& file_name,
                    Index const& paired_index) const {
        FILE* file = fopen((file_name + ".prd").c_str(), "w");
        DEBUG("Saving paired info, " << file_name <<" created");
        VERIFY(file != NULL);

        size_t comp_size = 0;
        for (auto I = component_.e_begin(), E = component_.e_end(); I != E; ++I) {
            EdgeId e1 = *I;
            auto inner_map = paired_index.GetHalf(e1);
            for (auto entry : inner_map) {
                if (component_.contains(entry.first)) { // if the second edge also lies in the same component
                    comp_size += entry.second.size();
                    continue;
                }
            }
        }

        fprintf(file, "%zu\n", comp_size);

        for (auto I = component_.e_begin(), E = component_.e_end(); I != E; ++I) {
            EdgeId e1 = *I;
            const auto& inner_map = paired_index.GetHalf(e1);
            std::map<typename Graph::EdgeId, typename Index::HistProxy> ordermap(inner_map.begin(), inner_map.end());
            for (auto entry : ordermap) {
                EdgeId e2 = entry.first;
                if (component_.contains(e2))
                    for (auto point : entry.second)
                        SerializePoint(file, e1.int_id(), e2.int_id(), point);
            }
        }

        fclose(file);
    }

    void SavePositions(const string& file_name,
                       EdgesPositionHandler<Graph> const& ref_pos) const {
        ofstream file((file_name + ".pos").c_str());
        DEBUG("Saving edges positions, " << file_name << " created");
        VERIFY(file.is_open());
        file << component_.e_size() << endl;
        for (auto it = component_.e_begin(); it != component_.e_end(); ++it) {
            vector<omnigraph::EdgePosition> pos_it = ref_pos.GetEdgePositions(*it);
            file << it->int_id() << " " << pos_it.size() << endl;
            for (size_t i = 0; i < pos_it.size(); i++) {
                file << "    " << pos_it[i].contigId << " " << pos_it[i].mr << endl;
            }
        }
    }

  private:
    void Save(FILE* file, EdgeId eid) const {
        fprintf(file, "%s\n", ToPrint(eid).c_str());
    }

    void Save(FILE* file, VertexId vid) const {
        fprintf(file, "%s\n", ToPrint(vid).c_str());
    }

    const GraphComponent<Graph> component_;

    virtual std::string ToPrint(VertexId v) const = 0;
    virtual std::string ToPrint(EdgeId e) const = 0;

  protected:

    //todo optimize component copy
//    DataPrinter(const GraphComponent<Graph>& component) :
//            component_(component) {
//    }

    DataPrinter(GraphComponent<Graph>&& component) :
            component_(std::move(component)) {
    }

    const GraphComponent<Graph>& component() const {
        return component_;
    }

  public:
    virtual ~DataPrinter() {
    }
};

template<class Graph>
class ConjugateDataPrinter: public DataPrinter<Graph> {
    typedef DataPrinter<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
  public:
    ConjugateDataPrinter(Graph const& g) :
            base(GraphComponent<Graph>::WholeGraph(g)) {
    }

    ConjugateDataPrinter(const GraphComponent<Graph>& graph_component) :
            base(GraphComponent<Graph>(graph_component, true)) {
    }

    template<class VertexIt>
    ConjugateDataPrinter(const Graph& g, VertexIt begin, VertexIt end) :
            base(GraphComponent<Graph>::FromVertices(g, begin, end, true)) {
    }

    std::string ToPrint(VertexId v) const {
        stringstream ss;
        ss
                << "Vertex "
                << v.int_id()
                << " ~ "
                << this->component().g().conjugate(v).int_id() << " .";
        return ss.str();
    }

    std::string ToPrint(EdgeId e) const {
        stringstream ss;
        ss
                << "Edge "
                << e.int_id()
                << " : "
                << this->component().g().EdgeStart(e).int_id()
                << " -> "
                << this->component().g().EdgeEnd(e).int_id()
                << ", l = "
                << this->component().g().length(e)
                << " ~ "
                << this->component().g().conjugate(e).int_id() << " .";
        return ss.str();
    }

};

template<class Graph>
class DataScanner {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    template<class T>
    void LoadEdgeAssociatedInfo(std::function<void (EdgeId, T)> setting_f, istream& in) const {
        size_t cnt;
        in >> cnt;
        for (size_t i = 0 ; i < cnt; ++i) {
            size_t edge_id;
            T t;
            string delim;
            in >> edge_id;
            in >> t;
            in >> delim;
            VERIFY(delim == ".");
            VERIFY(this->edge_id_map().find(edge_id) != this->edge_id_map().end());
            setting_f(this->edge_id_map()[edge_id], t);
        }
    }

    template<class T>
    void LoadEdgeAssociatedInfo(T& t, istream& in) const {
        size_t cnt;
        in >> cnt;
        for (size_t i = 0 ; i < cnt; ++i) {
            size_t edge_id;
            in >> edge_id;
            VERIFY(this->edge_id_map().find(edge_id) != this->edge_id_map().end());
            EdgeId eid = this->edge_id_map().find(edge_id)->second;
            t.Load(eid, in);
            string delim;
            in >> delim;
            VERIFY(delim == ".");
        }
    }

//    template<class C>
//    void LoadEdgeAssociatedInfo(const C& c, ostream& out) const {
//        SaveEdgeAssociatedInfo<decltype(C::operator[])>(boost::bind(&C::operator[], c, _1), out);
//    }

  public:
    virtual void LoadGraph(const string& file_name) = 0;

    void LoadCoverage(const string& file_name) {
        INFO("Reading coverage from " << file_name);
        ifstream in(file_name + ".cvr");
        LoadEdgeAssociatedInfo(g_.coverage_index(), in);
    }

    bool LoadFlankingCoverage(const string& file_name, FlankingCoverage<Graph>& flanking_cov) {
        if (!path::FileExists(file_name + ".flcvr")) {
            INFO("Flanking coverage saves are absent");
            return false;
        }
        INFO("Reading flanking coverage from " << file_name);
        ifstream in(file_name + ".flcvr");
        LoadEdgeAssociatedInfo(flanking_cov, in);
        return true;
    }

    template<typename Index>
    void LoadPaired(const string& file_name,
                    Index& paired_index,
                    bool force_exists = true) {
        typedef typename Graph::EdgeId EdgeId;
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
            size_t first_real_id, second_real_id;

            typename Index::Point point;
            DeserializePoint(file, first_real_id, second_real_id, point);

            TRACE(first_real_id << " " << second_real_id << " " << point);
            VERIFY(this->edge_id_map().find(first_real_id) != this->edge_id_map().end())
            EdgeId e1 = this->edge_id_map()[first_real_id];
            EdgeId e2 = this->edge_id_map()[second_real_id];
            if (e1 == EdgeId(NULL) || e2 == EdgeId(NULL))
                continue;
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

    bool LoadPositions(const string& file_name,
                       EdgesPositionHandler<Graph>& edge_pos) {
        FILE* file = fopen((file_name + ".pos").c_str(), "r");
        if (file == NULL) {
            INFO("No positions were saved");
            return false;
        }
        VERIFY(!edge_pos.IsAttached());
        edge_pos.Attach();
        INFO("Reading edges positions, " << file_name <<" started");
        VERIFY(file != NULL);
        size_t pos_count;
        int read_count = fscanf(file, "%zu\n", &pos_count);
        VERIFY(read_count == 1);
        for (size_t i = 0; i < pos_count; i++) {
            size_t edge_real_id, pos_info_count;
            char contigId[500];
            char cur_str[500];
            read_count = fscanf(file, "%zu %zu\n", &edge_real_id, &pos_info_count);
            VERIFY(read_count == 2);
            //    INFO(  edge_real_id);
            for (size_t j = 0; j < pos_info_count; j++) {
                int start_pos, end_pos;
                int m_start_pos, m_end_pos;
                read_count = fscanf(file, "%[^\n]s", cur_str);
                read_count = fscanf(file, "\n");
                read_count = sscanf(cur_str, "%s [%d - %d] --> [%d - %d]", contigId,
                                    &start_pos, &end_pos, &m_start_pos, &m_end_pos);
                //      INFO(cur_str);
                //      INFO (contigId<<" "<< start_pos<<" "<<end_pos);
                //      VERIFY(read_count == 3);
                VERIFY(read_count == 5);
                VERIFY(this->edge_id_map().find(edge_real_id) != this->edge_id_map().end());
                EdgeId eid = this->edge_id_map()[edge_real_id];
                edge_pos.AddEdgePosition(eid, string(contigId), start_pos - 1, end_pos, m_start_pos - 1, m_end_pos);
            }
        }
        fclose(file);
        return true;
    }

  private:
    Graph& g_;
    //  int edge_count_;
    map<size_t, EdgeId> edge_id_map_;
    map<size_t, VertexId> vertex_id_map_;

  protected:
    DataScanner(Graph &g) : g_(g) {
        INFO("Creating of scanner started");
        //    edge_count_ = 0;
    }

    Graph& g() {
        return g_;
    }

    map<size_t, EdgeId> &edge_id_map() {
        return edge_id_map_;
    }

    map<size_t, VertexId> &vertex_id_map() {
        return vertex_id_map_;
    }

    const map<size_t, EdgeId> &edge_id_map() const {
        return edge_id_map_;
    }

    const map<size_t, VertexId> &vertex_id_map() const {
        return vertex_id_map_;
    }

  public:
    virtual ~DataScanner() {

    }
};

template<class Graph>
class ConjugateDataScanner: public DataScanner<Graph> {
    typedef DataScanner<Graph> base;
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
private:
    restricted::IdSegmentStorage CreateIdStorage(const string& file_name) {
        FILE* file = fopen((file_name + ".gid").c_str(), "r");
        //This is to support compatibility to old saves. Will be removed soon
        if(file == NULL) {
            return this->g().GetGraphIdDistributor().ReserveUpTo(1000000000);
        }
        VERIFY_MSG(file != NULL, "Couldn't find file " << (file_name + ".gid"));
        size_t max;
        int flag = fscanf(file, "%zu\n", &max);
        VERIFY(flag == 1);
        fclose(file);
        return this->g().GetGraphIdDistributor().ReserveUpTo(max);
    }

  public:
    /*virtual*/
    void LoadGraph(const string& file_name) {
        restricted::IdSegmentStorage id_storage = CreateIdStorage(file_name);
        INFO("Trying to read conjugate de bruijn graph from " << file_name << ".grp");
        FILE* file = fopen((file_name + ".grp").c_str(), "r");
        VERIFY_MSG(file != NULL, "Couldn't find file " << (file_name + ".grp"));
        FILE* sequence_file = fopen((file_name + ".sqn").c_str(), "r");
        VERIFY_MSG(file != NULL, "Couldn't find file " << (file_name + ".sqn"));
        INFO("Reading conjugate de bruijn  graph from " << file_name << " started");
        size_t vertex_count;
        size_t edge_count;
        int flag = fscanf(file, "%zu %zu \n", &vertex_count, &edge_count);
        VERIFY(flag == 2);
        for (size_t i = 0; i < vertex_count; i++) {
            size_t vertex_real_id, conjugate_id;
            flag = fscanf(file, "Vertex %zu ~ %zu .\n", &vertex_real_id, &conjugate_id);
            TRACE("Vertex "<<vertex_real_id<<" ~ "<<conjugate_id<<" .");
            VERIFY(flag == 2);

            if (this->vertex_id_map().find((int) vertex_real_id) == this->vertex_id_map().end()) {
                size_t ids[2] = {vertex_real_id, conjugate_id};
                auto id_distributor = id_storage.GetSegmentIdDistributor(ids, ids + 2);
                VertexId vid = this->g().AddVertex(typename Graph::VertexData(), id_distributor);
                VertexId conj_vid = this->g().conjugate(vid);

                this->vertex_id_map()[vertex_real_id] = vid;
                this->vertex_id_map()[conjugate_id] = conj_vid;
            }
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

        const size_t longstring_size = 1000500; // TODO: O RLY magic constant? => Can't load edges >= 1Mbp
        char longstring[longstring_size];
        for (size_t i = 0; i < edge_count; i++) {
            size_t e_real_id, start_id, fin_id, length, conjugate_edge_id;
            flag = fscanf(file, "Edge %zu : %zu -> %zu, l = %zu ~ %zu .\n",
                          &e_real_id, &start_id, &fin_id, &length, &conjugate_edge_id);
            VERIFY(flag == 5);
            VERIFY(length < longstring_size);
            if (fasta) {
                flag = fscanf(sequence_file, ">%zu\n%s\n", &e_real_id, longstring);
            }
            else {
                flag = fscanf(sequence_file, "%zu %s .", &e_real_id, longstring);
            }
            VERIFY(flag == 2);
            TRACE("Edge " << e_real_id << " : " << start_id << " -> "
                  << fin_id << " l = " << length << " ~ " << conjugate_edge_id);
            if (this->edge_id_map().find((int) e_real_id) == this->edge_id_map().end()) {
                size_t ids[2] = {e_real_id, conjugate_edge_id};
                auto id_distributor = id_storage.GetSegmentIdDistributor(ids, ids + 2);
                Sequence tmp(longstring);
                EdgeId eid = this->g().AddEdge(this->vertex_id_map()[start_id], this->vertex_id_map()[fin_id], tmp, id_distributor);
                this->edge_id_map()[e_real_id] = eid;
                this->edge_id_map()[conjugate_edge_id] = this->g().conjugate(eid);
            }
        }
        fclose(file);
        fclose(sequence_file);
    }
  public:
    ConjugateDataScanner(Graph& g) :
            base(g) {
    }
};

inline std::string MakeSingleReadsFileName(const std::string& file_name,
                                    size_t index) {
    return file_name + "_paths_" + ToString(index) + ".mpr";
}

//helper methods
// todo think how to organize them in the most natural way

template<class Graph>
void PrintBasicGraph(const string& file_name, DataPrinter<Graph>& printer) {
    printer.SaveGraph(file_name);
    printer.SaveEdgeSequences(file_name);
    printer.SaveCoverage(file_name);
}

template<class graph_pack>
void PrintGraphPack(const string& file_name,
                    DataPrinter<typename graph_pack::graph_t>& printer,
                    const graph_pack& gp) {
    PrintBasicGraph(file_name, printer);
    //  printer.SavePaired(file_name + "_et", gp.etalon_paired_index);
    if (gp.edge_pos.IsAttached())
        printer.SavePositions(file_name, gp.edge_pos);
    if (gp.index.IsAttached())
        SaveEdgeIndex(file_name, gp.index.inner_index());
    if (gp.kmer_mapper.IsAttached())
        SaveKmerMapper(file_name, gp.kmer_mapper);
    if (gp.flanking_cov.IsAttached())
        printer.SaveFlankingCoverage(file_name, gp.flanking_cov);
}

template<class graph_pack>
void PrintGraphPack(const string& file_name, const graph_pack& gp) {
    ConjugateDataPrinter<typename graph_pack::graph_t> printer(gp.g);
    PrintGraphPack(file_name, printer, gp);
}

template<class Graph>
void PrintPairedIndex(const string& file_name, DataPrinter<Graph>& printer,
                      const PairedInfoIndexT<Graph>& paired_index) {
    printer.SavePaired(file_name, paired_index);
}

template<class Graph>
void PrintUnclusteredIndex(const string& file_name, DataPrinter<Graph>& printer,
                           const UnclusteredPairedInfoIndexT<Graph>& paired_index) {
    printer.SavePaired(file_name, paired_index);
}

template<class Graph>
void PrintClusteredIndex(const string& file_name, DataPrinter<Graph>& printer,
                         const PairedInfoIndexT<Graph>& clustered_index) {
    PrintPairedIndex(file_name + "_cl", printer, clustered_index);
}

template<class Graph>
void PrintScaffoldingIndex(const string& file_name, DataPrinter<Graph>& printer,
                         const PairedInfoIndexT<Graph>& clustered_index) {
    PrintPairedIndex(file_name + "_scf", printer, clustered_index);
}

template<class Graph>
void PrintScaffoldIndex(const string& file_name, DataPrinter<Graph>& printer,
    const PairedInfoIndexT<Graph>& scaffold_index) {
    PrintPairedIndex(file_name + "_scf", printer, scaffold_index);
}

template<class Graph>
void PrintUnclusteredIndices(const string& file_name, DataPrinter<Graph>& printer,
                             const UnclusteredPairedInfoIndicesT<Graph>& paired_indices) {
    for (size_t i = 0; i < paired_indices.size(); ++i)
        PrintUnclusteredIndex(file_name + "_" + ToString(i), printer, paired_indices[i]);
}

template<class Graph>
void PrintClusteredIndices(const string& file_name, DataPrinter<Graph>& printer,
                           const PairedInfoIndicesT<Graph>& paired_indices) {
    for (size_t i = 0; i < paired_indices.size(); ++i)
        PrintClusteredIndex(file_name  + "_" + ToString(i), printer, paired_indices[i]);
}

template<class Graph>
void PrintScaffoldingIndices(const string& file_name, DataPrinter<Graph>& printer,
                           const PairedInfoIndicesT<Graph>& paired_indices) {
    for (size_t i = 0; i < paired_indices.size(); ++i)
        PrintScaffoldingIndex(file_name  + "_" + ToString(i), printer, paired_indices[i]);
}

template<class graph_pack>
void PrintWithPairedIndex(const string& file_name,
                          DataPrinter<typename graph_pack::graph_t>& printer,
                          const graph_pack& gp,
                          const PairedInfoIndexT<typename graph_pack::graph_t>& paired_index,
                          bool clustered_index = false) {

    PrintGraphPack(file_name, printer, gp);
    if (!clustered_index) {
        PrintPairedIndex(file_name, printer, paired_index);
    } else {
        PrintClusteredIndex(file_name, printer, paired_index);
    }
}

template<class graph_pack>
void PrintWithClusteredIndex(const string& file_name,
                             DataPrinter<typename graph_pack::graph_t>& printer,
                             const graph_pack& gp,
                             const PairedInfoIndexT<typename graph_pack::graph_t>& paired_index) {
    PrintWithPairedIndex(file_name, printer, gp, paired_index, true);
}

template<class graph_pack>
void PrintWithPairedIndices(const string& file_name,
                            DataPrinter<typename graph_pack::graph_t>& printer,
                            const graph_pack& gp,
                            const PairedInfoIndicesT<typename graph_pack::graph_t>& paired_indices,
                            bool clustered_index = false) {
    PrintGraphPack(file_name, printer, gp);
    if (!clustered_index)
        PrintPairedIndices(file_name, printer, paired_indices);
    else
        PrintClusteredIndices(file_name, printer, paired_indices);
}

template<class graph_pack>
void PrintWithClusteredIndices(const string& file_name,
                               DataPrinter<typename graph_pack::graph_t>& printer,
                               const graph_pack& gp,
                               const PairedInfoIndicesT<typename graph_pack::graph_t>& paired_indices) {
    PrintWithPairedIndices(file_name, printer, gp, paired_indices, true);
}

template<class Graph>
void PrintSingleLongReads(const string& file_name, const LongReadContainer<Graph>& single_long_reads) {
    for (size_t i = 0; i < single_long_reads.size(); ++i){
        single_long_reads[i].DumpToFile(MakeSingleReadsFileName(file_name, i));
    }
}

template<class graph_pack>
void PrintAll(const string& file_name, const graph_pack& gp) {
    ConjugateDataPrinter<typename graph_pack::graph_t> printer(gp.g, gp.g.begin(), gp.g.end());
    PrintGraphPack(file_name, printer, gp);
    PrintUnclusteredIndices(file_name, printer, gp.paired_indices);
    PrintClusteredIndices(file_name, printer, gp.clustered_indices);
    PrintScaffoldingIndices(file_name, printer, gp.scaffolding_indices);
    PrintSingleLongReads(file_name, gp.single_long_reads);
    gp.ginfo.Save(file_name + ".ginfo");
}

template<class graph_pack, class VertexIt>
void PrintWithPairedIndex(const string& file_name, const graph_pack& gp,
                          VertexIt begin, VertexIt end,
                          const PairedInfoIndexT<typename graph_pack::graph_t>& paired_index,
                          bool clustered_index = false) {
    ConjugateDataPrinter<typename graph_pack::graph_t> printer(gp.g,
                                                                          begin, end);
    PrintWithPairedIndex(file_name, printer, gp, paired_index, clustered_index);
}

template<class graph_pack, class VertexIt>
void PrintWithClusteredIndex(const string& file_name, const graph_pack& gp,
                             VertexIt begin, VertexIt end,
                             const PairedInfoIndexT<typename graph_pack::graph_t>& clustered_index) {
    ConjugateDataPrinter<typename graph_pack::graph_t> printer(gp.g,
                                                                          begin, end);
    PrintWithPairedIndex(file_name, printer, gp, clustered_index, true);
}

template<class graph_pack>
void PrintWithPairedIndex(const string& file_name, const graph_pack& gp,
                          const PairedInfoIndexT<typename graph_pack::graph_t>& paired_index,
                          bool clustered_index = false) {
    PrintWithPairedIndex(file_name, gp, gp.g.begin(), gp.g.end(), paired_index,
                         clustered_index);
}

template<class graph_pack, class VertexIt>
void PrinGraphPack(const string& file_name, const graph_pack& gp,
                   VertexIt begin, VertexIt end) {
    ConjugateDataPrinter<typename graph_pack::graph_t> printer(gp.g,
                                                                          begin, end);
    PrintGraphPack(file_name, printer, gp);
}

template<class graph_pack>
void PrintWithClusteredIndex(const string& file_name, const graph_pack& gp,
                             const PairedInfoIndexT<typename graph_pack::graph_t>& clustered_index) {
    PrintWithPairedIndex(file_name, gp, clustered_index, true);
}

template<class graph_pack>
void PrintWithPairedIndices(const string& file_name, const graph_pack& gp,
                            const PairedInfoIndicesT<typename graph_pack::graph_t>& paired_indices,
                            bool clustered_index = false) {

    ConjugateDataPrinter<typename graph_pack::graph_t> printer(gp.g, gp.g.begin(), gp.g.end());

    PrintWithPairedIndices(file_name, printer, gp, paired_indices, clustered_index);
}

template<class graph_pack>
void PrintWithClusteredIndices(const string& file_name, const graph_pack& gp,
                               const PairedInfoIndicesT<typename graph_pack::graph_t>& paired_indices) {
    PrintWithPairedIndices(file_name, gp, paired_indices, true);
}

template<class Graph>
void ScanBasicGraph(const string& file_name, DataScanner<Graph>& scanner) {
    scanner.LoadGraph(file_name);
    scanner.LoadCoverage(file_name);
}

template<class graph_pack>
void ScanGraphPack(const string& file_name,
                   DataScanner<typename graph_pack::graph_t>& scanner, graph_pack& gp) {
    ScanBasicGraph(file_name, scanner);
    gp.index.Attach();
    if (LoadEdgeIndex(file_name, gp.index.inner_index())) {
        gp.index.Update();
    } else {
        WARN("Cannot load edge index, kmer coverages will be missed");
        gp.index.Refill();
    }
    //  scanner.LoadPaired(file_name + "_et", gp.etalon_paired_index);
    scanner.LoadPositions(file_name, gp.edge_pos);
    //load kmer_mapper only if needed
    if (gp.kmer_mapper.IsAttached())
        if (!LoadKmerMapper(file_name, gp.kmer_mapper)) {
            WARN("Cannot load kmer_mapper, information on projected kmers will be missed");
        }
    if (!scanner.LoadFlankingCoverage(file_name, gp.flanking_cov)) {
        WARN("Cannot load flanking coverage, flanking coverage will be recovered from index");
        gp.flanking_cov.Fill(gp.index.inner_index());
    }
}

template<class Graph>
void ScanPairedIndex(const string& file_name, DataScanner<Graph>& scanner,
                     UnclusteredPairedInfoIndexT<Graph>& paired_index,
                     bool force_exists = true) {
    scanner.LoadPaired(file_name, paired_index, force_exists);
}

template<class Graph>
void ScanClusteredIndex(const string& file_name, DataScanner<Graph>& scanner,
                        PairedInfoIndexT<Graph>& clustered_index,
                        bool force_exists = true) {
    scanner.LoadPaired(file_name + "_cl", clustered_index, force_exists);
}

template<class Graph>
void ScanScaffoldingIndex(const string& file_name, DataScanner<Graph>& scanner,
                          PairedInfoIndexT<Graph>& clustered_index,
                          bool force_exists = true) {
    scanner.LoadPaired(file_name + "_scf", clustered_index, force_exists);
}

template<class Graph>
void ScanPairedIndices(const std::string& file_name, DataScanner<Graph>& scanner,
                       UnclusteredPairedInfoIndicesT<Graph>& paired_indices,
                       bool force_exists = true) {
    for (size_t i = 0; i < paired_indices.size(); ++i)
        ScanPairedIndex(file_name  + "_" + ToString(i), scanner, paired_indices[i], force_exists);
}

template<class Graph>
void ScanClusteredIndices(const std:: string& file_name, DataScanner<Graph>& scanner,
                          PairedInfoIndicesT<Graph>& paired_indices,
                          bool force_exists = true) {
    for (size_t i = 0; i < paired_indices.size(); ++i)
        ScanClusteredIndex(file_name  + "_" + ToString(i), scanner, paired_indices[i], force_exists);
}

template<class Graph>
void ScanScaffoldingIndices(const std:: string& file_name, DataScanner<Graph>& scanner,
                            PairedInfoIndicesT<Graph>& paired_indices,
                            bool force_exists = true) {
    for (size_t i = 0; i < paired_indices.size(); ++i)
        ScanScaffoldingIndex(file_name  + "_" + ToString(i), scanner, paired_indices[i], force_exists);
}

template<class Graph>
void ScanScaffoldIndices(const string& file_name, DataScanner<Graph>& scanner,
        PairedInfoIndicesT<Graph>& scaffold_indices) {

    for (size_t i = 0; i < scaffold_indices.size(); ++i) {
        ScanScaffoldIndex(file_name  + "_" + ToString(i), scanner, scaffold_indices[i]);
    }
}

template<class graph_pack>
void ScanWithPairedIndex(const string& file_name,
                         DataScanner<typename graph_pack::graph_t>& scanner, graph_pack& gp,
                         PairedInfoIndexT<typename graph_pack::graph_t>& paired_index,
                         bool clustered_index = false) {
    ScanGraphPack(file_name, scanner, gp);
    if (!clustered_index) {
        ScanPairedIndex(file_name, scanner, paired_index);
    } else {
        ScanClusteredIndex(file_name, scanner, paired_index);
    }
}

template<class graph_pack>
void ScanWithPairedIndices(const string& file_name,
                           DataScanner<typename graph_pack::graph_t>& scanner, graph_pack& gp,
                           PairedInfoIndicesT<typename graph_pack::graph_t>& paired_indices,
                           bool clustered_index = false) {

    ScanGraphPack(file_name, scanner, gp);
    if (!clustered_index) {
        ScanPairedIndices(file_name, scanner, paired_indices);
    } else {
        ScanClusteredIndices(file_name, scanner, paired_indices);
    }
}

template<class graph_pack>
void ScanWithPairedIndex(const string& file_name, graph_pack& gp,
                         PairedInfoIndexT<typename graph_pack::graph_t>& paired_index,
                         bool clustered_index = false) {
    ConjugateDataScanner<typename graph_pack::graph_t> scanner(gp.g);
    ScanWithPairedIndex(file_name, scanner, gp, paired_index, clustered_index);
}

template<class graph_pack>
void ScanWithClusteredIndex(const string& file_name, graph_pack& gp,
                            PairedInfoIndexT<typename graph_pack::graph_t>& clustered_index) {
    ScanWithPairedIndex(file_name, gp, clustered_index, true);
}

template<class graph_pack>
void ScanWithClusteredIndices(const string& file_name,
                              DataScanner<typename graph_pack::graph_t>& scanner, graph_pack& gp,
                              PairedInfoIndicesT<typename graph_pack::graph_t>& paired_indices) {
    ScanWithPairedIndices(file_name, scanner, gp, paired_indices, true);
}

template<class graph_pack>
void ScanWithPairedIndices(const string& file_name, graph_pack& gp,
                           PairedInfoIndicesT<typename graph_pack::graph_t>& paired_indices,
                           bool clustered_index = false) {
    ConjugateDataScanner<typename graph_pack::graph_t> scanner(gp.g);
    ScanWithPairedIndices(file_name, scanner, gp, paired_indices, clustered_index);
}


template<class graph_pack>
void ScanWithClusteredIndices(const string& file_name, graph_pack& gp,
                              PairedInfoIndicesT<typename graph_pack::graph_t>& paired_indices) {
    ConjugateDataScanner<typename graph_pack::graph_t> scanner(gp.g);
    ScanGraphPack(file_name, scanner, gp);
    ScanClusteredIndices(file_name, scanner, paired_indices, false);
}

template<class Graph>
void ScanBasicGraph(const string& file_name, Graph& g) {
    ConjugateDataScanner<Graph> scanner(g);
    ScanBasicGraph<Graph>(file_name, scanner);
}

template<class Graph>
void ScanSingleLongReads(const string& file_name, LongReadContainer<Graph>& single_long_reads) {
    for (size_t i = 0; i < single_long_reads.size(); ++i){
        single_long_reads[i].LoadFromFile(MakeSingleReadsFileName(file_name, i), false);
    }
}

template<class graph_pack>
void ScanGraphPack(const string& file_name, graph_pack& gp) {
    ConjugateDataScanner<typename graph_pack::graph_t> scanner(gp.g);
    ScanGraphPack(file_name, scanner, gp);
}

template<class graph_pack>
void ScanAll(const std::string& file_name, graph_pack& gp,
             bool force_exists = true) {
    ConjugateDataScanner<typename graph_pack::graph_t> scanner(gp.g);
    ScanGraphPack(file_name, scanner, gp);
    ScanPairedIndices(file_name, scanner, gp.paired_indices, force_exists);
    ScanClusteredIndices(file_name, scanner, gp.clustered_indices, force_exists);
    ScanScaffoldingIndices(file_name, scanner, gp.scaffolding_indices, force_exists);
    ScanSingleLongReads(file_name,  gp.single_long_reads);
    gp.ginfo.Load(file_name + ".ginfo");
}
}
}
