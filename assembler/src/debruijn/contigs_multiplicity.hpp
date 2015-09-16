#pragma once

#include <array>
#include <string>
#include <iostream>
#include "path_extend/bidirectional_path.hpp"

using namespace path_extend;

class ContigsMultiplicity {
private:
    static const size_t COUNT_SAMPLES = 4;

    typedef typename std::array<uint32_t, COUNT_SAMPLES> MplStoringValue;
    typedef typename debruijn_graph::conj_graph_pack conj_graph_pack;
    typedef typename DefaultStoring::immutant_inverter<MplStoringValue> inverter;

public:
    static void Calculate(conj_graph_pack& graph_pack, const PathContainer& paths,
                          const std::string& kmers_mpl_file,
                          const std::string& contigs_mpl_file) {

        INFO("Create map");
        PerfectHashMap<conj_graph_pack::seq_t,
        MplStoringValue,
        kmer_index_traits<conj_graph_pack::seq_t>,
        DefaultStoring> kmers_mpl(graph_pack.index.k(),
                                   graph_pack.index.inner_index().workdir(),
                                   graph_pack.index.inner_index().index_ptr_);
        INFO("FillMpl");
        FillMplMap(kmers_mpl, graph_pack, kmers_mpl_file);
        INFO("CalculateContigs");
        CalculateContigsMpl(kmers_mpl, graph_pack, paths, contigs_mpl_file);
        INFO("Calculate end 31");
    }

private:
    template<typename T, size_t size>
    static void AddToArray(std::array<T, size>& res, const std::array<T, size>& from) {
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] += from[i];
        }
    };

    template<typename T>
    static void WriteBinary(std::ofstream& out, const T& value) {
        out.write((char*) &value, sizeof(T));
    }

    const static bool MY_INFO = false;

    template<typename K, typename T, typename S>
    static void FillMplMap(PerfectHashMap<K, MplStoringValue, T, S>& map, const conj_graph_pack& graph_pack,
                           const std::string& kmers_mpl_file) {
//        freopen("/home/toxa31/cool_log.txt", "w", stdout);
        //INFO("Fill start");
        for (auto it = map.value_begin(); it != map.value_end(); ++it) {
            it->fill(0);
        }
//        INFO("Fill zero end");
//        std::ifstream kmers_in(kmers_mpl_file + ".kmer", std::ios::binary);
//        std::ifstream kmers_mpl_in(kmers_mpl_file + ".mpl", std::ios::binary);
        std::ifstream kmers_in(kmers_mpl_file + ".kmer");
        std::ifstream kmers_mpl_in(kmers_mpl_file + ".mpl");
        if (MY_INFO) INFO("Files opened");
        size_t kmer_len = graph_pack.g.k() + 1;
        INFO("kmer_len = " << kmer_len);
        while (true) {
            string kmer_str;
            kmers_in >> kmer_str;
            if (kmers_in.fail()) {
                break;
            }
            if (MY_INFO) INFO("kmer create " << kmer_str);
            conj_graph_pack::seq_t kmer(kmer_len, kmer_str.c_str());
            MplStoringValue kmer_mpl;
            if (MY_INFO) INFO("kmer_mpl create");
//            kmer.BinRead(kmers_in);
            if (MY_INFO) INFO("kmer.BinRead");
            for (size_t i = 0; i < COUNT_SAMPLES; ++i) {
//                kmers_mpl_in.read((char*) &kmer_mpl[i], sizeof(uint32_t));
                kmers_mpl_in >> kmer_mpl[i];
                if (MY_INFO) INFO("kmer_mpl read");
            }
            if (kmers_in.fail() || kmers_mpl_in.fail()) {
                if (MY_INFO) INFO("one failed");
                break;
            }
            VERIFY(kmer_str.length() == kmer_len);
            if (MY_INFO) INFO("noone failed");
//            kmer = graph_pack.kmer_mapper.Substitute(kmer);
//            if (MY_INFO) INFO("substitude");
            if (!graph_pack.index.contains(kmer)) {
                if (MY_INFO) INFO("not contains");
                continue;
            }
            if (MY_INFO) INFO("contains");
            auto key_with_hash = map.ConstructKWH(kmer);
            if (MY_INFO) INFO("construct kwh");
            if (MY_INFO) INFO("construct kwh2");
            auto inverter_ = DefaultStoring::immutant_inverter<MplStoringValue>();
            if (MY_INFO) INFO("inverter");
            map.put_value(key_with_hash, kmer_mpl, inverter_);
            if (MY_INFO) INFO("put value");
        }
        INFO("Fill done");
//        fclose(stdout);
//        stdout = fdopen(1, "w");
    }

    template<typename K, typename T, typename S>
    static void CalculateContigsMpl(PerfectHashMap<K, MplStoringValue, T, S>& kmers_mpl, conj_graph_pack& graph_pack,
                                    const PathContainer& paths, const string& contigs_mpl_file) {
        std::ofstream output_id(contigs_mpl_file + ".id");
        std::ofstream output_mpl(contigs_mpl_file + ".mpl");

        size_t kmer_len = graph_pack.g.k() + 1;

        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            const BidirectionalPath& path = *iter.get();
            if (path.Length() <= 0) {
                continue;
            }
            MplStoringValue contig_cnt;
            contig_cnt.fill(0);
            for (size_t i = 0; i < path.Size(); ++i) {
                const Sequence& seq = graph_pack.g.EdgeNucls(path[i]);
                if (seq.size() < kmer_len) {
                    break;
                }
                auto kwh = kmers_mpl.ConstructKWH(conj_graph_pack::seq_t(kmer_len, seq));
                const std::array<uint32_t, COUNT_SAMPLES> &kmer_cnt = kmers_mpl.get_value(kwh, inverter());
                AddToArray(contig_cnt, kmer_cnt);

                for (size_t j = kmer_len; j < seq.size(); ++j) {
                    kwh <<= seq[j];
                    const MplStoringValue &kmer_cnt = kmers_mpl.get_value(kwh, inverter());
                    AddToArray(contig_cnt, kmer_cnt);
                }

                for (size_t i = 0; i < COUNT_SAMPLES; ++i) {
                    contig_cnt[i] /= seq.size() - kmer_len + 1;
                } 
            }
            output_id << path.GetId() << "\n";
            for (size_t i = 0; i < COUNT_SAMPLES; ++i) {
                if (i != 0) {
                    output_mpl << " ";
                }
                output_mpl << contig_cnt[i];
            }
            output_mpl << "\n";
        }
        output_id.close();
        output_mpl.close();
    }
};
