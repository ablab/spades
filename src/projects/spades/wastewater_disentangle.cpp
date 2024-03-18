//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "wastewater_disentangle.hpp"
#include "io/reads/file_reader.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "assembly_graph/core/graph_iterators.hpp"
#include "assembly_graph/core/graph.hpp"
#include "common/sequence/range.hpp"
#include "ssw/ssw_cpp.h"

#include <tuple>
#include <numeric>

struct Variation {
    size_t position_;
    char new_char_;
    size_t len_;
    std::string type_;
    double coverage_;
    Variation(size_t position, char new_char, size_t len, std::string type, double coverage)
    : position_(position), new_char_(new_char), len_(len), type_(type), coverage_(coverage) {

    }

    bool operator <(const Variation& x) const {
        if (this->position_ == x.position_) {
            return this->new_char_ < x.new_char_;
        }
        return this->position_ < x.position_;
    }

    bool operator ==(const Variation& x) {
        return (this->position_ == x.position_) && (this->type_ == x.type_) && (this->new_char_ == x.new_char_);
    }
};

struct Stats {
    std::string linaige_;
    unsigned match_;
    unsigned mismatch_;
    std::vector<Variation> found_variations_;

    Stats(const std::string &linaige, unsigned match, unsigned mismatch, const std::vector<Variation> &found_variations)
            : linaige_(linaige), match_(match), mismatch_(mismatch), found_variations_(found_variations) {

    }
    std::string ToString() const {
        return linaige_ + ": Matches - " + std::to_string(match_) + ", Mismatches - " + std::to_string(mismatch_);
    }
};


namespace debruijn_graph {
    // for string delimiter
    std::vector<std::string> split(std::string s, std::string delimiter) {
        size_t pos_start = 0, pos_end, delim_len = delimiter.length();
        std::string token;
        std::vector<std::string> res;

        while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
            token = s.substr (pos_start, pos_end - pos_start);
            pos_start = pos_end + delim_len;
            res.push_back (token);
        }

        res.push_back (s.substr (pos_start));
        return res;
    }

    Variation MakeVariation(const std::string &var_string) {
        return Variation(std::stoi(var_string.substr(1, var_string.length() - 2)), var_string[var_string.length() - 1], 1, "Mismatch", 0.0);
    }

    bool IsAligned(const std::vector<Range> &imperfect_alignment_union, size_t position) {
        for (auto r : imperfect_alignment_union) {
            if (r.start_pos <= position && r.end_pos >= position) {
                return true;
            }
        }
        return false;
    }


    std::map<std::string, std::vector<Variation>> ReadMatrix(const std::string &path_to_file, const std::vector<Range> &imperfect_alignment_union) {
        INFO(path_to_file);
        std::ifstream in(path_to_file);
        std::string delimeter = ",";
        std::string header;
        std::map<std::string, std::vector<Variation>> result;
        in >> header;
        auto header_vector = split(header, ",");
        std::map<size_t, std::string> id_to_var;
        for(size_t i = 0; i < header_vector.size(); ++i) {
            id_to_var[i] = header_vector[i];
        }

        while (!in.eof()) {
            std::string line;
            getline(in, line);
            auto line_vector = split(line, ",");
            std::string linaige = line_vector[0];
            for (size_t i = 1; i < line_vector.size(); ++i) {
                if (line_vector[i] == "1.0" && IsAligned(imperfect_alignment_union, std::stoi(header_vector[i].substr(1, header_vector[i].length() - 2)))) {
                    result[linaige].push_back(MakeVariation(header_vector[i]));
                }
            }
        }

        return result;
    }


    inline std::vector<Range> MergeIntervals(std::vector<Range> intervals) {
        if (intervals.size() == 0)
            return std::vector<Range>();
        std::sort(intervals.begin(), intervals.end());
        std::vector<Range> merged_intervals;
        merged_intervals.push_back(intervals[0]);
        for (size_t i = 1; i < intervals.size(); ++i) {
            if (!intervals[i].Intersect(merged_intervals.back())) {
                merged_intervals.push_back(intervals[i]);
            } else {
                Range new_range = Range(merged_intervals.back().start_pos, std::max(merged_intervals.back().end_pos, intervals[i].end_pos));
                merged_intervals.pop_back();
                merged_intervals.push_back(new_range);
            }
        }
        return merged_intervals;
    }

    std::vector<Stats> EvaluateSNPs(const std::vector<Variation> &variation_vector, const std::map<std::string, std::vector<Variation>> &matrix) {
        std::vector<Stats> stats_vector;
        for (auto p : matrix) {
            int match = 0;
            int mismatch = 0;
            std::vector<Variation> found_variations;
            for (auto &var : p.second) {
                bool found = false;
                for (auto &var2 : variation_vector) {
                    if (var == var2) {
                        var.coverage_ = var2.coverage_;
                        found = true;
                    }
                }
                if (found) {
                    match++;
                    found_variations.push_back(var);
                }
                else
                    mismatch++;
            }
            stats_vector.push_back(Stats(p.first, match, mismatch, found_variations));
        }
        return stats_vector;
    }


    void WastewaterDisentangle::SelectLinaiges(const std::vector<Stats> &stats, std::vector<Stats> &curated_linaiges, std::set<Variation> &unused_variations) {
        for (auto stat : stats) {
            if (stat.match_ < stat.mismatch_)
                continue;
            unsigned new_variations = 0;
            std::set<Variation> temp_variations;
            for (auto var : stat.found_variations_) {
                if (!unused_variations.count(var)) {
                    temp_variations.insert(var);
                    new_variations++;
                }
            }

            if (new_variations > 5) {
                curated_linaiges.push_back(stat);
                for (auto var : temp_variations)
                    unused_variations.insert(var);
            }
        }
    }

    std::map<std::string, double> WastewaterDisentangle::AssignRelativeCoverages(const std::vector<Stats> &curated_linaiges, const std::vector<double> &alpha_coverages,
                                                                                 const std::vector<std::map<char, double>> &meaningful_coverages) {
        unsigned total_strains = (unsigned)curated_linaiges.size();
        std::vector<std::map<char, double>> meaningful_ratios;
        meaningful_ratios.resize(meaningful_coverages.size());
        for (size_t i = 0; i < meaningful_coverages.size(); ++i) {
            double sum = 0.0;
            for (auto p : meaningful_coverages[i]) {
                sum += p.second;
            }
            if (math::eq(sum, 0.0))
                continue;

            for (auto p : meaningful_coverages[i]) {
                meaningful_ratios[i][p.first] = p.second / sum;
            }

        }

        double alpha_coverage = 0.0;
        if (alpha_coverages.size()) {
            alpha_coverage = std::accumulate(alpha_coverages.begin(), alpha_coverages.end(), 0.0) / float(alpha_coverages.size());
            total_strains++;
        }
        std::map<std::string, double> result;
        std::map<std::string, std::vector<double>> coverages;
        std::map<std::string, double> avg_coverages;
        for (auto linaige : curated_linaiges) {
            for (auto var : linaige.found_variations_) {
                bool found = false;
                for (auto linaige2 : curated_linaiges) {
                    if (found)
                        break;
                    if (linaige.linaige_ == linaige2.linaige_)
                        continue;
                    for (auto var2 : linaige2.found_variations_) {
                        if (var == var2) {
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    if (!(total_strains > 1 && math::eq(meaningful_ratios[var.position_][var.new_char_], 1.0))) {
                        INFO(meaningful_ratios[var.position_][var.new_char_]);
                        coverages[linaige.linaige_].push_back(meaningful_ratios[var.position_][var.new_char_]);
                    }
                }
            }
        }

        if (math::gr(alpha_coverage, 0.0)) {
            coverages["A"].push_back(meaningful_ratios[3037]['C']);
        }

        for (auto p : coverages) {
            double total = 0.0;
            for (auto cov : p.second)
                total += cov;
            avg_coverages[p.first] = total / (double)p.second.size();
        }

        double total_coverage = 0.0;
        for (auto p : avg_coverages)
            total_coverage += p.second;

        for (auto p : avg_coverages)
            result[p.first] = p.second/total_coverage * 100;

        return result;
    }

    void WastewaterDisentangle::run(graph_pack::GraphPack &gp, const char*) {
        const auto &graph = gp.get<Graph>();
        std::vector<size_t> trusted_contigs;
        for (size_t lib_id = 0; lib_id < cfg::get().ds.reads.lib_count(); ++lib_id) {
            if (cfg::get().ds.reads[lib_id].type() == io::LibraryType::TrustedContigs)
                trusted_contigs.push_back(lib_id);
        }
        io::SingleRead reference;
        for (auto lib_id : trusted_contigs) {
            auto stream = io::single_easy_reader(cfg::get().ds.reads[lib_id], false, false);
            while (!stream.eof()) {
                stream >> reference;
                break;
            }
        }
        auto ref_string = reference.GetSequenceString();
        if (ref_string.empty()) {
            INFO("Reference was not provided");
            return;
        }
        StripedSmithWaterman::Aligner aligner(  2, 5, 10, 5);
        StripedSmithWaterman::Filter filter;
        std::shared_ptr<StripedSmithWaterman::Alignment> alignment1(new StripedSmithWaterman::Alignment());
        std::shared_ptr<StripedSmithWaterman::Alignment> alignment2(new StripedSmithWaterman::Alignment());

        aligner.SetReferenceSequence(ref_string.c_str(), (int)ref_string.length());

        std::vector<Variation> variation_vector;
        std::vector<Range> perfect_alignments;
        std::vector<Range> imperfect_alignments;

        std::vector<double> alpha_coverages;
        std::vector<std::map<char, double>> meaningful_coverages;
        meaningful_coverages.resize(ref_string.size() + 1);

        for (const auto e1 : graph.edges()) {
            EdgeId e(e1);
            if (graph.coverage(e)< 20.0) {
                continue;
            }
            auto seq = graph.EdgeNucls(e);
            auto seq1_str = seq.str();
            aligner.Align(seq1_str.c_str(), filter, alignment1.get());
            auto seq2 = graph.EdgeNucls(graph.conjugate(e));
            auto seq2_str = seq2.str();
            aligner.Align(seq2_str.c_str(), filter, alignment2.get());
            auto alignment = alignment1->sw_score > alignment2->sw_score ? alignment1 : alignment2;
            auto seq_str = alignment1->sw_score > alignment2->sw_score ? seq1_str : seq2_str;

            size_t total = 0;
            size_t current = 0;
            if (alignment1->sw_score > alignment2->sw_score) {
                if (graph.OutgoingEdgeCount(graph.EdgeEnd(e)) == 0) {
                    total = seq1_str.length();
                } else {
                    total = seq1_str.length() - graph.k();
                }
            } else {
                if (graph.OutgoingEdgeCount(graph.EdgeEnd(graph.conjugate(e))) == 0) {
                    total = seq1_str.length();
                } else {
                    total = seq1_str.length() - graph.k();
                }
            }

            auto cigar = alignment->cigar;
            auto start_pos_query = alignment->query_begin;
            auto start_pos_reference = alignment->ref_begin;

            if(alignment->sw_score < 30)
                continue;

            auto first_type = (cigar.front() & 0x0F);
            auto last_type = (cigar.back() & 0x0F);
            auto first_len = (cigar.front() >> 4);
            auto last_len = (cigar.back() >> 4);

            if ((first_type == 4 && first_len > 10) || (last_type == 4 && last_len > 10)) // is S is too large we don't are likely to align incorrectly
                continue;


            DEBUG(alignment->cigar_string);
            if (first_type == 4) {
                for (size_t i = 0; i < first_len; ++i) {
                    if (first_len  > alignment->ref_begin + i)
                        continue;
                    meaningful_coverages[alignment->ref_begin + i - first_len + 1][seq_str[i]] += graph.coverage(e);
                    current++;
                    if (current == total)
                        continue;
                }
            }


            for (auto cigar_block : cigar) {
                if (current >= total)
                    break;
                auto len = (cigar_block >> 4);
                auto type = (cigar_block & 0x0F);
                DEBUG(cigar_block);
                DEBUG((cigar_block >> 4));
                DEBUG((cigar_block & 0x0F));
                DEBUG(len << type);

                switch (type) {
                    case 8:
                        for (size_t k = start_pos_query; k < start_pos_query + len; ++k) {
                            imperfect_alignments.push_back(Range(start_pos_reference, start_pos_reference + len - 1));
                            variation_vector.push_back(Variation(start_pos_reference + (k - start_pos_query) + 1, seq_str[k], 1, "Mismatch", graph.coverage(e)));
                            meaningful_coverages[start_pos_reference + (k - start_pos_query) + 1][seq_str[k]] += graph.coverage(e);
                            current++;
                            if (current == total)
                                break;
                        }
                        start_pos_reference += len;
                        start_pos_query += len;
                        break;
                    case 7:
                        perfect_alignments.push_back(Range(start_pos_reference, start_pos_reference + len - 1));
                        imperfect_alignments.push_back(Range(start_pos_reference, start_pos_reference + len - 1));
                        for (size_t k = start_pos_query; k < start_pos_query + len; ++k) {
                            meaningful_coverages[start_pos_reference + (k - start_pos_query) +
                                                 1][ref_string[start_pos_reference + (k - start_pos_query)]] += graph.coverage(e);
                            current++;
                            if (current == total)
                                break;

                        }
                        if (len > 140)
                            alpha_coverages.push_back(graph.coverage(e));
                        start_pos_reference += len;
                        start_pos_query += len;
                        break;
                    case 4: //S
                        break;
                    case 1: //I
                        variation_vector.push_back(Variation(start_pos_reference + 1, seq_str[start_pos_query], len, "Insertion", graph.coverage(e)));
                        start_pos_query += len;
                        current+= len;
                        if (current >= total)
                            break;
                        break;
                    case 2: //D
                        variation_vector.push_back(Variation(start_pos_reference + 1, ref_string[start_pos_reference], len, "Deletion", graph.coverage(e)));
                        imperfect_alignments.push_back(Range(start_pos_reference, start_pos_reference + len - 1));
                        start_pos_reference += len;
                        break;
                    default:
                        break;

                }

            }
            if (current >= total)
                continue;
            if (last_type == 4) {
                for (size_t i = 0; i < last_len; ++i) {
                    meaningful_coverages[alignment->ref_end + i + 1][seq_str[alignment->query_end + i]] += graph.coverage(e);
                    current++;
                    if (current == total)
                        continue;
                }
            }

        }

        std::vector<Range> alignment_union = MergeIntervals(perfect_alignments);

        std::vector<Range> imperfect_alignment_union = MergeIntervals(imperfect_alignments);

        size_t covered_bases_perfect = 0;
        size_t covered_bases_imperfect = 0;

        for (auto range : alignment_union) {
            covered_bases_perfect += range.size() + 1;
        }
        for (auto range : imperfect_alignment_union) {
            covered_bases_imperfect += range.size() + 1;
        }

        bool alpha_is_present = false;
        if (covered_bases_imperfect - covered_bases_perfect  < 4) {
             alpha_is_present = true;
        }

        std::sort(variation_vector.begin(), variation_vector.end());

        std::vector<Variation> new_variation_vector;
        new_variation_vector.push_back(variation_vector.front());
        for (size_t i = 1; i < variation_vector.size(); ++i) {
            if (variation_vector[i] == new_variation_vector.back()) {
                new_variation_vector.back().coverage_ += variation_vector[i].coverage_;
            } else {
                new_variation_vector.push_back(variation_vector[i]);
            }
        }
        std::swap(new_variation_vector, variation_vector);
        INFO("Printing variation vector");
        for (auto variation : variation_vector) {
            INFO(variation.type_ << " " << variation.position_ << " " << variation.len_);
        }
        auto matrix = ReadMatrix(cfg::get().sewage_matrix, imperfect_alignment_union);//"/home/dmm2017/Desktop/sewage/usher_barcodes.csv");
        INFO("Matrix parsing finished");
        auto stats_vector = EvaluateSNPs(variation_vector, matrix);
        INFO("Top results");
        std::sort(stats_vector.begin(), stats_vector.end(), [](const Stats &a, const Stats &b) -> bool {return a.match_ > b.match_;});
        std::vector<Stats> curated_linaiges;
        std::set<Variation> unused_variations;
        if (alpha_is_present) {
            DEBUG("Alpha variant is found in the sample");
        }
        else {
            DEBUG("Alpha variant is not found in the sample");
            alpha_coverages.clear();
        }

        SelectLinaiges(stats_vector, curated_linaiges, unused_variations);
        for (size_t i = 0; i < curated_linaiges.size(); ++i) {
            INFO("Top match: " << curated_linaiges[i].linaige_ << " " << curated_linaiges[i].match_ << " " << curated_linaiges[i].mismatch_);
        }
        auto coverages = AssignRelativeCoverages(curated_linaiges, alpha_coverages, meaningful_coverages);
        INFO("Estimated coverages:");
        std::ofstream out_stream(cfg::get().output_dir / "lineages.csv");
        out_stream << "Lineage,Abundance" << std::endl;
        for (auto p : coverages) {
            INFO("Linaige " << p.first << ": " << p.second << "%");
            out_stream << p.first << "," << p.second << std::endl;
        }

    }

}