//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <io/osequencestream.hpp>
#include "standard.hpp"
#include "graph_pack.hpp"
#include "config_struct.hpp"

#ifndef PROJECT_BWA_PAIR_INFO_FILLER_HPP_H
#define PROJECT_BWA_PAIR_INFO_FILLER_HPP_H

namespace bwa_pair_info {

class BWAPairInfoFiller {

private:
    debruijn_graph::conj_graph_pack& gp_;

    string bwa_path_;

    string work_dir_;

    size_t min_contig_len_;

    size_t nthreads_;

    string index_base_;

    bool index_constructed_;

    void OutputEdges(const string& filename) const {
        io::osequencestream_simple oss(filename);
        for (auto it = gp_.g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            debruijn_graph::EdgeId e = *it;
            if (gp_.g.length(e) >= min_contig_len_) {
                oss.set_header(ToString(gp_.g.int_id(e)));
                oss << gp_.g.EdgeNucls(e);
            }
        }
    }

    bool CreateIndex(const string& contigs) {
        int run_res = 0;
        string err_log = path::append_path(work_dir_, "index.err");
        string index_line = bwa_path_ + string(" index ") + "-a " + "is " + contigs + " 2>" + err_log;
        INFO("Running bwa index ...: " << index_line);
        run_res = system(index_line.c_str());
        if (run_res != 0) {
            ERROR("bwa index failed, cannot align reads");
            return false;
        }
        return true;
    }

    bool CreateIndexIfNeeded(const string& contigs) {
        bool result = false;
        if (index_constructed_ || CreateIndex(contigs)) {
            index_constructed_ = true;
            index_base_ = contigs;
            result = true;
        }
        return result;
    }


    bool RunBWA(const string& reads_file, const string& out_sam_file) const {
        string last_line = bwa_path_ + " mem " + index_base_ + " "  + reads_file + "  > " + out_sam_file + " 2>"
            + out_sam_file + ".txt";
        INFO("Running bwa mem ...:" << last_line);
        int run_res = system(last_line.c_str());
        if (run_res != 0) {
            ERROR("bwa index failed, cannot align reads");
            return false;
        }
        return true;
    }

    bool AlignLib(const io::SequencingLibrary<debruijn_graph::debruijn_config::DataSetData> &lib,
                  const string& sam_file_base,
                  vector<pair<string, string>>& resulting_sam_files) {

        VERIFY(index_constructed_);
        resulting_sam_files.clear();
        size_t file_index = 0;
        bool any_aligned = false;

        for (auto iter = lib.paired_begin(); iter != lib.paired_end(); iter++) {
            string left_reads = iter->first;
            string left_sam = sam_file_base + "_1_" + ToString(file_index) + ".sam";
            bool res = RunBWA(left_reads, left_sam);
            if (!res) {
                WARN("Failed to align left reads " << left_reads);
                continue;
            }
            string right_reads = iter->second;
            string right_sam = sam_file_base + "_2_" + ToString(file_index) + ".sam";
            res = RunBWA(right_reads, right_sam);
            if (!res) {
                WARN("Failed to align right reads " << right_reads);
                continue;
            }

            resulting_sam_files.push_back(make_pair(left_sam, right_sam));
            any_aligned = true;
        }
        return any_aligned;
    }

    bool ProcessPairedRead(const string& line) {

    }

    bool ParseSAMFile(const string& left_sam, const string& right_sam) {

    }


public:

    BWAPairInfoFiller(debruijn_graph::conj_graph_pack& gp,
                      const string& bwa_path,
                      const string& work_dir,
                      size_t min_contig_len = 1000,
                      size_t nthreads = 1):
        gp_(gp), bwa_path_(bwa_path), work_dir_(work_dir),
        min_contig_len_(min_contig_len), nthreads_(nthreads),
        index_base_(""), index_constructed_(false) {
    }

    bool Init() {
        path::make_dir(work_dir_);
        string edges_file = path::append_path(work_dir_, "long_edges.fasta");
        OutputEdges(edges_file);
        return CreateIndexIfNeeded(edges_file);
    }

    ~BWAPairInfoFiller() {
        //path::remove_if_exists(work_dir_);
    }

    bool ProcessLib(size_t lib_index) {
        string lib_dir =  path::append_path(work_dir_, ToString(lib_index));
        path::make_dir(lib_dir);
        vector<pair<string, string>> sam_files;

        if (!AlignLib(cfg::get().ds.reads[lib_index], path::append_path(lib_dir, "single"), sam_files)) {
            WARN("Failed to align lib #" << lib_index);
            return false;
        }



    }
};

};

#endif //PROJECT_BWA_PAIR_INFO_FILLER_HPP_H
